#!/usr/bin/env python
import argparse
import csv
import json
import logging
import os
import re
import sys

# Configure logging
_log = logging.getLogger("bdb")


ANISOU_PAT     = re.compile(r"^ANISOU")
"""
                e.g. 1czi, 1ent, 1smr, 4ape (entries refined with RESTRAIN)
REMARK   3  OTHER REFINEMENT REMARKS: THE QUANTITY GIVEN IN THE TEMPERATURE
REMARK   3  FACTOR FIELD OF THE *ATOM* AND *HETATM* RECORDS BELOW IS U**2,
REMARK   3  WHICH IS THE MEAN-SQUARE AMPLITUDE OF ATOMIC VIBRATION. THE
REMARK   3  TEMPERATURE FACTOR, B, CAN BE DERIVED BY THE FOLLOWING RELATION -
REMARK   3  B = 8 * (PI)**2 * U**2. IT IS AN INDICATION OF POSSIBLE ERRORS IN
REMARK   3  THE REFINEMENT THAT SOME ARE SLIGHTLY NEGATIVE.

                e.g. 1eed
REMARK   3  THE ATOMIC TEMPERATURE FACTORS IN THIS ENTRY ARE GIVEN AS
REMARK   3  U VALUES NOT B VALUES.  B VALUES MAY BE CALCULATED BY THE
REMARK   3  THE FOLLOWING: BISO (BISO = 8 PI==2== UISO).

                e.g. 5pep
REMARK   3  ISOTROPIC UISO VALUES ARE PROVIDED IN THE FIELD THAT
REMARK   3  USUALLY CONTAINS B VALUES.
"""
B_TYPE_PAT     = re.compile(r"^REMARK   3   B VALUE TYPE : (\w+( \w+)?)\s*$")
B_MSQAV_PAT    = re.compile(r"(MEAN-SQUARE AMPLITUDE OF ATOMIC VIBRATION|"
                             "U\*\*2|UISO)")
EXPDTA_PAT     = re.compile(r"^EXPDTA")
PDB_ID_PAT     = re.compile(r"^[0-9a-zA-Z]{4}$")
PROGRAM_PAT    = re.compile(r"^REMARK   3   PROGRAM     : "
                             "(?!NULL|NONE|NO REFINEMENT)")
REMARK_3_PAT   = re.compile(r"^REMARK   3")
REFMARKS_PAT   = re.compile(r"^REMARK   3  OTHER REFINEMENT REMARKS: "
                             "(?!NULL|NONE)")
REMARK_4_PAT   = re.compile(r"^REMARK   4")
FORMAT_PAT     = re.compile(r"^REMARK   4 [0-9A-Z]{4} COMPLIES WITH FORMAT V. "
                             "(?P<version>[\d.]+), "
                             "(?P<date>\d{2}-[A-Z]{3}-\d{2})")
N_TLS_PAT      = re.compile(r"^REMARK   3   NUMBER OF TLS GROUPS  : (\d+)\s*$")
TLS_REMARK_PAT = re.compile(r"^REMARK   3   "
                             "ATOM RECORD CONTAINS RESIDUAL B FACTORS ONLY")
RESIDUAL_PAT   = re.compile(r"RESIDUAL\s+([BU]-?\s*(FACTORS?|VALUES?)\s+)?"
                             "ONLY")
RESIDUAL_PAT2  = re.compile(r"ATOMIC\s+[BU]-?\s*(FACTORS?|VALUES?)\s+ARE"
                             "\s+RESIDUALS(\s+FROM\s+TLS(\s+REFINEMENT)?)?")
RESIDUAL_PAT3  = re.compile(r"[BU]-?\s*(FACTORS?|VALUES?)\s+ARE\s+RESIDUAL"
                             "\s+[BU]-?\s*(FACTORS?|VALUES?),?"
                             "(\(?WHICH\s+DO\s+NOT\s+"
                             "INCLUDE\s+THE\s+CONTRIBUTION\s+FROM\s+THE\s+TLS"
                             "\s+PARAMETERS\)?)?(.?\s+USE\s+TLSANL(\s+\(?\s*"
                             "DISTRIBUTED\s+WITH\s+CCP4\)?)?\s+TO\s+OBTAIN\s+"
                             "THE\s+(FULL\s+)?[BU]-?\s*(FACTORS?|VALUES?)"
                             ")?")
# eg. 2ix9:
RESIDUAL_PAT4  = re.compile(r"([BU]-?\s*(FACTORS?|VALUES?)\s+CORRESPOND\s+"
                             "TO\s+(THE\s+)?OVERALL?"
                             "[BU]-?\s*(FACTORS?|VALUES?)\s+EQUAL\s+TO\s+"
                             "THE\s+)?RESIDUAL[S]?\s+PLUS\s+(THE\s+)?TLS\s+"
                             "(COMPONENT)?")
SUM_TLS_RES    = re.compile(r"^REMARK   3   ATOM RECORD CONTAINS SUM OF TLS"
                             "AND RESIDUAL B FACTORS")
SUM_PAT        = re.compile(r"SUM\s+OF\s+TLS\s+AND\s+RESIDUAL\s+"
                             "[BU]-?\s*(FACTORS?|VALUES?)")
SUM_PAT2       = re.compile(r"[BU]-?\s*(FACTORS?|VALUES?)\s*:?\s+WITH\s+"
                             "TLS\s+ADDED")
# e.g. 2x2s 2wsp:
SUM_PAT3  = re.compile(r"(GLOBAL\s+)?[BU]-?\s*(FACTORS?|VALUES?),?\s*"
                             "(CONTAINING\s+)?RESIDUALS?\s+(AND|\+)\s+TLS\s+"
                             "COMPONENTS?(HAVE\s+BEEN\s+DEPOSITED)?")

def get_bdb_entry_outdir(root, pdb_id):
    out_dir = os.path.join(root, pdb_id[1:3], pdb_id)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    return out_dir

def get_b_value_type_from_string(s):
    """Return the B VALUE TYPE if it is present in this record.

    The returned B VALUE TYPE string can be "residual", "unverified" or None.
    """
    b_type = None
    #b_type = record.rstrip()[28:]
    m = re.search(B_TYPE_PAT, s)
    if m:
        b_type = m.group(1)
        if b_type == "LIKELY RESIDUAL":
            b_type = "residual"
        elif b_type == "UNVERIFIED":
            b_type = "unverified"
        else:
            _log.error("Unexpected B VALUE TYPE found: {1:s}", b_type)
    return b_type

def get_expdta_from_record(record):
    """Return the EXPDTA value if it is present in this record."""
    method = None
    if re.search(EXPDTA_PAT, record):
        method = record.rstrip()[10:]
    return method

def get_n_tls_groups(record):
    """Return the number of TLS groups if present in this record."""
    n_tls = None
    m = re.search(N_TLS_PAT, record)
    if m:
        n_tls = m.group(1)
        n_tls = None if n_tls == 0 else n_tls
        try:
            n_tls = int(n_tls)
        except exceptions.ValueError:
            _log.error("Unexpected value encountered for "
                       "NUMBER OF TLS GROUPS: {1:s}. None returned", n_tls)
            n_tls = None
    return n_tls

def get_pdb_header(pdb_file_path):
    """Return the PDB-file header records as a list.

    We assume the PDB file has the following composition:
    Header records
    [MODEL]
    ATOM
    (ANISOU)
    (SIGUIJ)
    (HETATM)
    TER
    [ENDMDL]
    Trailer records
    END
    """
    header = list()
    try:
        with open(pdb_file_path, "r") as pdb:
            for record in pdb:
                if re.search(r"^(MODEL|ATOM|HETATM)", record):
                    return header
                else:
                    header.append(record[0:79]) # keep trailing whitespace
    except IOError as ex:
        _log.error(ex)
    return header

def get_pdb_header_and_trailer(pdb_file_path):
    """Return the PDB-file header and trailer records as two lists.

    We assume the PDB file has the following composition:
    Header records
    [MODEL]
    ATOM
    (ANISOU)
    (SIGUIJ)
    (HETATM)
    TER
    [ENDMDL]
    Trailer records
    END
    """
    header = list()
    trailer = list()
    head_records = True
    try:
        with open(pdb_file_path, "r") as pdb:
            for record in pdb:
                if re.search(r"^(MODEL|ATOM|HETATM)", record):
                    head_records = False
                if head_records:
                    header.append(record[0:79]) # keep trailing whitespace
                elif not re.search(r"^(MODEL|ATOM|HETATM|ANISOU|SIGUIJ|"\
                                    "TER|ENDMDL|END\s+)", record):
                    trailer.append(record[0:79])
    except IOError as ex:
        _log.error(ex)
    return header, trailer

def get_pdb_trailer(pdb_file_path):
    """Return the PDB-file trailer records as a list.

    We assume the PDB file has the following composition:
    Header records
    [MODEL]
    ATOM
    (ANISOU)
    (SIGUIJ)
    (HETATM)
    TER
    [ENDMDL]
    Trailer records
    END
    """
    trailer = list()
    header = True
    try:
        with open(pdb_file_path, "r") as pdb:
            for record in pdb:
                if re.search(r"^(MODEL|ATOM|HETATM)", record):
                    header = False
                elif not header and not \
                        re.search(r"^(MODEL|ATOM|HETATM|ANISOU|SIGUIJ|"\
                                   "TER|ENDMDL|END\s+)", record):
                    trailer.append(record[0:79]) # keep trailing whitespace
    except IOError as ex:
        _log.error(ex)
    return trailer

def get_refmark_from_string(s):
    """Get the text written in REMARK 3, OTHER REFINEMENT REMARKS.

    Note: it is assumed the string is actually part
          of OTHER REFINEMENT REMARKS."""
    refmark = None
    # First line
    if re.search(REFMARKS_PAT, s):
        refmark = s.rstrip()[38:]
    elif re.search(REMARK_3_PAT, s):
        refmark = s.rstrip()[12:]
    return refmark

def get_refprog_from_string(s):
    """Return the refinement program if present in s."""
    program = None
    if re.search(PROGRAM_PAT, s):
        program = s.rstrip()[27:]
    return program

def get_raw_pdb_info(pdb_file_path):
    """Find the information necessary for bdb pipelin in a PDB file.

    Return a dictionary
    "expdta"       : the value of the EXPDTA record as a string.
    "b_msqav"      : True if the B-factor file contains U**2 (mean-square
                     amplitude of atomic vibration) instead of
                     8 * PI**2 * U**2 according to REMARK 3.
    "b_type"       : B VALUE TYPE. Can be "residual", "unverified" or None.
    "format_date"  : PDB format date.
    "format_vers"  : PDB format version.
    "has_anisou"   : True if ANISOU records are present in the PDB file.
    "other_refinement_remarks"
                   : OTHER REFINEMENT REMARKS in REMARK 3
                     as a single string.
    "refprog"      : the value of the refinement program record as a string.
    "tls_groups"   : number of TLS groups as integer. If not present: None.
    "tls_residual" : True if it is mentioned in the TLS details or elsewhere
                     that the ATOM records contain residual B-factors only.
    "tls_sum"      : True if it is mentioned somewhere in REMARK 3 that the
                     ATOM records contain the sum of TLS and residual B-factors
    """
    # TODO Note: what happens if programs are listed in multiple records?
    #            Try to parse from other part of header?
    b_type, expdta, format_date, format_vers, \
            other_refinement_remarks, refprog, tls_groups = (None,)*7
    b_msqav, has_anisou, refmarks, tls_residual, tls_sum = (False,)*5
    try:
        with open(pdb_file_path, "r") as pdb:
            for record in pdb:
                if re.search(EXPDTA_PAT, record):
                    expdta = get_expdta_from_record(record)
                elif re.search(REMARK_3_PAT, record):
                    if re.search(PROGRAM_PAT, record):
                        refprog = get_refprog_from_string(record)
                    elif re.search(B_TYPE_PAT, record):
                        b_type = get_b_value_type_from_string(record)
                    elif re.search(N_TLS_PAT, record):
                        tls_groups = get_n_tls_groups(record)
                    elif re.search(TLS_REMARK_PAT, record):
                        tls_residual = True
                    #elif re.search(RESIDUAL_PAT, record):
                    #    tls_residual = True
                    #elif re.search(RESIDUAL_PAT2, record):
                    #    tls_residual = True
                    #elif re.search(RESIDUAL_PAT3, record):
                    #    tls_residual = True
                    #elif re.search(RESIDUAL_PAT4, record):
                    #    tls_residual = True
                    elif re.search(SUM_TLS_RES, record):
                        tls_sum = True
                    #elif re.search(SUM_PAT, record):
                    #    tls_sum = True
                    #elif re.search(SUM_PAT2, record):
                    #    tls_sum = True
                    #elif re.search(SUM_PAT3, record):
                    #    tls_sum = True
                    # This should always be the final REMARK 3 remark...
                    if re.search(REFMARKS_PAT, record):
                        refmarks = True
                        other_refinement_remarks = \
                                get_refmark_from_string(record)
                    # ... so the next REMARK 3 lines can be appended
                    if refmarks and not re.search(REFMARKS_PAT, record):
                        other_refinement_remarks = other_refinement_remarks +\
                            " " + get_refmark_from_string(record)
                elif re.search(REMARK_4_PAT, record):
                    m = re.search(FORMAT_PAT, record)
                    if m:
                        format_date = m.group("date")
                        format_vers = m.group("version")
                        try:
                            format_vers = float(format_vers)
                        except ValueError:
                            _log.error("Unexpected value encountered for REMARK 4"\
                                    " FORMAT VERSION: {1:f}. None returned",
                                    format_vers)
                            format_vers = None
                if re.search(ANISOU_PAT, record):
                    has_anisou = True
    except IOError as ex:
        _log.error(ex)
    # Maybe the other refinement remarks contain something interesting
    if other_refinement_remarks:
        if re.search(RESIDUAL_PAT, other_refinement_remarks):
            tls_residual = True
        elif re.search(RESIDUAL_PAT2, other_refinement_remarks):
            tls_residual = True
        elif re.search(RESIDUAL_PAT3, other_refinement_remarks):
            tls_residual = True
        elif re.search(RESIDUAL_PAT4, other_refinement_remarks):
            tls_residual = True
        elif re.search(SUM_PAT, other_refinement_remarks):
            tls_sum = True
        elif re.search(SUM_PAT2, other_refinement_remarks):
            tls_sum = True
        elif re.search(SUM_PAT3, other_refinement_remarks):
            tls_sum = True
        elif re.search(B_MSQAV_PAT, other_refinement_remarks):
            b_msqav = True
    return {"b_msqav"                  : b_msqav,
            "b_type"                   : b_type,
            "expdta"                   : expdta,
            "has_anisou"               : has_anisou,
            "format_date"              : format_date,
            "format_vers"              : format_vers,
            "other_refinement_remarks" : other_refinement_remarks,
            "refprog"                  : refprog,
            "tls_groups"               : tls_groups,
            "tls_residual"             : tls_residual,
            "tls_sum"                  : tls_sum
            }

#class SingleLevelFilter(logging.Filter):
#    def __init__(self, passlevel, reject):
#        self.passlevel = passlevel
#        self.reject = reject
#
#    def filter(self, record):
#        if self.reject:
#            return (record.levelno != self.passlevel)
#        else:
#            return (record.levelno == self.passlevel)

def init_bdb_logger(pdb_id, root=".", global_log=False, append=False):
    """Configure logging for a single entry or globally."""
    logger = logging.getLogger("bdb")
    if not global_log:
        log_name = pdb_id + ".log"
        out_dir = get_bdb_entry_outdir(root, pdb_id)
        log_file = os.path.join(out_dir, log_name)
    else:
        log_file = os.path.join(root, "bdb.log")

    # File logger
    log_file = logging.FileHandler(log_file, "a" if append else "w")
    logger.addHandler(log_file)
    form1 = logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s")
    log_file.setFormatter(form1)

    # Console logger stdout
    stdout = logging.StreamHandler(sys.stdout)
    #f1 = SingleLevelFilter(logging.ERR, False)
    #stdout.addFilter(f1)
    form2 = logging.Formatter("%(message)s")
    stdout.setFormatter(form2)
    logger.addHandler(stdout)

    # Console logger stderr
    #stderr = logging.StreamHandler(sys.stderr)
    #f2 = SingleLevelFilter(logging.ERR, True)
    #stderr.addFilter(f2)
    #stderr.setFormatter(form1)
    #logger.addHandler(stderr)

    # Default log level
    logger.setLevel(logging.INFO)
    return logger

def is_valid_directory(parser, arg):
    """ Check if directory exists."""
    if not os.path.isdir(arg):
        parser.error("The directory {} does not exist!".format(arg))
    else:
        # File exists so return the directory
        return arg

def is_valid_file(parser, arg):
    """ Check if file exists and is not empty."""
    if not os.path.isfile(arg):
        parser.error("The file {} does not exist!".format(arg))
    elif not os.stat(arg).st_size > 0:
        parser.error("The file {} is empty!".format(arg))
    else:
        # File exists and is not empty so return the filename
        return arg

def is_valid_pdbid(parser, arg):
    """ Check if this is a valid PDB identifier (anno 2014)."""
    if not re.search(PDB_ID_PAT, arg):
        parser.error("Not a valid PDB ID: {} !".format(arg))
    else:
        return arg

def write_dict_json(d, filename="info.json", pretty=False):
    """Dump the dictionary to a JSON file."""
    try:
        with open(filename, "w") as f:
            if pretty:
                json.dump(d, f,
                          sort_keys=True, indent=4)
            else:
                json.dump(d, f)
    except IOError as ex:
        _log.error(ex)

def write_whynot(pdb_id, reason, filename=None, directory="."):
    """Create a WHY NOT file.

    Return a Boolean.
    """
    filename = pdb_id + ".whynot" if not filename else filename
    _log.warn(("{0:4s} | Writing WHY NOT entry.").format(pdb_id))
    try:
        with open(os.path.join(directory, filename), "w") as whynot:
            whynot.write("COMMENT: " + reason + "\n" +\
                         "BDB," + pdb_id + "\n")
            return True
    except IOError as ex:
        _log.error(ex)
        return False

if __name__ == "__main__":
    """Write WHY NOT entry."""
    parser = argparse.ArgumentParser(
            description="Create a WHY NOT entry or test BDB utils")
    parser.add_argument("-v", "--verbose", help="show verbose output",
                        action="store_true")
    subparsers = parser.add_subparsers(help="sub-command help")
    test = subparsers.add_parser("test", help="Test PDB file header parser")
    test.add_argument(
        "pdb_file_path",
        help="PDB file location.",
        type=lambda x: is_valid_file(parser, x))
    run = subparsers.add_parser("why_not", help="Create a WHY NOT entry")
    run.add_argument(
        "pdbid",
        help="PDB accession code.",
        type=lambda x: is_valid_pdbid(parser, x))
    run.add_argument(
        "output_dir",
        help="Directory where the output files should be saved.",
        type=lambda x: is_valid_directory(parser, x))
    run.add_argument(
        "--output_file",
        help="WHY NOT file name",
        default="whynot.txt")
    run.add_argument(
        "whynot",
        help="WHY NOT message.")
    args = parser.parse_args()
    if args.pdb_file_path:
        # Test mode
        _log = init_bdb_logger(args.pdb_file_path, global_log=True)
        if args.verbose:
            _log.setLevel(logging.DEBUG)
        d = get_raw_pdb_info(args.pdb_file_path)
        print json.dumps(d, sort_keys=True, indent=4)
    else:
        _log = init_bdb_logger(args.pdbin, global_log=True)
        if args.verbose:
            _log.setLevel(logging.DEBUG)
        if write_whynot(args.pdbid, args.whynot,
                    args.output_file, args.output_dir):
            sys.exit(0)
        else:
            sys.exit(1)
