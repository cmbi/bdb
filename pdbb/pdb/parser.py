#    BDB: A databank of PDB entries with full isotropic B-factors.
#    Copyright (C) 2014  Wouter G. Touw  (<wouter.touw@radboudumc.nl>)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License in the
#    LICENSE file that should have been included as part of this package.
#    If not, see <http://www.gnu.org/licenses/>.
import logging
_log = logging.getLogger(__name__)

import datetime
import os
import re


RE_BTYPE = re.compile(r"^  3   B VALUE TYPE : (?P<btype>.*)")
RE_REF_PROG = re.compile(r"^  3   PROGRAM     : (?P<refprogs>.*)")
RE_REF_REMARKS = re.compile(r"^  3  OTHER REFINEMENT REMARKS: (?!NULL|NONE)")
RE_B_MSQAV = re.compile(r"(MEAN-SQUARE AMPLITUDE OF ATOMIC VIBRATION|"
                        "U\*\*2|UISO)")
RE_FORMAT = re.compile(r"^  4 [0-9A-Z]{4} COMPLIES WITH FORMAT V. "
                       "(?P<version>[\d.]+), "
                       "(?P<date>\d{2}-[A-Z]{3}-\d{2})")
RE_TLS_GROUPS = re.compile(r"^  3   NUMBER OF TLS GROUPS  :\s*(\d+)\s*$")
RE_TLS_SEL = re.compile(r"""
        ^\s+3\s+RESIDUE\sRANGE\s:\s+
        (?P<ch_1>[0-9a-zA-Z])\s+
        (?P<rn_1>-?\d+)
        (?P<ic_1>[a-zA-Z]?)\s+
        (?P<ch_2>[0-9a-zA-Z])\s+
        (?P<rn_2>-?\d+)
        (?P<ic_2>[a-zA-Z]?)\s*$
        """, re.VERBOSE)
RE_TLS_RES = re.compile(r"^  3   ATOM RECORD CONTAINS RESIDUAL B FACTORS ONLY")
RE_TLS_RES_1 = re.compile(r"RESIDUAL\s+([BU]-?\s*(FACTORS?|VALUES?)\s+)?ONLY")
RE_TLS_RES_2 = re.compile(r"ATOMIC\s+[BU]-?\s*(FACTORS?|VALUES?)\s+(SHOWN\s+)?"
                          "ARE\s+RESIDUALS(\s+FROM\s+TLS(\s+REFINEMENT)?)?")
RE_TLS_RES_3 = re.compile(r"[BU]-?\s*(FACTORS?|VALUES?)\s+ARE\s+RESIDUAL"
                          "\s+[BU]-?\s*(FACTORS?|VALUES?),?"
                          "(\s+\(?WHICH\s+DO\s+NOT\s+"
                          "INCLUDE\s+THE\s+CONTRIBUTION\s+FROM\s+THE\s+TLS"
                          "\s+PARAMETERS\)?)?.?(\s+USE\s+TLSANL(\s+\(?\s*"
                          "DISTRIBUTED\s+WITH\s+CCP4\)?)?\s+TO\s+OBTAIN\s+"
                          "THE\s+(FULL\s+)?[BU]-?\s*(FACTORS?|VALUES?))?")
RE_TLS_SUM = re.compile(r"^  3   ATOM RECORD CONTAINS SUM OF TLS "
                        "AND RESIDUAL B FACTORS")
RE_TLS_SUM_1 = re.compile(r"SUM\s+OF\s+TLS\s+AND\s+RESIDUAL\s+"
                          "[BU]-?\s*(FACTORS?|VALUES?)")
RE_TLS_SUM_2 = re.compile(r"[BU]-?\s*(FACTORS?|VALUES?)\s*:?\s+WITH\s+"
                          "TLS\s+ADDED")
RE_TLS_SUM_3 = re.compile(r"(GLOBAL\s+)?[BU]-?\s*(FACTORS?|VALUES?),?\s*"
                          "(CONTAINING\s+)?RESIDUALS?\s+(AND|\+)\s+TLS\s+"
                          "COMPONENTS?(HAVE\s+BEEN\s+DEPOSITED)?")
RE_TLS_SUM_4 = re.compile(r"([BU]-?\s*(FACTORS?|VALUES?)\s+CORRESPOND\s+"
                          "TO\s+(THE\s+)?OVERALL?"
                          "[BU]-?\s*(FACTORS?|VALUES?)\s+EQUAL\s+TO\s+"
                          "THE\s+)?RESIDUAL[S]?\s+PLUS\s+(THE\s+)?TLS"
                          "(\s+COMPONENT)?")
RE_TLS_SUM_5 = re.compile(r"[BU]-?\s*(FACTORS?|VALUES?)\s*CONTAINS?\s+"
                          "(BOTH\s+)?TLS\s+AND\s+RESIDUALS?\s+COMPONENTS?.?")
RE_TLS_SUM_6 = re.compile(r"""
        (
            (ANISOTROPIC\s+)?
            [BU]-?\s*(FACTORS?|VALUES?)\s+
            (THAT\s+RESULT\s+FROM\s+)?
        )?
        (THE\s+)?
        COMBINATION\s+
        (OF\s+)?(THE\s+)?
        TLS\s+COMPONENTS?\s+
        (WITH\s+)?(THE\s+)?
        (RESIDUAL\s+)?
        (INDIVIDUAL\s+)?
        [BU]-?\s*(FACTORS?|VALUES?)\s*.?\s*
        """, re.VERBOSE)


def parse_pdb_file(pdb_file_path):
    """
    Parses the given pdb file, returning a dict where the key is the
    record name (e.g. 'ATOM   ') and the value is a list of all lines of that
    record name type.

    No validation is performed on the content of the pdb file.

    If the file at pdb_file_path doesn't exist, a ValueError is raised.
    """
    _log.info("Parsing pdb file {}".format(pdb_file_path))

    if not os.path.exists(pdb_file_path):
        _log.error("'{}' not found".format(pdb_file_path))
        raise ValueError("'{}' not found".format(pdb_file_path))

    with open(pdb_file_path) as pdb_file:
        records = {}
        for record in pdb_file:
            record_name = record[0:6]

            # If this is the first occurrence of a record name, initialise
            # the value with an empty list.
            if record_name not in records:
                records[record_name] = []
            records[record_name].append(record[7:])
        _log.debug("Parsed {0} records".format(len(records)))
        return records


def parse_dep_date(pdb_records):
    """
    Parses the deposition date from the pdb HEADER record, returning a
    date.

    The HEADER record is mandatory and spans a single line.

    If no HEADER records are found, a ValueError is raised.
    If multiple HEADER records are found, a ValueError is raised.
    If the date can not be parsed, a ValueError is raised.
    """
    _log.debug("Parsing deposition date from HEADER record")
    if "HEADER" not in pdb_records:
        _log.error("No HEADER record found")
        raise ValueError("No HEADER record found")

    if len(pdb_records["HEADER"]) > 1:
        _log.error("Multiple HEADER records found")
        raise ValueError("Multiple HEADER records found")

    dep_date = None
    record = pdb_records["HEADER"][0]
    try:
        dep_date = datetime.datetime.strptime(record[43:52], "%d-%b-%y")
    except ValueError as e:
        _log.error("{}".format(e))
        raise ValueError("Error parsing deposition date")
    return dep_date


def parse_exp_methods(pdb_records):
    """
    Parses the experiment methods from the pdb EXPDTA records, returning a
    list.

    The EXPDTA record is mandatory and can span multiple lines. If more than
    one experiment method is used, they are separated by a semi-colon.

    Although only a limited number of experiment methods are allowed, this
    function performs no validation of the values provided.

    If no EXPDTA records are found, a ValueError is raised.
    If no experimental method is found, an empty list is returned.
    """
    _log.info("Parsing experiment methods from EXPDTA records")
    if "EXPDTA" not in pdb_records:
        _log.error("No EXPDTA records found")
        raise ValueError("No EXPDTA records found")

    exp_methods = []
    for record in pdb_records["EXPDTA"]:
        for method in record.split(";"):
            exp_methods.append(method.strip())
    _log.debug("Found {} experiment methods".format(len(exp_methods)))
    return exp_methods


def parse_btype(pdb_records):
    """
    Parses the btype from the pdb REMARK records, returning a string.

    If no btype value is found, None is returned.
    """
    for record in pdb_records["REMARK"]:
        m = RE_BTYPE.search(record)
        if m is not None:
            btype = m.group("btype")
            btype = btype.rstrip()
            if btype == "LIKELY RESIDUAL":
                return "residual"
            elif btype == "UNVERIFIED":
                return "unverified"
            else:
                raise ValueError("Unexpected B VALUE TYPE found: {0:s}", btype)
    return None


def parse_other_ref_remarks(pdb_records):
    """
    Parses the other refinement remarks from the pdb REMARK records and returns
    a string of all entries concatenated, including newline characters.

    If no other refinement remarks are found, None is returned.
    """
    ref_rem = None
    for i, record in enumerate(pdb_records["REMARK"]):
        # If the regular expression matches, extract all of the following
        # REMARK 3 lines.
        if RE_REF_REMARKS.search(record):
            ref_rem = record[31:].rstrip()
            j = 1
            while (i+j < len(pdb_records["REMARK"]) and
                   pdb_records["REMARK"][i+j][0:3] == "  3"):
                ref_rem = ref_rem + " " + \
                    pdb_records["REMARK"][i+j][5:].rstrip()
                j = j + 1
            return ref_rem
    return None


def parse_format_date_version(pdb_records):
    """
    Parses the format date and version from the pdb REMARK records and returns
    the tuple (version, date) where date is a string and version is a float.

    If either the date or format are not found, None is returned for both.
    """
    for record in pdb_records["REMARK"]:
        m = RE_FORMAT.search(record)
        if m:
            format_date = m.group("date")
            format_vers = m.group("version")
            try:
                format_vers = float(format_vers)
            except (ValueError):
                _log.error("Unexpected value encountered for REMARK 4 FORMAT "
                           "VERSION: %s. None returned", format_vers)
                format_vers = None
            return format_vers, format_date
    return None, None


def parse_num_tls_groups(pdb_records):
    """
    Parses the number of tls groups from the pdb REMARK records, returning an
    integer of the number of groups.

    If the number of tls groups is not found, None is returned.
    If the number of tls groups is NULL, None is returned.
    """
    for record in pdb_records["REMARK"]:
        m = RE_TLS_GROUPS.search(record)
        if m is not None:
            n_tls = m.group(1)
            return int(n_tls)
    return None


def parse_tls_selection(pdb_records):
    """
    Parses the residue selection for a tls group, returning a
    dict.

    If the tls range selection is not found, an empty list is returned.
    """
    selections = []
    for record in pdb_records["REMARK"]:
        m = RE_TLS_SEL.search(record)
        if m is not None:
            chain_1 = m.group("ch_1")
            num_1 = m.group("rn_1")
            ic_1 = m.group("ic_1")
            chain_2 = m.group("ch_2")
            num_2 = m.group("rn_2")
            ic_2 = m.group("ic_2")
            selection = {"chain_1": chain_1,
                         "num_1": int(num_1),
                         "ic_1": None if ic_1 == '' else ic_1,
                         "chain_2": chain_2,
                         "num_2": int(num_2),
                         "ic_2": None if ic_2 == '' else ic_2}
            selections.append(selection)
    return selections


def parse_ref_prog(pdb_records):
    """
    Parses the refinement program from the pdb REMARK records, returning it as
    a string.

    If no refinement program is found, None is returned.
    """
    for record in pdb_records["REMARK"]:
        m = RE_REF_PROG.search(record)
        if m is not None:
            refprog = m.group("refprogs")
            refprog = refprog.rstrip()
            if (refprog == "NULL" or refprog == "NONE" or
                    refprog == "NO REFINEMENT"):
                return None
            return refprog
    return None


def is_bmsqav(other_refinement_remarks):
    """
    True if the B-factor file contains U**2 (mean-square amplitude of atomic
    vibration) instead of 8 * PI**2 * U**2 according to REMARK 3.
    """
    if other_refinement_remarks is None:
        return False

    if RE_B_MSQAV.search(other_refinement_remarks):
        return True
    return False


def is_tls_residual(pdb_records):
    """
    True if it is mentioned in the TLS details or elsewhere that the ATOM
    records contain residual B-factors only.

    First the REMARK 3 records are checked for conventional messages. If no
    evidence is found, the other refinement remarks are checked.
    """
    for record in pdb_records["REMARK"]:
        if RE_TLS_RES.search(record):
            return True

    other_refinement_remarks = parse_other_ref_remarks(pdb_records)
    if other_refinement_remarks is None:
        return False

    if RE_TLS_RES_1.search(other_refinement_remarks):
        return True
    if RE_TLS_RES_2.search(other_refinement_remarks):
        return True
    if RE_TLS_RES_3.search(other_refinement_remarks):
        return True
    return False


def is_tls_sum(pdb_records):
    """
    True if it is mentioned somewhere in REMARK 3 that the ATOM records contain
    the sum of TLS and residual B-factors

    First the REMARK 3 records are checked for conventional messages. If no
    evidence is found, the other refinement remarks are checked.
    """
    for record in pdb_records["REMARK"]:
        if RE_TLS_SUM.search(record):
            return True

    other_refinement_remarks = parse_other_ref_remarks(pdb_records)
    if other_refinement_remarks is None:
        return False

    if RE_TLS_SUM_1.search(other_refinement_remarks):
        return True
    if RE_TLS_SUM_2.search(other_refinement_remarks):
        return True
    if RE_TLS_SUM_3.search(other_refinement_remarks):
        return True
    if RE_TLS_SUM_4.search(other_refinement_remarks):
        return True
    if RE_TLS_SUM_5.search(other_refinement_remarks):
        return True
    if RE_TLS_SUM_6.search(other_refinement_remarks):
        return True
    return False


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
    header = []
    trailer = []
    head_records = True
    with open(pdb_file_path, "r") as pdb:
        for record in pdb:
            if re.search(r"^(MODEL|ATOM|HETATM)", record):
                head_records = False
            if head_records:
                header.append(record[0:80])  # keep trailing whitespace
            elif not re.search(r"^(MODEL|ATOM|HETATM|ANISOU|SIGUIJ|"
                               "TER|ENDMDL|END\s+)", record):
                trailer.append(record[0:80])
    return header, trailer
