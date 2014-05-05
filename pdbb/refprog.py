#!/usr/bin/env python
from __future__ import print_function

import argparse
import logging
import os
import re
import sys

from pdbb.pdb.parser import (parse_pdb_file, parse_other_ref_remarks, is_bmsqav,
                            parse_btype, parse_format_date_version,
                            parse_num_tls_groups, parse_ref_prog,
                            is_tls_residual, is_tls_sum)
from pdbb.bdb_utils import write_whynot
from pdbb.check_beq import check_beq, report_beq


_log = logging.getLogger(__name__)


RE_BEXCEPT = re.compile(r"""
                             \s+SUM\s+|
                             (
                                (?<!MAXIMUM LIKELIHOOD\s)
                                RESIDUALS?
                                (?!\s+(FEATURES?|ELECTRON|DENSITY))
                                \s+
                             )
                          """, re.VERBOSE)


def is_bdb_includable_refprog(refprog):
    """Check if the refinement program can be included in the bdb.

    Warning: always return True now, might be changed later after inspection.

    Return a Boolean.
    """
    # useful_pat = re.compile(r"""
    #     ^(
    #       REFMAC
    #     | CNS           # Anisotropic B-factors nor TLS groups implemented
    #     | X-PLOR        # Anisotropic B-factors nor TLS groups implemented
    #     | BUSTER
    #     )$
    # """, re.VERBOSE)
    # return re.search(useful_pat, refprog)
    return True


def decide_refprog_restrain(pdb_info):
    """Determine whether a BDB file can be created for this RESTRAIN PDB file.

    Only applicable for structures refined by RESTRAIN.

    The decision is based on the refinement program interpreted from the
    PDB file. Furthermore, several remarks and details in the
    header are used.

    WARNING: this code assumes the Beq values from ANISOU records are not
             identical to the reported corresponding B-factors

    Return a tuple of Booleans and the decision message:
    useful       : True if a BDB file can probably be created for this PDB file
    assume_iso   : whether the PDB file should be assumed to have
                   total isotropic B-factors
    req_tlsanl   : True if running TLSANL is required
    """

    # Check and declare

    assert pdb_info["prog_last"][0] == "RESTRAIN"

    (useful, assume_iso, req_tlsanl) = (False, False, False)
    msg = "{}: ".format(pdb_info["prog_last"][0])

    # If ANISOU records are present, make sure Beq values are not OK
    if pdb_info["has_anisou"]:
        assert pdb_info["beq_identical"] <= 0.9999

    RE_U = re.compile(r"""
                        THE\s+QUANTITY\s+PRESENTED\s+
                        IN\s+THE\s+TEMPERATURE\s+FACTOR\s+
                        FIELD\s+IS\s+U\.
                        """, re.VERBOSE)

    # Decide

    if pdb_info["b_msqav"]:
        useful = True
        msg = "{}B-factor field contains mean square atomic displacement "\
              "(REMARK 3)".format(msg)
    elif re.search(RE_U, pdb_info["other_refinement_remarks"]):
        """ e.g. 3cms, 4cms """
        msg = "{}B-factor field contains \"U\" (REMARK 3)".format(msg)
    elif pdb_info["b_type"] or pdb_info["has_anisou"] or \
            pdb_info["tls_groups"] or pdb_info["tls_residual"] or \
            pdb_info["tls_sum"]:
        msg = "{}unexpected content cannot (yet) be handled".format(msg)
    else:
        msg = "{}probably full B-factors".format(msg)
        useful = True
        assume_iso = True

    # Report

    if not useful:
        write_whynot(pdb_info["pdb_id"], msg)
        _log.warn("{}.".format(msg))
    else:
        _log.info("{}".format(msg))

    return useful, assume_iso, req_tlsanl, msg


def decide_refprog_remediation(pdb_info):
    """Determine whether a PDB has been flagged by a PDB remediation.

    The wwPDB has reviewed the type of ADPs during the 2011
    remediation:
    http://www.wwpdb.org/documentation/2011remediation_overview-061711.pdf

    However we also find this type of remark for earlier versions...
    earlier remediations?

    Return a tuple of three Booleans and the decision message.
    useful       : True if a BDB file can probably be created for this PDB file
    assume_iso   : whether the PDB file should be assumed to have
                   total isotropic B-factors
    req_tlsanl   : True if running TLSANL is required (also set for other
                   programs than REFMAC)
    """

    # Check and declare

    (useful, assume_iso, req_tlsanl) = (False, False, False)
    msg = "{}: ".format(pdb_info["prog_last"][0])

    # Decide

    if pdb_info["b_type"] == "residual":
        req_tlsanl = True
        msg = "{}residual B-factors (wwPDB remediation{})".format(msg,
                " 2011" if pdb_info["format_vers"] == 3.30 else "")
        if pdb_info["prog_last"][0] == "REFMAC":
            useful = True
    elif pdb_info["b_type"] == "unverified":
        msg = "{}B-factor type could not be determined (wwPDB remediation{})".\
        format(msg, " 2011" if pdb_info["format_vers"] == 3.30 else "")
    elif pdb_info["prog_last"][0] == "REFMAC" and \
            pdb_info["format_vers"] == 3.30:
        # The REFMAC entry was found to contain full B-factors
        useful = True
        assume_iso = True
        msg = "{}full B-factors (wwPDB remediation 2011)".format(msg)
    else:
        msg = "{}unexpected B-factor type annotation (wwPDB remediation)".\
            format(msg)

    # Report

    if not useful:
        write_whynot(pdb_info["pdb_id"], msg)
        _log.warn("{}.".format(msg))
    else:
        _log.info("{}".format(msg))

    return useful, assume_iso, req_tlsanl, msg


def decide_refprog_tls(pdb_info):
    """Determine whether a BDB entry can be created for this TLS-refined PDB
    file.

    Only applicable for PDB files with defined TLS groups.

    Return a tuple of three Booleans and the decision message.
    useful       : True if a BDB file can probably be created for this PDB file
    assume_iso   : whether the PDB file should be assumed to have
                   total isotropic B-factors
    req_tlsanl   : True if running TLSANL is required (also set for other
                   programs than REFMAC)
    """

    # Check and declare

    (useful, assume_iso, req_tlsanl) = (False, False, False)
    msg = "{}: ".format(pdb_info["prog_last"][0])
    bneq_msg = "Not enough B-factors could be reproduced from ANISOU records"

    # If ANISOU records are present, make sure Beq values are not OK
    if pdb_info["has_anisou"]:
        assert pdb_info["beq_identical"] <= 0.9999

    # Decide

    if pdb_info["tls_residual"]:  # Mentioned residual B-factors
        if pdb_info["has_anisou"]:
            msg = "{}TLS group(s), residual B-factors (REMARK 3) "\
                  "and ANISOU records. {}".format(msg, bneq_msg)
        else:
            msg = "{}TLS group(s), residual B-factors (REMARK 3) "\
                  "without ANISOU records".format(msg)
            if pdb_info["prog_last"][0] == "REFMAC":
                useful = True
            # else: # inspect first
                # useful = True
            req_tlsanl = True
    elif pdb_info["tls_sum"]:  # Mentioned full B-factors
        if pdb_info["has_anisou"]:
            # (probably, TLSANL was run)
            msg = "{}TLS group(s), full B-factors (REMARK 3) "\
                  "and ANISOU records. {}".format(msg, bneq_msg)
        else:
            """ REFMAC: e.g 2aen, 2wo7, 2wwj, 2wyb, 2wyc, 2wye, 2x00, 2x5s,
            2xgq, 2xgx, 2xhk, 2xko, 2xkp, 2xr9
            """
            msg = "{}TLS group(s), full B-factors (REMARK 3) "\
                  "without ANISOU records".format(msg)
            assume_iso = True
            if pdb_info["prog_last"][0] == "REFMAC":
                useful = True
            # else: # inspect first
                # useful = True
    elif re.search(RE_BEXCEPT, pdb_info["other_refinement_remarks"]) and \
            re.search("[BU]-?\s*(FACTORS?|VALUES?)",
                    pdb_info["other_refinement_remarks"]): # any exceptions?
        msg = "{}TLS group(s) and, possibly, residual or full B-factors "\
              "(REMARK 3, unrecognized format)".format(msg)
    else:  # TLS refinement without hints about B-value type
        if pdb_info["has_anisou"]:
            """ REFMAC: e.g. 1oj7, 1pm7, 2gcl, 3c2s, 3l6r, 3o1a """
            msg = "{}TLS group(s), no B-factor type details (REMARK 3). {}".\
                format(msg, bneq_msg)
        else:
            # This will be a large group, inspect first
            msg = "{}TLS group(s), no B-factor type details (REMARK 3) "\
                  "and no ANISOU records".format(msg)

    # Report

    if not useful:
        write_whynot(pdb_info["pdb_id"], msg)
        _log.warn("{}.".format(msg))
    else:
        _log.info("{}".format(msg))

    return useful, assume_iso, req_tlsanl, msg


def decide_refprog_notls(pdb_info):
    """Determine whether a BDB entry can be created for this PDB file without
    TLS groups.

    Only applicable for PDB files without defined TLS groups.

    Return a tuple of three Booleans and the decision message.
    useful       : True if a BDB file can probably be created for this PDB file
    assume_iso   : whether the PDB file should be assumed to have
                   total isotropic B-factors
    req_tlsanl   : True if running TLSANL is required (also set for other
                   programs than REFMAC)
    """

    # Check and declare

    (useful, assume_iso, req_tlsanl) = (False, False, False)
    msg = "{}: ".format(pdb_info["prog_last"][0])
    bneq_msg = "Not enough B-factors could be reproduced from ANISOU records"

    # If ANISOU records are present, make sure Beq values are not OK
    if pdb_info["has_anisou"]:
        assert pdb_info["beq_identical"] <= 0.9999


    # Decide

    if re.search(r"TLS", pdb_info["other_refinement_remarks"]):
        """ Any exceptions?
        REFMAC: e.g. 1fse, 1gmm, 1gqq, 1h3g, 1jnx, 1krh, 1muu, 1oc0, 1oiq,
        1oir, 1oit, 1ux9, 1uzl, 1zca, 2hwy, 2jbm, 2oeu, 2pe4, 2wnl, 2wvi,
        2zze, 2zzf, 2zzg, 3dem, 3hbt
        """
        msg = "{}TLS remark without TLS group(s)".format(msg)
        if pdb_info["has_anisou"]:
            msg = "{}. {}".format(msg, bneq_msg)

        if pdb_info["tls_sum"]:
            assume_iso = True
            msg = "{}. Full B-factors (REMARK 3)".format(msg)
            if pdb_info["prog_last"][0] == "REFMAC":
                useful = True
            # else: # inspect first
                # useful = True
        elif pdb_info["tls_residual"]:
            req_tlsanl = True
            msg = "{}. Residual B-factors (REMARK 3)".format(msg)
            if pdb_info["prog_last"][0] == "REFMAC":
                useful = True
    elif pdb_info["tls_residual"]:
        """ REFMAC: e.g. 3ch0, 2pq7, 3h3z """
        msg = "{}residual B-factors without TLS group(s) (REMARK 3)".\
            format(msg)
        if pdb_info["has_anisou"]:
            msg = "{}. {}".format(msg, bneq_msg)
    elif pdb_info["tls_sum"]:
        msg = "{}full B-factors without TLS group(s) (REMARK 3)".format(msg)
        if pdb_info["has_anisou"]:
            msg = "{}. {}".format(msg, bneq_msg)
    elif re.search(RE_BEXCEPT, pdb_info["other_refinement_remarks"]) and \
            re.search("[BU]-?\s*(FACTORS?|VALUES?)",
                    pdb_info["other_refinement_remarks"]): # any exceptions?
        msg = "{}possibly, residual or full B-factors (REMARK 3, unrecognized"\
              " format). No TLS groups".format(msg)
        if pdb_info["has_anisou"]:
            msg = "{}. {}".format(msg, bneq_msg)
    elif pdb_info["has_anisou"]:
        msg = "{}probably full/mixed anisotropic refinement. {}".format(
                msg, bneq_msg)
    else:
        # Probably refinement without any type of
        # anisotopic displacement parameters
        useful = True
        assume_iso = True
        msg = "{}probably full B-factors".format(msg)

    # Report

    if not useful:
        write_whynot(pdb_info["pdb_id"], msg)
        _log.warn("{}.".format(msg))
    else:
        _log.info("{}".format(msg))

    return useful, assume_iso, req_tlsanl, msg


def decide_refprog(pdb_info):
    """Determine whether refinement program can be used in the bdb project.

    The decision is based on the refinement program interpreted from the
    PDB file. Furthermore, several remarks and details in the
    header are used.

    WARNING: this code assumes the Beq values from ANISOU records are not
             identical to the reported corresponding B-factors

    Return a tuple of Booleans and the decision message:
    useful       : True if a BDB file can probably be created for this PDB file
    assume_iso   : whether the PDB file should be assumed to have
                   total isotropic B-factors
    req_tlsanl   : True if running TLSANL is required (also set for other
                   programs than REFMAC)
    """

    # Check and declare
    assert isinstance(pdb_info["prog_last"], list)
    (useful, assume_iso, req_tlsanl) = (False, False, False)

    # If ANISOU records are present, make sure Beq values are not OK
    if pdb_info["has_anisou"]:
        assert pdb_info["beq_identical"] <= 0.9999

    # There must be only a single refinement program. If not, return False for
    # all values in the return tuple.
    if len(pdb_info["prog_last"]) > 1:
        refprog = [str(p) for p in pdb_info["prog_last"]]
        msg = "Combination of refinement programs cannot (yet) be "\
                  "included in the bdb: {}".format(" and ".join(refprog))
        write_whynot(pdb_info["pdb_id"], msg)
        _log.warn("{}.".format(msg))
        return False, False, False, msg
    elif len(pdb_info["prog_last"]) == 0:
        """ e.g. 3cw1 """
        msg = "Program(s) in REMARK 3 not interpreted as refinement "\
                  "program(s)"
        write_whynot(pdb_info["pdb_id"], msg)
        _log.error("{}.".format(msg))
        return False, False, False, msg

    assert isinstance(pdb_info["prog_last"][0], str)
    _log.info("Interpreted last-used refinement program: {}.".format(
        pdb_info["prog_last"][0]))

    msg = "{}: ".format(pdb_info["prog_last"][0])

    # Check if the refinement program is supported. If not, return False for
    # all values in return tuple.
    if not is_bdb_includable_refprog(pdb_info["prog_last"][0]):
        msg = "{}this program cannot (yet) be included".format(msg)
        write_whynot(pdb_info["pdb_id"], msg)
        _log.warn("{}.".format(msg))
        return False, False, False, msg

    # The orr may have to be searched
    if pdb_info["other_refinement_remarks"] is None:
        pdb_info["other_refinement_remarks"] = ""

    # Start deciding

    # Special treatment for RESTRAIN files
    if pdb_info["prog_last"][0] == "RESTRAIN":
        return decide_refprog_restrain(pdb_info)

    # REFMAC and other refprogs:

    # We trust the wwPDB remediation decision
    if (pdb_info["format_vers"] == 3.30 and \
            pdb_info["prog_last"][0] == "REFMAC") \
        or pdb_info["b_type"]:
        return decide_refprog_remediation(pdb_info)

    # B-factors cannot be residual and full at the same time
    if pdb_info["tls_residual"] and pdb_info["tls_sum"]:
        msg = "{}residual and full B-factors (REMARK 3)".format(msg)
        write_whynot(pdb_info["pdb_id"], msg)
        _log.warn("{}.".format(msg))
        return False, False, False, msg

    # TLS groups defined in the PDB file..
    if pdb_info["tls_groups"]:
        return decide_refprog_tls(pdb_info)
    else: # ..or not (?)
        return decide_refprog_notls(pdb_info)

def except_refprog_warn():
    message = "Pre-defined exceptional refinement program case found"
    _log.warn("{}.".format(message))


def filter_progs(pin, pv):
    """Remove non-refinement programs.

    Return lists of program(s) and version(s)
    """
    assert isinstance(pin, list)
    assert isinstance(pv, list)
    assert len(pin) == len(pv)
    pop_progs = [
        "ARP",             # model building
        "ARP/WARP",        # model building
        "CHAIN",           # graphics
        "COOT",            # model building, graphics
        "DM",              # density modification
        "FRODO",           # graphics
        "HKL-3000",        # data collection
        "LAFIRE",          # uses CNS, Refmac5, phenix.refine or autoBUSTER
        "MOLPROBITY",      # validation
        "O",               # model building, graphics
        "OOPS",            # model building, graphics
        "PIKSOL",
        "PRODRG",          # small-mol topology
        "PROTEIN",         # crytallographic data analysis
        "QUANTA",          # model building, workbench
        "SCWRL",           # predict side-chains
        "SFALL",           # structure factor calculation
        "SOLVE/RESOLVE",   # phasing
        "TOM",             # model building, graphics
        "TOM/FRODO",       # model building, graphics
        "XFIT",            # peak fitting
        "XPLEO",           # real space model completion
        "XTALVIEW"         # workbench
        ]
    for rp in pin:
        if rp in pop_progs:
            pop = pin.index(rp)
            pin = pin[:pop] + pin[(pop + 1):]
            pv = pv[:pop] + pv[(pop + 1):]
    return pin, pv


def get_refi_data(pdb_records, structure, pdb_id):
    """Determine whether this PDB file can be used in the bdb project.

    The decision is based on refinement details parsed from the header.

    If entries have ANISOU records, Beq values are compared with
    the B-factor values in the ATOM records.
    Warning: it is assumed that B-factors are full and isotropic if they can be
    reproduced from Beq values.

    Return a dict containing refinement info.
    "assume_iso"   : whether the PDB file should be assumed to have
                     total isotropic B-factors
    "b_msqav"      : True if the B-factor file contains U**2 (mean-square
                     amplitude of atomic vibration) instead of
                     8 * PI**2 * U**2 according to REMARK 3.
    "b_type"       : B VALUE TYPE. Can be "residual", "unverified" or None.
    "beq_identical": Percentage of Beq values (from ANISOU records) identical
                     to B-factor values (ATOM records). None without ANISOU.
    "correct_uij"  : False if a non-standard combination of the Uij values in
                     the ANISOU records was necessary to reproduce the
                     B-factors.
    "format_date"  : PDB format date
    "format_vers"  : PDB format version
    "has_anisou"   : True if ANISOU records are present in the PDB file.
    "other_refinement_remarks"
                   : OTHER REFINEMENT REFMARKS in REMARK 3
                     as a single string.
    "prog_inter"   : a list with interpreted refinement programs(s) if not None
    "prog_last"    : a list with the guessed last-used refinement program. Can
                     be empty if no refinement programs were found in REMARK 3.
                     If the list is longer than 1 a decision about the programs
                     in the list could not be made.
    "prog_vers"    : a list with version of the
                     refinement programs(s) if not None
    "ref_prog"     : a list with refinement programs(s) if not None
    "refprog"      : the value of the refinement program record as a string.
    "is_bdb_includable"
                   : True if this PDB file is useful
    "req_tlsanl    : True if running TLSANL is required (also set for other
                     programs than REFMAC)
    "tls_groups"   : number of TLS groups as integer. If not present: None.
    "tls_residual" : True if it is mentioned in the TLS details or elsewhere
                     that the ATOM records contain residual B-factors only.
    "tls_sum"      : True if it is mentioned somewhere in REMARK 3 that the
                     ATOM records contain the sum of TLS and residual B-factors
    """
    is_bdb_includable = False
    assume_iso = False
    message = None
    _log.debug("Parsing refinement program...")

    # Parse the pdb records for refinement data
    format_vers, format_date = parse_format_date_version(pdb_records)
    other_refinement_remarks = parse_other_ref_remarks(pdb_records)
    pdb_info = {
            "pdb_id": pdb_id,
            "b_type": parse_btype(pdb_records),
            "has_anisou": "ANISOU" in pdb_records,
            "format_date": format_date,
            "format_vers": format_vers,
            "other_refinement_remarks": other_refinement_remarks,
            "b_msqav": is_bmsqav(other_refinement_remarks),
            "refprog": parse_ref_prog(pdb_records),
            "tls_groups": parse_num_tls_groups(pdb_records),
            "tls_residual": is_tls_residual(pdb_records),
            "tls_sum": is_tls_sum(pdb_records)}

    _log.debug("Interpreting PDB file...")

    # For entries that have ANISOU records..
    reproduced = {"beq_identical": None, "correct_uij": None}
    if pdb_info["has_anisou"]:
        reproduced = check_beq(structure)
        report_beq(reproduced)
        # ..we assume we can save time
        if reproduced["beq_identical"] > 0.9999:
            is_bdb_includable = True
            assume_iso = True
            message = "Assuming full isotropic B-factors because enough "\
                      "B-factors could be reproduced from the ANISOU records"
            _log.info(message)
    pdb_info.update(reproduced)

    prog = pdb_info["refprog"]
    prog_inter = None
    version = None
    prog_last = None
    req_tlsanl = False
    # Programs mentioned in REMARK 3
    if prog:
        # Interpret refinement program(s)
        prog, prog_inter, version = parse_refprog(prog)
        # Decide the final structure's most likely refprog signature
        prog_last = last_used(prog_inter, version)
        pdb_info.update({"prog_inter": prog_inter,
                         "prog_last": prog_last,
                         "prog_vers": version})

        for p, i, v in zip(prog, prog_inter, version):
            _log.debug("Refinement program: {} - interpreted as: {} - version:"
                       " {} - last used: {}.".format(p, i, v, prog_last))

        # Interpreted refinement program
        if prog:
            # (not enough) Beq values from ANISOU records are equal to B
            if not assume_iso:
                # interpret refinement data
                is_bdb_includable, assume_iso, req_tlsanl, message = \
                        decide_refprog(pdb_info)
            else:
                _log.info("{}: probably full B-factors.".format(
                    " and ".join(prog_last)))
        else:
            # we should not end up here under normal circumstances
            message = "Refinement program parse error"
            write_whynot(pdb_id, message)
            _log.error("{}.".format(message))
    else:
        message = "No refinement program found"
        _log.warn("{}.".format(message))
        if not assume_iso:
            write_whynot(pdb_id, message)

    more_refprog = {"assume_iso": assume_iso,
                    "decision": message,
                    "prog_inter": prog_inter,
                    "prog_last": prog_last,
                    "prog_vers": version,
                    "is_bdb_includable": is_bdb_includable,
                    "ref_prog": prog,
                    "req_tlsanl": req_tlsanl}
    pdb_info.update(more_refprog)
    return pdb_info


def last_used(pin, pv):
    """Make an educated guess about the refinement program that was used last.

    This function is primarily based on the following references:
    - Bauke Dijkstra, personal communication 2010, 2014 and
    - Supplementary material Touw & Vriend 2010 (Acta Cryst D66, 1341-50)
    and accompanying papers for rare program combinations

    Return an empty list if the input does not contain any refinement programs.
    Return a list of length 1 if the last program can be determined.
    Return a list with multiple programs if we cannot choose between them.

    Versions are not yet interpreted.
    """
    assert isinstance(pin, list)
    assert isinstance(pv, list)
    assert len(pin) == len(pv)
    # Remove non-refinement programs
    pin, pv = filter_progs(pin, pv)
    # e.g. X-PLOR 3.1, X-PLOR 3.8 in 1c4k
    pin = list(set(pin))
    # General rules
    if "SHELX" in pin and len(pin) > 1:  # multiple programs
        """
        If SHELXL/H has been used, it has probably been used last. (SHELXH was
        a version of SHELXL for very large structures that has become obsolete
        since newer versions of SHELXL)
        Examples were found of
        1. SHELX (without L/H) used for phasing only
        2. SHELXL was used for refinement and the depositers only mentioned
           SHELX (without L/H)
        There are more examples of situation 1. We therefore
        explicitly check for SHELXL/H and filter out SHELX if multiple programs
        are mentioned
        """
        shelx = pin.index("SHELX")
        if not re.search("[LH]", pv[shelx]):
            pin = pin[:shelx] + pin[(shelx + 1):]
        else:
            return ["SHELX"]
    if "CORELS" in pin and len(pin) > 1:
        """
        If CORELS has been used, it has probably been used first, unless
        it is the only program
        """
        corels = pin.index("CORELS")
        pin = pin[:corels] + pin[(corels + 1):]
    # Decide about some (sub)combinations
    #    First       Last
    known_combis = [
        ("CCP4", "CNS"),
        ("CCP4", "X-PLOR"),
        ("CNS", "BUSTER"),
        ("CNS", "PHENIX.REFINE"),
        ("CNS", "REFMAC"),
        ("CNS", "TNT"),
        ("CNX", "BUSTER"),
        ("CNX", "PHENIX.REFINE"),
        ("CNX", "REFMAC"),
        ("CNX", "TNT"),
        ("EREF", "PROLSQ"),
        ("EREF", "X-PLOR"),
        ("PROLSQ", "TNT"),
        ("PROLSQ", "REFMAC"),
        ("REFMAC", "MAIN"),  # According to currently available papers
        ("RESTRAIN", "X-PLOR"),
        ("TWIN_LSQ", "CNS"),
        ("TNT", "REFMAC"),
        ("X-PLOR", "BUSTER"),
        ("X-PLOR", "CNS"),
        ("X-PLOR", "PHENIX.REFINE"),
        ("X-PLOR", "PROLSQ"),  # First simulated annealing, then least squares
        ("X-PLOR", "REFMAC"),
        ("X-PLOR", "TNT"),
        ]
    for combi in known_combis:
        pin = one_of_the_two(pin, loser=combi[0], winner=combi[1])
    return pin


def one_of_the_two(lst, loser, winner):
    """Remove the loser from the list.

    Only if both loser and winner are present in the list lst.

    Return a list
    """
    assert isinstance(lst, list)
    if loser in lst and winner in lst and not loser == winner:
        li = lst.index(loser)
        lst = lst[:li] + lst[(li + 1):]
    return lst


def parse_refprog(refprog):
    """Parse and sort the refinement program(s) found in the PDB file.

    Return a tuple of program(s), interpreted program(s) and version(s).

    Multiple refinement programs
    - ";" (3uec) and "+" (3fyd, 4e9j) have been used as delimiters
    - "," is normally used to split multiple programs
      exceptions: 1lsh: CNS 0.9,1.0,1.1
                  1ksj: BUSTER, BETA VERSION
                  3cw1: O, VERSION 9.0.7
    - "&" is almost always used to split multiple programs
      exceptions: 3rvz, 3rvy, 3rw0: CNS 1.1 & 1.3
    - "AND" is almost always used to split multiple programs
      exception: 1qjx: X-PLOR 3.1 AND 3.85
    - "/" is sometimes used to split multiple programs
      exceptions: ARP/WARP, SOLVE/RESOLVE, BUSTER/TNT (normally BUSTER-TNT),
      1amq/1amr: TOM/FRODO
      1etl/1etn/1etm: XPLOR, FMLS/VP
      1onq: REFMAC 5.1.24/TLS
      1w5r: REFMAC 5.2.005 24/04/2001
      2jjy: REFMAC 5.2.019 24/04/2001
    - "_" has been normally been used in PHENIX.REFINE version strings and only
      once for multiple programs:
      3h86: PHENIX.REFINE_REFMAC 5.5.0070
    - seperated by space (exception):
      1cq1: REFMAC X-PLOR 3.843

    - Versions: see below

    - Programs: see below

    Note: it is assumed PHENIX means PHENIX.REFINE, other types such as
          ENSEMBLE_REFINEMENT are stored seperately if we know about them
          (otherwise they will appear as version string)
    Note: it is assumed the same refinement program is listed only once

    """
    # Exceptions
    if refprog == "NULL" or refprog == "NONE" or \
            refprog == "NO REFINEMENT":
        except_refprog_warn()
        return ([None], [None], [None])
    elif refprog == "X-PLOR 3.1 AND 3.85":
        except_refprog_warn()
        return (["X-PLOR 3.85"], ["X-PLOR"], ["3.85"])
    elif refprog == "X-PLOR 3.1, 3.816":
        except_refprog_warn()
        return (["X-PLOR 3.816"], ["X-PLOR"], ["3.816"])
    elif refprog == "CNS 1.1 & 1.3":
        except_refprog_warn()
        return (["CNS 1.3"], ["CNS"], ["1.3"])
    elif refprog == "CNS 0.4, O, OOPS":
        except_refprog_warn()
        return (["CNS 0.4", "O", "OOPS"],
                ["CNS", "O", "OOPS"],
                ["0.4", "-", "-"])
    elif refprog == "CNS 0.1-0.4":
        except_refprog_warn()
        return (["CNS 0.4"], ["CNS"], ["0.4"])
    elif refprog == "CNS 0.9,1.0,1.1":
        except_refprog_warn()
        return (["CNS 1.1"], ["CNS"], ["1.1"])
    elif refprog == "CNS 1.3 WITH DEN REFINEMENT":
        except_refprog_warn()
        return (["CNS 1.3"], ["CNS"], ["1.3"])
    elif refprog == "CNS 1.2 (USING XTAL_TWIN UTILITIES)":
        except_refprog_warn()
        return (["CNS 1.2"], ["CNS"], ["1.2"])
    elif refprog == "PHENIX.REFINE_REFMAC 5.5.0070":
        except_refprog_warn()
        return (["PHENIX.REFINE", "REFMAC 5.5.0070"],
                ["PHENIX.REFINE", "REFMAC"],
                ["-", "5.5.0070"])
    elif refprog == "PHENIX (CCI APPS 2007_04_06_1210)":
        except_refprog_warn()
        return (["PHENIX (PHENIX.REFINE: 2007_04_06_1210)"], ["PHENIX.REFINE"],
                ["2007_04_06_1210"])
    elif refprog == "PHENIX VERSION 1.8_1069 (PHENIX.REFINE)":
        except_refprog_warn()
        return (["PHENIX (PHENIX.REFINE: 1.8_1069)"], ["PHENIX.REFINE"],
                ["1.8_1069"])
    elif refprog == "PHENIX 1.6.2_432 - REFINE":
        except_refprog_warn()
        return (["PHENIX (PHENIX.REFINE: 1.6.2_432)"], ["PHENIX.REFINE"],
                ["1.6.2_432"])
    elif refprog == "PHENIX REFINE" or refprog == "PHENIX, REFINE" or \
            refprog == "PHENIX AUTOREFINE":
        except_refprog_warn()
        return (["PHENIX.REFINE"], ["PHENIX.REFINE"], ["-"])
    elif refprog == "REFMAC 5.1.24/TLS":
        except_refprog_warn()
        return (["REFMAC 5.1.24"], ["REFMAC"], ["5.1.24"])
    elif refprog == "REFMAC 5.2.0005 24/04/2001":
        except_refprog_warn()
        return (["REFMAC 5.2.0005"], ["REFMAC"], ["5.2.0005"])
    elif refprog == "REFMAC 5.2.0019 24/04/2001":
        except_refprog_warn()
        return (["REFMAC 5.2.0019"], ["REFMAC"], ["5.2.0019"])
    elif refprog == "REFMAC X-PLOR 3.843":
        except_refprog_warn()
        return (["REFMAC", "X-PLOR 3.843"],
                ["REFMAC", "X-PLOR"],
                ["-", "3.843"])
    elif refprog == "REFMAC5 5.2.0019":
        except_refprog_warn()
        return (["REFMAC 5.2.0019"], ["REFMAC"], ["5.2.0019"])
    elif refprog == "REFMAC 5.5.0109 (AND PHENIX)":
        except_refprog_warn()
        return (["REFMAC 5.5.0109", "PHENIX.REFINE"],
                ["REFMAC", "PHENIX.REFINE"],
                ["5.5.0109", "-"])
    elif refprog == "BUSTER, BETA VERSION":
        except_refprog_warn()
        return (["BUSTER BETA"], ["BUSTER"], ["BETA"])
    elif refprog == "TNT BUSTER/TNT":
        except_refprog_warn()
        return (["TNT", "BUSTER/TNT"], ["TNT", "BUSTER"], ["-", "-"])
    elif refprog == "O, VERSION 9.0.7":
        except_refprog_warn()
        return (["O 9.0.7"], ["O"], ["9.0.7"])
    prog_pat = re.compile(r"""
        ,
        |&
        |;
        |\+
        (?!SVN)       # 4ow3 has PHENIX (PHENIX.REFINE: DEV_1549+SVN)
        |AND
        |
        (?<!ARP)      # Note that if the program(-part)s before or after
        (?<!SOLVE)    # the slash have been combined with a different part
        (?<!BUSTER)   # after or before the slash (respectively) a split
        (?<!TOM)      # is made as well.
        (?<!FMLS)
        /             # We need negative lookaheads and -behinds for a / split
        (?!WARP
        |RESOLVE
        |TNT
        |FRODO
        |VP)
        """, re.VERBOSE)
    prog = prog_pat.split(refprog)
    prog = filter(None, prog)
    prog = [p.strip(" ").upper() for p in prog]
    # Do not sort here, we might want to use the given order of progs later
    # prog = sorted(prog)
    prog_inter = [None] * len(prog)
    vers = [None] * len(prog)
    no_versions = [
        "ARP",
        "BILDER",
        "CCP4",
        "CEDAR",
        "CHAIN",
        "CORELS",
        "DM",
        "DYNAMIX",
        "EREF",
        "FFX",
        "FMLS/VP",
        "FRODO",
        "HIPHOP",
        "IMPLOR",
        "LAFIRE",
        "LALS",
        "MAIN",
        "MOLPROBITY",
        "MOLLY",
        "MOPRO",
        "NCNS",
        "NCNS-TINKER",
        "NMREF",
        "PIKSOL",
        "PHASER",
        "PMB",
        "POLYVISION",
        "PRIMEX",
        "PRODRG",
        "PROTEIN",
        "QUANTA",
        "RESTRAIN",
        "SCWRL",
        "SHARP",
        "SFALL",
        "SOLVE",
        "SOLVE/RESOLVE",
        "TIBBITTS",
        "TOM",
        "TOM/FRODO",
        "XFIT",
        "XPLEO",
        "XTALVIEW"]
    for i, p in enumerate(prog):
        # ... and yes, we also find REFAMC (3m1o) and REFMEC (3e9q)...
        if re.match("[REFMAC]{6}", p):
            prog_inter[i] = "REFMAC"
            # if one of the programs is REFMAC, we expect
            # "REFMAC <version>"
            # "REFMAC<version>"
            # "REFMAC V <version>"
            # ...
            s0 = re.search(r"^[REFMAC]{6} ?(?:V )?([.0-9A-Z]+)$", p)
            # ... or "REFMAC"
            s1 = re.search(r"^[REFMAC]{6}\s*$", p)
            if s0:
                if s0.group(1) != "NULL":
                    vers[i] = s0.group(1)
            elif s1:
                vers[i] = "-"
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.match("CNS", p):
            prog_inter[i] = "CNS"
            # if one of the programs is CNS, we expect
            # "CNS (<version>)"
            # ...
            s0 = re.search(r"^CNS \(([.0-9A-Z]+)\)$", p)
            # ... or (more commonly)
            # "CNS <version>"
            # "CNS<version>"
            # "CNS V. <version>"
            # "CNS-<version>"
            # (for safety we keep this regex seperate from the previous one)
            s1 = re.search(r"^CNS(?:[ -])?(?:V. )?([.0-9A-Z]+)$", p)
            # ... or "CNS"
            s2 = re.search(r"^CNS\s*$", p)
            if s0:
                if s0.group(1) != "NULL":
                    vers[i] = s0.group(1)
            elif s1:
                if s1.group(1) != "NULL":
                    vers[i] = s1.group(1)
            elif s2:
                vers[i] = "-"
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.match("TWIN_LSQ", p):
            prog_inter[i] = "CNS"
            vers[i] = "TWIN_LSQ" if re.search(r"^TWIN_LSQ$", p) else "np"
        elif re.match("CNX", p):
            prog_inter[i] = "CNX"
            # (ACCELRYS) is not useful
            p = re.sub(" ?\(ACCELRYS\)", "", p)
            # if one of the programs is CNX, we expect
            # "CNX <version>"
            # "CNX<version>"
            s0 = re.search(r"^CNX ?([.0-9\-]+)$", p)
            # ... or "CNS"
            s1 = re.search(r"^CNX\s*$", p)
            if s0:
                if s0.group(1) != "NULL":
                    vers[i] = s0.group(1)
            elif s1:
                vers[i] = "-"
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.match("X-?PLOR", p):
            prog_inter[i] = "X-PLOR"
            # (ONLINE) is not useful
            p = re.sub(" ?\(ONLINE\)", "", p)
            # if one of the programs is X-PLOR, we expect
            # "X-PLOR <version>"
            s0 = re.search(r"^X-PLOR ([.0-9A-Z]+)\s*$", p)
            # ... or without version
            # "X-PLOR"
            # "XPLOR"
            s1 = re.search(r"^X-?PLOR\s*$", p)
            if s0:
                if s0.group(1) != "NULL":
                    vers[i] = s0.group(1)
            elif s1:
                vers[i] = "-"
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.search(r"PHENIX.ENSEMBLE_REFINEMENT", p):
            prog_inter[i] = "PHENIX.ENSEMBLE_REFINEMENT"
            # if one of the programs is PHENIX.ENSEMBLE_REFINEMENT, we expect
            # "PHENIX (PHENIX.ENSEMBLE_REFINEMENT: <version>)"
            s = re.search(
                r"^PHENIX \(PHENIX.ENSEMBLE_REFINEMENT: ([.\-_0-9A-Z]+)\)$",
                p)
            if s:
                if s.group(1) != "NULL":
                    vers[i] = s.group(1)
            # ... and otherwise we cannot (yet) handle this version format
            # and don't want to
            else:
                vers[i] = "np"
        elif re.match("PHENIX", p):
            prog_inter[i] = "PHENIX.REFINE"
            # if one of the programs is PHENIX.REFINE, we expect
            # "PHENIX <version>"
            # "PHENIX.REFINE: <version>"
            # ...
            s0 = re.search(r"^PHENIX(?:.REFINE:)? ([.\-_0-9A-Z]+)$", p)
            # ... or
            # "PHENIX (<version>)"
            # "PHENIX (PHENIX.REFINE: <version>)"
            s2 = re.search(
                r"^PHENIX \((?:PHENIX.REFINE: )?([.\-_0-9A-Z]+)\)$", p)
            # ... or without version
            # "PHENIX"
            # "PHENIX.REFINE"
            # "PHENIX.REFINEMENT"
            # "PHENIX (PHENIX.REFINE)"
            # "PHENIX.REFINE (PHENIX.REFINE)" is also matched
            # "PHENIX (PHENIX)" as well
            s1 = re.search(
                r"^PHENIX(?:.REFINE(?:MENT)?)?(?: \(PHENIX(.REFINE)?\))?$", p)
            if s0:
                if s0.group(1) != "NULL":
                    vers[i] = s0.group(1)
            elif s1:
                vers[i] = "-"
            elif s2:
                if s2.group(1) != "NULL":
                    vers[i] = s2.group(1)
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.search("BUSTER", p):
            prog_inter[i] = "BUSTER"
            # if one of the programs is BUSTER, we expect
            # "AUTOBUSTER <version>"
            # "BUSTER <version>"
            # "BUSTER-TNT <version>"
            # "BUSTER-TNT V. <version>"
            # "BUSTER-TNT BUSTER <version>"
            # combinations are also matched
            # <version> : [0-9.X]
            s0 = re.search(r"^(?:AUTO)?BUSTER(?:-TNT)? (?:BUSTER )?"
                           "(?:V. )?([0-9.X]+)$", p)
            # ... or without version
            # "AUTOBUSTER"
            # "BUSTER"
            # "BUSTER TNT"
            # "BUSTER-TNT"
            # "BUSTER/TNT"
            s1 = re.search(r"^(?:AUTO)?BUSTER(?:[ \-/]TNT)?\s*$", p)
            if s0:
                if s0.group(1) != "NULL":
                    vers[i] = s0.group(1)
            elif s1:
                vers[i] = "-"
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.match("TNT", p):
            prog_inter[i] = "TNT"
            # if one of the programs is TNT, we expect
            # "TNT <version>"
            # "TNT V. <version>"
            # "TNT V. <version> PRERELEASE"
            # "TNT <version> PRERELEASE" is also matched
            s0 = re.search(r"^TNT (?:V. )?([\-.0-9A-Z]+)( PRERELEASE)?\s*$", p)
            # ... or without version
            # "TNT"
            s1 = re.search(r"^TNT\s*$", p)
            if s0:
                if s0.group(1) != "NULL":
                    if s0.group(2):
                        vers[i] = s0.group(1) + s0.group(2)
                    else:
                        vers[i] = s0.group(1)
            elif s1:
                vers[i] = "-"
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.match("SHELX", p):
            prog_inter[i] = "SHELX"
            # if one of the programs is SHELX, we expect
            # without version:
            # "SHELX"
            # ...
            s0 = re.search(r"^SHELX\s*$", p)
            # ... or
            # "SHELX<version>":
            # SHELX[HLS]
            # SHELX-[0-9]{2,4}
            # SHELXL-[0-9]{2,4}
            # SHELXH-[0-9]{2,4}
            # SHELX-97-1
            # SHELX 97-1
            # SHELXL 97
            # SHELXL97
            # SHELX-L
            s1 = re.search(r"^SHELX[- ]?([HLS]?[\- ]?([0-9\-]{2,4})?)$", p)
            if s0:
                vers[i] = "-"
            elif s1:
                if s1.group(1) != "NULL":
                    vers[i] = s1.group(1)
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.match("ARP/WARP", p):
            prog_inter[i] = "ARP/WARP"
            # if one of the programs is ARP/WARP, we expect
            # "ARP/WARP V. <version>"
            # ...
            s0 = re.search(r"^ARP/WARP V. ([.0-9]+)$", p)
            # ... or without version
            # "ARP/WARP"
            s1 = re.search(r"^ARP/WARP$", p)
            if s0:
                if s0.group(1) != "NULL":
                    vers[i] = s0.group(1)
            elif s1:
                vers[i] = "-"
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.match("COOT", p):
            prog_inter[i] = "COOT"
            # if one of the programs is COOT, we expect
            # "COOT <version>"
            # "COOT<version>"
            # "COOT V. <version>"
            # <version> : [0-9.\-]+(-PRE-[0-9]+)?
            s0 = re.search(r"^COOT ?(?:V. )?([0-9.\-]+(-PRE-[0-9]+)?)\s*$", p)
            # ... or without version
            # "COOT"
            s1 = re.search(r"^COOT\s*$", p)
            if s0:
                if s0.group(1) != "NULL":
                    vers[i] = s0.group(1)
            elif s1:
                vers[i] = "-"
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.match("O", p):
            prog_inter[i] = "O"
            # if one of the programs is O, we expect
            # "O<version>"
            # <version> : [0-9.]+
            s0 = re.search(r"^O([0-9.]+)\s*$", p)
            # ... or without version
            # "O"
            s1 = re.search(r"^O\s*$", p)
            if s0:
                if s0.group(1) != "NULL":
                    vers[i] = s0.group(1)
            elif s1:
                vers[i] = "-"
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.match("PROFFT", p):
            """
            PROTIN and NUCLIN set up restraints for
            protein refinement with PROLSQ/PROFFT and
            nucleic acid refinement with NUCLSQ, respectively,

            All of these programs will be interpreted as PROLSQ with versions
            PROLSQ: PROLSQ
            PROFFT: PROFFT
            PROTIN: PROLSQ
            NUCLIN: NUCLSQ
            NUCLSQ: NUCLSQ
            GPRLSA: GPRLSA
            and possible modifications
            """
            prog_inter[i] = "PROLSQ"
            # if one of the programs is PROFFT, we expect
            # "PROFFT (MODIFIED BY Z.OTWINOWSKI)" (1 case)
            # ...
            s0 = re.search(r"^PROFFT \(MODIFIED BY Z.OTWINOWSKI\)$", p)
            # ... or without version
            # "PROFFT"
            s1 = re.search(r"^PROFFT$", p)
            if s0:
                vers[i] = "PROFFT MODIFIED BY Z.OTWINOWSKI"
            elif s1:
                vers[i] = "PROFFT"
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.match("PROLSQ", p):
            prog_inter[i] = "PROLSQ"
            # if one of the programs is PROLSQ, we expect
            # "PROLSQ (MODIFIED BY G.J.QUIGLEY)" (8 cases)
            # ...
            s0 = re.search(r"^PROLSQ \(MODIFIED BY G.J.QUIGLEY\)$", p)
            # ... or without version
            # "PROLSQ"
            s1 = re.search(r"^PROLSQ$", p)
            if s0:
                vers[i] = "PROLSQ MODIFIED BY G.J.QUIGLEY"
            elif s1:
                vers[i] = "PROLSQ"
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.match("PROTIN", p):
            prog_inter[i] = "PROLSQ"
            vers[i] = "PROLSQ" if re.search(r"^PROTIN$", p) else "np"
        elif re.match("NUCLSQ", p):
            prog_inter[i] = "PROLSQ"
            # if one of the programs is NUCLSQ, we expect
            # "NUCLSQ (MODIFIED BY G.J.QUIGLEY)" (1 case)
            # ...
            s0 = re.search(r"^NUCLSQ \(MODIFIED BY G.J.QUIGLEY\)$", p)
            # ... or without version
            # "NUCLSQ"
            s1 = re.search(r"^NUCLSQ$", p)
            if s0:
                vers[i] = "NUCLSQ MODIFIED BY G.J.QUIGLEY"
            elif s1:
                vers[i] = "NUCLSQ"
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.match("NUCLIN", p):
            prog_inter[i] = "PROLSQ"
            vers[i] = "NUCLSQ" if re.search(r"^NUCLIN$", p) else "np"
        elif re.match("GPRLSA", p):
            prog_inter[i] = "PROLSQ"
            vers[i] = "GPRLSA" if re.search(r"^GPRLSA$", p) else "np"
        elif re.match("DERIV", p):
            prog_inter[i] = "PROLSQ"
            vers[i] = "DERIV" if re.search(r"^DERIV$", p) else "np"
        elif re.match("CERIUS", p):
            prog_inter[i] = "CERIUS"
            # We match CERIUS <version>
            #      and CERIUS<version>
            # <version>: [0-9.\-]+
            # ...
            s0 = re.search(r"^CERIUS ?([0-9.\-]+)$", p)
            # ... or CERIUS without version
            s1 = re.search(r"^CERIUS$", p)
            if s0:
                if s0.group(1) != "NULL":
                    vers[i] = s0.group(1)
            elif s1:
                vers[i] = "-"
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        elif re.match("HKL-?3000", p):
            prog_inter[i] = "HKL-3000"
            # if one of the programs is HKL-3000, we expect
            # "HKL-3000"
            # "HKL3000"
            s = re.search(r"^HKL-?3000$", p)
            vers[i] = "-" if s else "np"
        elif re.match("GROMOS", p):
            prog_inter[i] = "GROMOS"
            # if one of the programs is GROMOS, we expect
            # "GROMOS<version>"
            # ...
            s0 = re.search(r"^GROMOS([0-9]+)$", p)
            # ... or without version
            # "GROMOS"
            s1 = re.search(r"^GROMOS$", p)
            if s0:
                if s0.group(1) != "NULL":
                    vers[i] = s0.group(1)
            elif s1:
                vers[i] = "-"
            # ... and otherwise we cannot (yet) handle this version format
            else:
                vers[i] = "np"
        else:
            other = True
            for r in no_versions:
                # Now we only distinguish RESTRAIN from RESTRAINED
                # Other exceptions will be logged
                if re.search(r"^" + r + "(?!ED)", p):
                    prog_inter[i] = r
                    vers[i] = "-" if re.search(r"^" + r + "$", p) else "np"
                    other = False
            if other:
                # AGARWAL FAST-FOURIER TRANSFORM LEAST-SQUARES
                # HENDRICKSON-KONNERT LEAST-SQUARES REFINEMENT
                # * OF T. A. JONES.  THE R VALUE IS 0.180.
                # RESTRAINED RECIPROCAL-SPACE LEAST-SQUARES
                # SIMULATED ANNEALING METHOD
                # CONSTRAINED RECIPROCAL-SPACE LEAST-SQUARES
                # FAST-FOURIER LEAST-SQUARES REFINEMENT
                # JACK-LEVITT
                # REAL-SPACE REFINEMENT
                prog_inter[i] = "OTHER"
                vers[i] = "np"
        # Report
        if prog_inter[i] == "OTHER":
            _log.warn("{}: program {} could not (yet) be parsed.".format(
                prog_inter[i], p))
        elif vers[i] == "np":
            _log.warn("{}: version could not (yet) be parsed.".format(
                prog_inter[i]))
        elif vers[i] == "-":
            _log.debug("{}: version not present.".format(prog_inter[i]))
    return prog, prog_inter, vers

