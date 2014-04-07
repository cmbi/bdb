#!/usr/bin/env python
from __future__ import print_function
from bdb_utils import get_raw_pdb_info, init_bdb_logger,\
unique, write_whynot, PDB_LOGFORMAT
import argparse
import logging
import os
import re
import sys

# Configure logging
_log = logging.getLogger("bdb")

def check_expmethod(method):
    """Check if the experimental method can be included in the bdb.

    Return a Boolean.
    """
    return(method == "X-RAY DIFFRACTION")

def decide_expdta(expdta, pdb_id, out_dir=".", global_files=False):
    """Determine whether this EXPDTA list can be included in the bdb.

    Return a Boolean.
    """
    useful = False
    if len(expdta) == 1:
        _log.info(("{0:" + PDB_LOGFORMAT + "} | Experimental method: {1:s}.").
                  format(pdb_id, expdta[0]))
        if check_expmethod(expdta[0]):
            useful = True
        else:
            message = "Experimental method: " + expdta[0]
            write_whynot(pdb_id, message, directory=out_dir)
            _log.warn(("{0:" + PDB_LOGFORMAT + "} | {1:s} cannot be "\
                       "included in the bdb.").format(pdb_id, message))
            if global_files:
                write_unsupported_expdta(expdta)
    elif len(expdta) > 1:
        message = "Experimental method: multiple (" + " and ".join(expdta) +\
                  ")"
        write_whynot(pdb_id, message, directory=out_dir)
        _log.warn(("{0:" + PDB_LOGFORMAT + "} | {1:s}.").format(pdb_id,
                                                               message))
        if global_files:
            write_unsupported_expdta(expdta)
    else: # we should not end up here
        message = "Experimental method: EXPDTA parse error"
        write_whynot(pdb_id, message, directory=out_dir)
        _log.error(("{0:" + PDB_LOGFORMAT + "} | {1:s}.").format(pdb_id,
                                                               message))
    return(useful)

def do_expdta(pdb_xyz, pdb_id=None, out_dir=".", global_files=False):
    """Determine whether this PDB file can be included in the bdb.

    Return a dictionary.
    "expdta_useful" : True if this PDB file is useful
    "expdta"        : a list with experimental method(s) if not None
    """
    pdb_id = pdb_xyz if pdb_id is None else pdb_id
    success = False
    _log.debug(("{0:" + PDB_LOGFORMAT + "} | "\
                "Parsing EXPDTA...").format(pdb_id))
    expdta = get_expdta(pdb_xyz)
    if expdta:
        expdta = parse_expdta(expdta)
        _log.debug(("{0:" + PDB_LOGFORMAT + "} | Experimental method(s):"\
                    "{1:s}.").format(pdb_id, " and ".join(expdta)))
        if expdta:
            if decide_expdta(expdta, pdb_id, out_dir, global_files):
                success = True
        else: # we should not end up here
            message = "Experimental method: EXPDTA parse error"
            write_whynot(pdb_id, message, directory=out_dir)
            _log.error(("{0:" + PDB_LOGFORMAT + "} | {1:s}.").format(pdb_xyz,
                                                                    message))
    else:
        message = "Experimental method: no EXPDTA record"
        write_whynot(pdb_id, message, directory=out_dir)
        _log.warn(("{0:" + PDB_LOGFORMAT + "} | {1:s}.").format(pdb_xyz,
                                                               message))
    return({"expdta_useful": success, "expdta":expdta})


def get_expdta(pdb_xyz):
    """Find the EXPDTA record in a PDB file.

    Return the value of the EXPDTA record as a string.
    """
    expdta = get_raw_pdb_info(pdb_xyz)["expdta"]
    if expdta:
        _log.debug(("{0:" + PDB_LOGFORMAT + "} | {1:s}.").
                   format(pdb_xyz, expdta))
    return(expdta)

def parse_expdta(expdta):
    """Parse and sort the experimental method(s) in the EXPDTA record.

    Return a list.
    """
    method = expdta.split(";")
    method = filter(None, method)
    method = [m.strip(" ").upper() for m in method]
    method = sorted(method)
    return(method)

def write_unsupported_expdta(expdta, expdta_file="unsupported_expdta.txt"):
    """Add EXPDTA to the unique list of unsupported EXPDTA if applicable."""
    unsupported = []
    try:
        if os.path.exists(expdta_file):
            with open(expdta_file, "r") as f:
                unsupported = [line.strip() for line in f]
        unsupported.append(" and ".join(expdta))
        unsupported = unique(unsupported)
        with open(expdta_file, "w") as f:
            for u in unsupported:
                f.write("{0:s}\n".format(u))
    except IOError as ex:
        _log.error(ex)

if __name__ == "__main__":
    """Process the EXPDTA record from a PDB file.

    Exit with an exit code of 1 if the PDB file has no EXPDTA record, if it
    cannot be parsed, if it contains multiple methods or if the method cannot
    be used in the bdb project.
    """
    parser = argparse.ArgumentParser(description="Parse EXPDTA")
    parser.add_argument("-v", "--verbose", help="verbose output",
                        action="store_true")
    parser.add_argument("--pdbid", help="PDB file name.")
    parser.add_argument("xyzin", help="Input coordinates in PDB format.")
    args = parser.parse_args()
    pdb_id = args.pdbid if args.pdbid is not None else args.xyzin
    _log = init_bdb_logger(pdb_id, global_log=True)
    if args.verbose:
        _log.setLevel(logging.DEBUG)
    if do_expdta(args.xyzin, pdb_id)["expdta_useful"]:
        sys.exit(0)
    else:
        sys.exit(1)
