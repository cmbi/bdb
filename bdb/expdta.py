from __future__ import print_function

import argparse
import logging
import os
import re
import sys

from bdb.pdb.parser import parse_exp_methods
from bdb.bdb_utils import init_bdb_logger, write_whynot


_log = logging.getLogger("bdb")


def check_exp_methods(pdb_records, pdb_id, out_dir=".", global_files=False):
    """
    Check if the experiment methods are suitable for adding the PDB file to
    the BDB.

    Returns a dict such that:
        "expdta_useful" : True if this PDB file is useful
        "expdta"        : a list with experimental method(s) if not None
    """

    _log.debug(("{0:4s} | Parsing EXPDTA...").format(pdb_id))

    # Parse the experiment methods. A ValueError is raised when no EXPDTA
    # records are found, which is fatal.
    try:
        exp_methods = parse_exp_methods(pdb_records)
    except ValueError:
        message = "Experimental method: EXPDTA parse error"
        write_whynot(pdb_id, message, directory=out_dir)
        _log.error(("{0:4s} | {1:s}.").format(pdb_id, message))
        return {"expdta_useful": False, "expdta": []}

    # Multiple experiment methods are not supported
    if len(exp_methods) > 1:
        message = "Experimental method: multiple (" + " and ".join(exp_methods) + ")"
        write_whynot(pdb_id, message, directory=out_dir)
        _log.warn(("{0:4s} | {1:s}.").format(pdb_id, message))
        if global_files:
            write_unsupported_expdta(exp_methods)
        return {"expdta_useful": False, "expdta": exp_methods}

    assert len(exp_methods) == 1

    _log.info(("{0:4s} | Experimental method:{1:s}.").format(
        pdb_id, exp_methods[0]))

    useful = False
    if "X-RAY DIFFRACTION" in exp_methods:
        useful = True
    else:
        message = "Experimental method: " + exp_methods[0]
        write_whynot(pdb_id, message, directory=out_dir)
        _log.warn(("{0:4s} | {1:s} cannot be included in the bdb.").format(
            pdb_id, message))
        if global_files:
            write_unsupported_expdta(exp_methods)

    return {"expdta_useful": useful, "expdta": exp_methods}


def write_unsupported_expdta(expdta, expdta_file="unsupported_expdta.txt"):
    """Add EXPDTA to the unique list of unsupported EXPDTA if applicable."""
    unsupported = []
    try:
        if os.path.exists(expdta_file):
            with open(expdta_file, "r") as f:
                unsupported = [line.strip() for line in f]
        unsupported.append(" and ".join(expdta))
        unsupported = set(unsupported)
        with open(expdta_file, "w") as f:
            for u in unsupported:
                f.write("{0:s}\n".format(u))
    except IOError as ex:
        _log.error(ex)
