from __future__ import print_function

import logging
import os

from pdbb.pdb.parser import parse_exp_methods
from pdbb.bdb_utils import write_whynot


_log = logging.getLogger(__name__)


def check_exp_methods(pdb_records, pdb_id):
    """
    Check if the experiment methods are suitable for adding the PDB file to
    the BDB.

    Returns a dict such that:
        "expdta_useful" : True if this PDB file is useful
        "expdta"        : a list with experimental method(s) if not None
    """

    _log.debug("Parsing EXPDTA...")

    # Parse the experiment methods. A ValueError is raised when no EXPDTA
    # records are found, which is fatal.
    try:
        exp_methods = parse_exp_methods(pdb_records)
    except ValueError:
        message = "Experimental method: EXPDTA parse error"
        write_whynot(pdb_id, message)
        _log.error("{}.".format(message))
        return {"expdta_useful": False, "expdta": []}

    # Multiple experiment methods are not supported
    if len(exp_methods) > 1:
        message = "Experimental method: multiple ({0})".format(
            " and ".join(exp_methods))
        write_whynot(pdb_id, message)
        _log.warn("{}.".format(message))
        return {"expdta_useful": False, "expdta": exp_methods}

    assert len(exp_methods) == 1

    _log.info("Experimental method: {0}.".format(exp_methods[0]))

    useful = False
    if "X-RAY DIFFRACTION" in exp_methods:
        useful = True
    else:
        message = "Experimental method: " + exp_methods[0]
        write_whynot(pdb_id, message)
        _log.warn("{} cannot be included in the bdb.".format(message))

    return {"expdta_useful": useful, "expdta": exp_methods}
