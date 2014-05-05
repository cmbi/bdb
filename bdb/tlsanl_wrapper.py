#!/usr/bin/env python
from __future__ import print_function

import argparse
import logging
import os
import subprocess
import sys

from bdb.bdb_utils import write_whynot


_log = logging.getLogger(__name__)


def run_tlsanl(pdb_file_path, xyzout, pdb_id, log_out_dir=".",
               verbose_output=False):
    """Run TLSANL.

    A REFMAC file with residual isotropic B-factors and proper TLS descriptions
    is expected. Total isotropic B-factors are written out in the ATOM and
    ANISOU records.

    WARNING: it is assumed that ATOM & HETATM records in the input PDB must
    first be sorted on chain ID and residue number before the TLS ranges
    can be interpreted.

    Detailed documentation for TLSANL can be found at
    http://www.ccp4.ac.uk/html/tlsanl.html.
    """
    _log.info("Preparing TLSANL run...")
    success = False
    keyworded_input = "BINPUT t\nBRESID t\nISOOUT FULL\nNUMERIC\nEND\n"
    p = subprocess.Popen(["tlsanl", "XYZIN", pdb_file_path, "XYZOUT", xyzout],
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate(input=keyworded_input)
    try:
        with open(os.path.join(log_out_dir, "tlsanl.log"), "w") as tlsanl_log:
            tlsanl_log.write(stdout)
            if verbose_output:
                print(stdout)
        with open(os.path.join(log_out_dir, "tlsanl.err"), "w") as tlsanl_err:
            tlsanl_err.write(stderr)
            if verbose_output:
                print(stderr)
    except IOError as ex:
        _log.error(ex)
    if p.returncode != 0:
        message = "TLSANL problem (exit code: {0:3d})".format(p.returncode)
        write_whynot(pdb_id, message)
        _log.error("{0:s}".format(message))
    elif os.stat(xyzout).st_size <= 2000:
        # from script at http://deposit.rcsb.org/adit/REFMAC.html
        message = "TLSANL problem"
        write_whynot(pdb_id, message)
        _log.error("{0:s}".format(message))
    elif os.stat(os.path.join(log_out_dir, "tlsanl.err")).st_size > 0:
        message = "TLSANL problem"
        write_whynot(pdb_id, message)
        _log.error("{0:s}".format(message))
    else:
        success = True
        _log.info("TLSANL ran without problems.")
    return success

