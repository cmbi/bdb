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
from __future__ import print_function

import logging
_log = logging.getLogger(__name__)

import os
import pyconfig
import re
import subprocess

from pdbb.bdb_utils import write_whynot


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
        with open(os.path.join(log_out_dir, pyconfig.get("TLSANL_LOG")),
                  "w") as tlsanl_log:
            tlsanl_log.write(stdout)
            if verbose_output:
                print(stdout)
        with open(os.path.join(log_out_dir, pyconfig.get("TLSANL_ERR")),
                  "w") as tlsanl_err:
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
    elif os.stat(os.path.join(log_out_dir,
                              pyconfig.get("TLSANL_ERR"))).st_size > 0:
        message = "TLSANL problem"
        write_whynot(pdb_id, message)
        _log.error("{0:s}".format(message))
    else:
        success = True
        _log.info("TLSANL ran without problems.")
    return success


def parse_skttls_summ(tlsanl_log):
    """Parse Skttles summary from TLSANL log file

    Return the total number of bonds between residues and the number of bonds
    beyond the 95th and 99th percentile for any residual as a dict of integers.
    """
    skttls = {"skttls_tot": None,
              "skttls_95th": None,
              "skttls_99th": None}

    RE_SKTTLS_TOT = re.compile(r"^#  Total number of bonds between residues:"
                               "\s*(?P<num>\d+)")
    RE_SKTTLS_95th = re.compile(r"^#  Number of bonds beyond 95th percentile "
                                "for any residual:\s*(?P<num>\d+)")
    RE_SKTTLS_99th = re.compile(r"^#  Number of bonds beyond 99th percentile "
                                "for any residual:\s*(?P<num>\d+)")

    with open(tlsanl_log, "r") as log:
        lines = [l.strip() for l in log]
        for l in lines:
            m = RE_SKTTLS_TOT.search(l)
            if m is not None:
                skttls["skttls_tot"] = int(m.group("num"))
            m = RE_SKTTLS_95th.search(l)
            if m is not None:
                skttls["skttls_95th"] = int(m.group("num"))
            m = RE_SKTTLS_99th.search(l)
            if m is not None:
                skttls["skttls_99th"] = int(m.group("num"))

    return skttls
