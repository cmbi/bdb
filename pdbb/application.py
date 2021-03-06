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

import argparse
import json
import os
import pyconfig
import shutil

from pdbb.bdb_utils import (is_valid_directory, is_valid_file, is_valid_pdbid,
                            get_bdb_entry_outdir, write_whynot)
from pdbb.check_beq import (determine_b_group, get_structure,
                            write_multiplied_8pipi)
from pdbb.expdta import check_exp_methods
from pdbb.pdb.parser import parse_pdb_file
from pdbb.refprog import get_refi_data
from pdbb.requirements import check_deps
from pdbb.tlsanl_wrapper import parse_skttls_summ, run_tlsanl


def init_logger(pdb_id, verbose):
    log_name = pdb_id + ".log"
    log_file_path = os.path.join(pyconfig.get("BDB_FILE_DIR_PATH"), log_name)

    fmt = "%(asctime)s | %(levelname)-7s | {0:4s} | %(message)s".format(pdb_id)

    logging.basicConfig(
        filename=log_file_path,
        filemode="w",
        level=logging.INFO if not verbose else logging.DEBUG, format=fmt)


def date_handler(obj):
    return obj.date().isoformat() if hasattr(obj, 'isoformat') else obj


def create_bdb_entry(pdb_file_path, pdb_id, verbose=False):
    """Create a bdb entry.

    Return True when a bdb has been created successfully.
    """

    _log.debug("Creating bdb entry...")

    # Parse the given pdb file into a dict...
    pdb_records = parse_pdb_file(pdb_file_path)

    # and a Biopython structure
    structure = get_structure(pdb_file_path, pdb_id, verbose)

    bdbd = {"pdb_id": pdb_id}
    expdta = check_exp_methods(pdb_records, pdb_id)
    bdbd.update(expdta)
    created_bdb_file = False
    if expdta["expdta_useful"]:
        refi_data = get_refi_data(pdb_records, structure, pdb_id)
        bdbd.update(refi_data)

        # Info about B-factor group type
        b_group = determine_b_group(structure)
        bdbd.update(b_group)

        # skttles outliers
        skttls = {"skttls_tot": None,
                  "skttls_95th": None,
                  "skttls_99th": None}
        bdbd.update(skttls)

        if refi_data["is_bdb_includable"]:
            bdb_file_dir = pyconfig.get("BDB_FILE_DIR_PATH")
            bdb_file_path = os.path.join(bdb_file_dir, pdb_id + ".bdb")
            if refi_data["req_tlsanl"]:
                if run_tlsanl(
                        pdb_file_path=pdb_file_path,
                        xyzout=bdb_file_path,
                        pdb_id=pdb_id,
                        log_out_dir=bdb_file_dir):
                    created_bdb_file = True
                    tlsanl_log = os.path.join(bdb_file_dir,
                                              pyconfig.get("TLSANL_LOG"))
                    skttls = parse_skttls_summ(tlsanl_log=tlsanl_log)
                    bdbd.update(skttls)

            elif refi_data["b_msqav"]:
                if write_multiplied_8pipi(
                        pdb_file_path=pdb_file_path,
                        xyzout=bdb_file_path,
                        pdb_id=pdb_id,
                        verbose=verbose):
                    created_bdb_file = True

            elif refi_data["assume_iso"]:
                shutil.copy(pdb_file_path, bdb_file_path)
                created_bdb_file = True

            else:
                message = "Unexpected bdb status"
                write_whynot(pdb_id, message)
                _log.error("{}.".format(message))

        # Write the bdb metadata to a json file
        try:
            with open(os.path.join(pyconfig.get("BDB_FILE_DIR_PATH"),
                      pdb_id + ".json"),
                      "w") as f:
                json.dump(bdbd, f, sort_keys=True, indent=4,
                          default=date_handler)
        except IOError as ex:
            _log.error(ex)
            return False

    return created_bdb_file


def main():
    """Create a bdb entry."""

    parser = argparse.ArgumentParser(
        description="Create a BDB entry. The BDB entry (.bdb) or WHY NOT\
        (.whynot) file, log and json files will be created in a separate\
        directory using the given pdb_id and the directory structure:\
        BDB_ROOT/ab/1abc/1abc.(bdb|whynot|log|json)")
    parser.add_argument(
        "-v", "--verbose",
        help="show verbose output",
        action="store_true")
    parser.add_argument(
        "bdb_root_path",
        help="Root directory of the bdb data.",
        type=lambda x: is_valid_directory(parser, x))
    parser.add_argument(
        "pdb_file_path",
        help="PDB file location.",
        type=lambda x: is_valid_file(parser, x))
    parser.add_argument(
        "pdb_id",
        help="PDB accession code.",
        type=lambda x: is_valid_pdbid(parser, x))
    args = parser.parse_args()

    pyconfig.set("BDB_FILE_DIR_PATH", get_bdb_entry_outdir(args.bdb_root_path,
                                                           args.pdb_id))
    init_logger(args.pdb_id, args.verbose)

    # Check that the system has the required programs and libraries installed
    check_deps()

    if create_bdb_entry(pdb_file_path=args.pdb_file_path, pdb_id=args.pdb_id,
                        verbose=args.verbose):
        _log.debug("Finished bdb entry.")
    # exit with status 0 when a BDB or a WHY NOT entry has been created
