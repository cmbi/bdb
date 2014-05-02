from __future__ import print_function

import argparse
import json
import logging
import os
import pyconfig
import shutil
import sys

from bdb.bdb_utils import (is_valid_directory, is_valid_file, is_valid_pdbid,
                           get_bdb_entry_outdir, write_whynot)
from bdb.check_beq import (determine_b_group, get_structure,
                           write_multiplied_8pipi)
from bdb.expdta import check_exp_methods
from bdb.pdb.parser import parse_pdb_file
from bdb.refprog import get_refi_data
from bdb.tlsanl_wrapper import run_tlsanl


_log = logging.getLogger(__name__)

def init_logger(pdb_id, verbose):
    log_name = pdb_id + ".log"
    log_file_path = os.path.join(pyconfig.get("BDB_FILE_DIR_PATH"), log_name)

    fmt = "%(asctime)s | %(levelname)-7s | {0:4s} | %(message)s".format(pdb_id)

    logging.basicConfig(
        filename=log_file_path,
        filemode="w",
        level=logging.INFO if not verbose else logging.DEBUG, format=fmt)


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
        b_group = determine_b_group(structure, pdb_id)
        bdbd.update(b_group)


        # Write the bdb metadata to a json file
        try:
            with open(os.path.join(pyconfig.get("BDB_FILE_DIR_PATH"),
                pdb_id + ".json"),
                   "w") as f:
                json.dump(bdbd, f, sort_keys=True, indent=4)
        except IOError as ex:
            _log.error(ex)
            return False

        if refi_data["is_bdb_includable"]:
            bdb_file_path = os.path.join(pyconfig.get("BDB_FILE_DIR_PATH"),
                    pdb_id + ".bdb")
            if refi_data["req_tlsanl"]:
                if run_tlsanl(
                        pdb_file_path=pdb_file_path,
                        xyzout=bdb_file_path,
                        pdb_id=pdb_id,
                        log_out_dir=pyconfig.get("BDB_FILE_DIR_PATH")):
                    created_bdb_file = True
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
    return created_bdb_file


def main():
    """Create a bdb entry."""

    parser = argparse.ArgumentParser(
        description="Create a bdb entry. The bdb entry (.bdb) or WHY NOT\
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

    pyconfig.set("BDB_FILE_DIR_PATH",get_bdb_entry_outdir(args.bdb_root_path,
                                                          args.pdb_id))
    init_logger(args.pdb_id, args.verbose)

    # Check that the system has the required programs and libraries installed
    # TODO: This should be moved to the `setup.py` file or at least be provided
    #       via a function.
    import requirements

    if create_bdb_entry(pdb_file_path=args.pdb_file_path, pdb_id=args.pdb_id,
            verbose=args.verbose):
        _log.debug("Finished bdb entry.")
        sys.exit(0)
    else:
        sys.exit(1)
