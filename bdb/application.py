from __future__ import print_function

import argparse
import json
import logging
import os
import shutil
import sys

from bdb.bdb_utils import (is_valid_directory, is_valid_file, is_valid_pdbid,
                           get_bdb_entry_outdir, write_whynot)
from bdb.check_beq import write_multiplied
from bdb.expdta import check_exp_methods
from bdb.pdb.parser import parse_pdb_file
from bdb.refprog import do_refprog
from bdb.tlsanl_wrapper import run_tlsanl


_log = logging.getLogger(__name__)


def init_logger(pdb_id, root, verbose):
    log_name = pdb_id + ".log"
    out_dir = get_bdb_entry_outdir(root, pdb_id)
    log_file_path = os.path.join(out_dir, log_name)

    fmt = "%(asctime)s | %(levelname)-7s | {0:4s} | %(message)s".format(pdb_id)

    logging.basicConfig(
        filename=log_file_path,
        filemode='w',
        level=logging.INFO if not verbose else logging.DEBUG,
        format=fmt)


def do_bdb(bdb_root_path, pdb_file_path, pdb_id, global_files):
    """Create a bdb entry.

    Return a Boolean."""

    _log.debug("Creating bdb entry...")

    # Parse the given pdb file into a dict.
    pdb_records = parse_pdb_file(pdb_file_path)

    out_dir = get_bdb_entry_outdir(bdb_root_path, pdb_id)
    bdbd = {"pdb_id": pdb_id}
    expdta = check_exp_methods(pdb_records, pdb_id, out_dir, global_files)
    bdbd.update(expdta)
    done = False
    if expdta["expdta_useful"]:
        refprog = do_refprog(pdb_file_path, pdb_id, out_dir, global_files)
        bdbd.update(refprog)

        # Write the bdb metadata to a json file
        try:
            with open(os.path.join(out_dir, pdb_id + ".json"), "w") as f:
                json.dump(bdbd, f, sort_keys=True, indent=4)
        except IOError as ex:
            _log.error(ex)
            return False

        if refprog["refprog_useful"]:
            bdb_file_path = os.path.join(out_dir, pdb_id + ".bdb")
            # TODO do we need extractor? or tlsextract (ccp4)
            if refprog["req_tlsanl"]:
                if run_tlsanl(
                        xyzin=pdb_file_path,
                        xyzout=bdb_file_path,
                        pdb_id=pdb_id,
                        log_out_dir=out_dir
                        ):
                    done = True
            elif refprog["b_msqav"]:
                if write_multiplied(
                        xyzin=pdb_file_path,
                        xyzout=bdb_file_path,
                        pdb_id=pdb_id
                        ):
                    done = True
            elif refprog["assume_iso"]:
                shutil.copy(pdb_file_path, bdb_file_path)
                done = True
            else:
                message = "Unexpected bdb status"
                write_whynot(pdb_id, message, directory=out_dir)
                _log.error("{}.".format(message))
    return done


def main():
    """Create a bdb entry.

    The bdb entry (.bdb) or WHY NOT (.whynot) file,
    log and json files will be created in a separate directory
    using the given pdbid and the directory structure:
    BDB_ROOT/ab/1abc/1abc.(bdb|whynot|log|json)
    """

    parser = argparse.ArgumentParser(
        description="Create a bdb entry")
    parser.add_argument(
        "-g", "--global_files",
        help="Create files with PDB-wide information. Useful for local bdb "
             "copies. WARNING: do not use in an embarassingly parallel "
             "setting!",
        action="store_true")
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

    init_logger(args.pdb_id, args.bdb_root_path, args.verbose)

    # Check that the system has the required programs and libraries installed
    # TODO: This should be moved to a `setup.py` file.
    import requirements

    if do_bdb(args.bdb_root_path,
              args.pdb_file_path,
              args.pdb_id, args.global_files):
        _log.debug("Finished bdb entry.")
        sys.exit(0)
    else:
        sys.exit(1)
