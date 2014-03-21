#!/usr/bin/env python
from __future__ import print_function
from bdb_utils import is_valid_directory, is_valid_file, is_valid_pdbid,\
get_bdb_entry_outdir, init_bdb_logger, write_dict_json, write_whynot,\
PDB_LOGFORMAT
from check_beq import write_multiplied
from expdta import do_expdta
from refprog import do_refprog
from tlsanl_wrapper import run_tlsanl
import argparse
import logging
import os
import re
import shutil
import sys

def do_bdb(root, pdb_xyz, pdb_id, global_files):
    """Create a bdb entry.

    Return a Boolean."""
    done = False
    greet(pdb_id)
    out_dir = get_bdb_entry_outdir(root, pdb_id)
    bdbd = {"pdb_id": pdb_id}
    expdta = do_expdta(pdb_xyz, pdb_id, out_dir, global_files)
    bdbd.update(expdta)
    if expdta["expdta_useful"]:
        refprog = do_refprog(pdb_xyz, pdb_id, out_dir, global_files)
        bdbd.update(refprog)
        write_dict_json(bdbd, os.path.join(out_dir, pdb_id + ".json"),
                pretty=True)
        if refprog["refprog_useful"]:
            bdb_xyz = os.path.join(out_dir, pdb_id + ".bdb")
            # TODO do we need extractor? or tlsextract (ccp4)
            if refprog["req_tlsanl"]:
                if run_tlsanl(
                        xyzin=pdb_xyz,
                        xyzout=bdb_xyz,
                        pdb_id=pdb_id,
                        log_out_dir=out_dir
                        ):
                    done = True
            elif refprog["b_msqav"]:
                if write_multiplied(
                        xyzin=pdb_xyz,
                        xyzout=bdb_xyz,
                        pdb_id=pdb_id
                        ):
                    done =True
            elif refprog["assume_iso"]:
                shutil.copy(pdb_xyz, bdb_xyz)
                done = True
            else:
                message = "Unexpected bdb status"
                write_whynot(pdb_id, message, directory=out_dir)
                _log.error(("{0:" + PDB_LOGFORMAT + "} | {1:s}.").
                        format(pdb_id, message))
    return done

def get_argparser():
    """Create and return an argument parser."""
    parser = argparse.ArgumentParser(
        description="Create a bdb entry")
    parser.add_argument(
        "-g", "--global_files",
        help="Create files with PDB-wide information. Useful for local bdb "\
             "copies. "\
             "WARNING: do not use in an embarassingly parallel setting!",
        action="store_true")
    parser.add_argument(
        "-v", "--verbose",
        help="show verbose output",
        action="store_true")
    parser.add_argument(
        "bdb_root",
        help="Root directory of the bdb data.",
        type=lambda x: is_valid_directory(parser, x))
    parser.add_argument(
        "xyzin",
        help="Input coordinates in PDB format.",
        type=lambda x: is_valid_file(parser, x))
    parser.add_argument(
        "pdbid",
        help="PDB file name.",
        type=lambda x: is_valid_pdbid(parser, x))
    return parser

def good_bye(pdb_xyz):
    """Say good_bye."""
    _log.debug(("{0:" + PDB_LOGFORMAT + "} | "\
                "Finished bdb entry.").format(pdb_xyz))

def greet(pdb_xyz):
    """Say hello."""
    _log.debug(("{0:" + PDB_LOGFORMAT + "} | "\
                "Creating bdb entry...").format(pdb_xyz))

if __name__ == "__main__":
    """Create a bdb entry.

    The bdb entry (.bdb) or WHY NOT (.whynot) file,
    log and json files will be created in a separate directory
    using the given pdbid and the directory structure:
    BDB_ROOT/ab/1abc/1abc.(bdb|whynot|log|json)
    """
    parser = get_argparser()
    args = parser.parse_args()
    _log = init_bdb_logger(args.pdbid, args.bdb_root)
    if args.verbose:
        _log.setLevel(logging.DEBUG)
    import requirements
    if do_bdb(args.bdb_root, args.xyzin, args.pdbid, args.global_files):
        good_bye(args.pdbid)
        sys.exit(0)
    else:
        sys.exit(1)
