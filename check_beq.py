#!/usr/bin/env python
from __future__ import division, print_function
from bdb_utils import get_pdb_header_and_trailer, init_bdb_logger,\
        write_whynot, PDB_LOGFORMAT
import argparse
import Bio.PDB
import itertools
import logging
import os
import re
import subprocess
import sys

# Configure logging
_log = logging.getLogger("bdb")

UF = 8*3.14159265359**2

def check_beq(pdb_xyz, pdb_id=None, verbose=False):
    """Determine if Beq values are the same as the reported B-factors.

    The margin is 0.015 Angstrom**2

    Return a dictionary with values that are None if ANISOU records are absent.
    beq_identical: float that indicates the percentage of B-factors in the
                   ATOM records that could be reproduced within the margin by
                   calculating the Beq values from the ANISOU records.
    correct_uij  : False if a non-standard combination of the Uij values in the
                   ANISOU records was necessary to reproduce the B-factors.
    """
    pdb_id = pdb_xyz if pdb_id is None else pdb_id
    margin = 0.015
    p = Bio.PDB.PDBParser(QUIET=not verbose)
    structure = p.get_structure(pdb_id, pdb_xyz)
    atoms = structure.get_atoms()
    has_anisou = False
    eq = 0
    ne = 0
    reproduced = 0.0
    correct_uij = True
    for atom in structure.get_atoms():
        if atom is not None:
            anisou = atom.get_anisou()
            if anisou is not None:
                has_anisou = True
                # Beq = 8*pi**2*Ueq
                # Ueq = 1/3<u.u> == 1/3<|u|**2> = 1/3(U11+U22+U33)
                beq = UF * sum(anisou[0:3]) / 3
                b   = atom.get_bfactor()
                if (b - margin) <= beq <= (b + margin):
                    eq = eq + 1
                elif check_combinations(anisou, b, margin, pdb_id):
                    """ e.g. 2a83, 2p6e, 2qik, 3bik, 3d95, 3d96, 3g5t
                    """
                    eq = eq + 1
                    correct_uij = False
                    _log.debug(("{0:" + PDB_LOGFORMAT + "} | "\
                                "B-factor reproduced by non-standard "\
                                "combination of Uij values in the ANISOU "\
                                "record of ATOM: {1:s}").
                                format(pdb_id, atom.get_full_id()))
                else:
                    """ e.g 1g8t, 1kr7, 1llr, 1mgr, 1o9g, 1pm1, 1q7l, 1qjp,
                    1s2p, 1si6, 1sxu, 1sxy, 1sy0, 1sy2, 1ug6, 1x9q, 2a83, 2acp,
                    2at5, 2bwi, 2ceu, 2fri, 2frj, 2hmn, 2htx, 2j73, 2p6e, 2p6f,
                    2p6g, 2qfn, 2qik, 2v0a, 2xgb, 2xl6, 2xle, 2xlw, 3bwo, 3dqy,
                    3fde, 3g5t, 3jql, 3nju, 3nna, 3oxp
                    """
                    ne = ne + 1
                    _log.debug(("{0:" + PDB_LOGFORMAT + "} | "\
                                "Beq not identical to B-factor in ATOM "\
                                "record: {1:s} {2:3.2f} {3:3.2f}").format(
                                    pdb_id,
                                    atom.get_full_id(),
                                    b,
                                    beq))
        else:
            break;
    reproduced = eq / (eq + ne) if has_anisou else None
    correct_uij = correct_uij if has_anisou else None
    return {"beq_identical": reproduced,
            "correct_uij"  : correct_uij}

def check_combinations(anisou, b, margin, pdb_id=None):
    """Check if the B-factor can be reproduced by non-standard U combinations.

    Standard: U11, U22, and U33 are the first three values in the ANISOU record
    """
    pdb_id = pdb_id if pdb_id else ""
    assert(len(anisou) == 6)
    reproduced = False
    for c in itertools.combinations(list(xrange(0, 6)), 3):
        if c == (0, 1, 2): # we have already calculated this
            pass
        beq = UF * (anisou[c[0]] + anisou[c[1]] + anisou[c[2]])/3
        if (b - margin) <= beq <= (b + margin):
            reproduced = True
            _log.debug(("{0:" + PDB_LOGFORMAT + "} | B-factor could only be "\
                       "reproduced by combining non-standard Uij values "\
                       "{1:d} {2:d} {3:d}.").format(pdb_id, c[0], c[1], c[2]))
            break
    return reproduced

def is_backbone(atom):
    """Return True if the atom is probably a backbone atom."""
    return atom.get_name() in ['N', 'CA', 'C', 'O']

def b_equal(b1, b2):
    """Return True if the two B-factor values are equal within a margin.

    Margin 0.02 Angstrom**2
    """
    margin = 0.02
    return (b1 - margin) <= b2 <= (b1 + margin)

def determine_b_group(pdb_xyz, pdb_id=None, verbose=False):
    """Determine the most likely B-factor parameterization.

    Output can be one of the strings
    overall (margin 0.02 Angstrom**2)
    residue
    two-back-and-side
    individual
    """
    # TODO take more atoms and residues into account,
    #      current approach is rather greedy
    # TODO filter out zero occupancy? (e.g. 1etu)
    # TODO also check on chain level?
    group = "individual"
    margin = 0.02
    pdb_id = pdb_xyz if pdb_id is None else pdb_id
    p = Bio.PDB.PDBParser(QUIET=not verbose)
    structure = p.get_structure(pdb_id, pdb_xyz)
    residues = structure.get_residues()
    b_res = list()
    i = 0
    max_res = 2
    # 2 useful residues should be sufficient to make a decision
    while (i < max_res):
        res = residues.next()
        if res.get_id()[0] == " ": # Exclude HETATM and waters
            b_back = list()
            b_side = list()
            for atom in res:
                if not re.match("H", atom.get_name()): # Exclude hydrogens
                    b = atom.get_bfactor()
                    _log.debug(("{0:" + PDB_LOGFORMAT + "} | {1:s} - B-factor"\
                                ": {2:3.2f}").format(
                                    pdb_id,
                                    atom.get_full_id(),
                                    b,
                                    ))
                    if is_backbone(atom):
                        b_back.append(atom.get_bfactor())
                        if len(b_back) > 1: # Whithin backbone
                            if not b_equal(b_back[0], b_back[1]):
                                #group = "individual"
                                i = max_res # stop comparing residues...
                                break; # or atoms
                    else:
                        b_side.append(atom.get_bfactor())
                        if len(b_side) > 1: # Whithin side-chain
                            if not b_equal(b_side[0], b_side[1]):
                                #group = "individual"
                                i = max_res
                                break;
                    if len(b_back) > 0 and len(b_side) > 0:
                        # Between backbone and side-chain
                        if not b_equal(b_back[0], b_side[0]):
                            group = "two-back-and-side"
                            i = max_res
                            break;
                        else:
                            group = "residue"
                            break;
            #b_back.extend(b_side)
            b_res.append(b_back)
            i = i + 1 # useful atoms in this residue
    if len(b_res) > 1 and b_equal(b_res[0][0],b_res[1][0]):
        group = "overall"
    _log.info(("{0:" + PDB_LOGFORMAT + "} | Most likely B-factor group type:"\
               " {1:s}.").format(pdb_id, group))
    return group

def get_check_beq_parser():
    """Create and return an argument parser."""
    parser = argparse.ArgumentParser(
            description="Check Beq, find B-factor model or "
                        "calculate B-factors"
            )
    parser.add_argument(
            "-v",
            "--verbose",
            help="show versbose output",
            action="store_true"
            )
    parser.add_argument(
            "--pdbid",
            help="PDB file name."
            )
    sub = parser.add_subparsers(
            help="sub-command help"
            )
    xyz_help = "Input coordinates in PDB format."
    calc = sub.add_parser(
            "calc",
            help="Multiply B-factor column with 8*pi**2"
            )
    calc.add_argument(
            "xyzin",
            help=xyz_help
            )
    calc.add_argument(
            "xyzout",
            help="Output coordinates in PDB format."
            )
    check = sub.add_parser(
            "check",
            help="Check Beq values or find B-factor model."
            )
    check.add_argument(
            "xyzin",
            help=xyz_help
            )
    check_mode = check.add_mutually_exclusive_group(required=True)
    check_mode.add_argument(
            "--beq",
            help="Check if Beq values calculated from the ANISOU records in a "
            "PDB file are equal to the B-factor values reported in the ATOM "
            "records.",
            action="store_true"
            )
    check_mode.add_argument(
            "--group",
            help="Determine most likely B-factor model parameterization "
            "(overall, one per residue, two per residue (backbone and "
            "side-chain, or individual). TLS groups are not "
            "taken into account.",
            action="store_true"
            )
    return parser

def greet(pdb_id, mode=None):
    """Say hello."""
    if mode == "beq":
        _log.info(("{0:" + PDB_LOGFORMAT + "} | "\
                   "Checking Beq values in ANISOU records...").format(pdb_id))
    elif mode == "group":
        _log.info(("{0:" + PDB_LOGFORMAT + "} | "\
                   "Determining most likely B-factor model...").format(pdb_id))
    elif mode == "calc":
        _log.info(("{0:" + PDB_LOGFORMAT + "} | "\
                   "Calculating B-factors from Uiso values...").format(pdb_id))

def multiply(structure):
    """Multiply B-factors with 8*pi**2."""
    for atom in structure.get_atoms():
        atom.set_bfactor(UF * atom.get_bfactor())
    return structure

def transfer_header_and_trailer(xyzin, xyzout):
    """Transfer header and trailer from xyzin to xyzout."""
    h, t = get_pdb_header_and_trailer(xyzin)
    records = list()
    records.extend(h)
    end = "END"
    try:
        with open(xyzout, "r") as pdb_out:
            for coord in pdb_out:
                if re.search(r"^END\s*$", coord):
                    end = coord.rstrip("\n")
                else:
                    records.append(coord.rstrip("\n"))
        records.extend(t)
        records.append(end)
        with open(xyzout + "2", "w") as pdb_out:
            for record in records:
                pdb_out.write("{0:s}\n".format(record))
    except IOError as ex:
        _log.error(ex)

def write_multiplied(xyzin, xyzout, pdb_id=None, verbose=False):
    """Multiply the B-factors in the input PDB file with 8*pi^2."""
    pdb_id = pdb_id if pdb_id else "usio"
    greet(pdb_id, mode="calc")
    p = Bio.PDB.PDBParser(QUIET=not verbose)
    structure = p.get_structure(pdb_id, xyzin)
    structure = multiply(structure)
    io = Bio.PDB.PDBIO()
    io.set_structure(structure)
    # Header and trailer records not present in this output file
    io.save(xyzout)
    transfer_header_and_trailer(xyzin, xyzout)

def report_beq(pdb_id, reproduced):
    """Report if Beqs are identical to B-factors."""
    if reproduced["beq_identical"] is None:
        _log.debug(("{0:" + PDB_LOGFORMAT + "} | No ANISOU records").
                format(pdb_id))
        return
    if not reproduced["correct_uij"]:
        _log.warn(("{0:" + PDB_LOGFORMAT + "} | One or more B-factors could "\
                   "only be reproduced by a non-standard combination of Uij "\
                   "values in the corresponding ANISOU record.").
                   format(pdb_id))
    if reproduced["beq_identical"] == 1:
        _log.info(("{0:" + PDB_LOGFORMAT + "} | The B-factors in the"\
                " ATOM records could all be reproduced within 0.015 A**2 by "\
                "calculating Beq from the corresponding ANISOU records.").
                format(pdb_id,100*(1-reproduced["beq_identical"])))
    elif reproduced["beq_identical"] < 1:
        _log.warn(("{0:" + PDB_LOGFORMAT + "} | {1:3.2f}% of the B-factors in"\
                " the ATOM records could not be reproduced within 0.015 A**2 "\
                "by calculating Beq from the corresponding ANISOU records.").
                format(pdb_id,100*(1-reproduced["beq_identical"])))
    else:
        _log.info(("{0:" + PDB_LOGFORMAT + "} | "\
                "No ANISOU records.").format(pdb_id))

def swear(pdb_id, message):
    """Report problems."""
    write_whynot(pdb_id, message)
    _log.error(("{0:" + PDB_LOGFORMAT + "} | {1:s}").format(pdb_id, message))

if __name__ == "__main__":
    """Run Beq check or multiply Uiso with 8*pi^2."""
    parser = get_check_beq_parser()
    args = parser.parse_args()
    print(args)
    pdb_id = args.pdbid if args.pdbid is not None else args.xyzin
    _log = init_bdb_logger(pdb_id, global_log=True)
    if args.verbose:
        _log.setLevel(logging.DEBUG)
    if args.beq:
        # Check Beq mode
        greet(pdb_id, mode="beq")
        report_beq(pdb_id, check_beq(args.xyzin, pdb_id))
    elif args.group:
        # Check group mode
        greet(pdb_id, mode="group")
        determine_b_group(args.xyzin, pdb_id, args.verbose)
    else:
        # Calc mode
        write_multiplied(args.xyzin, args.xyzout, pdb_id, args.verbose)
