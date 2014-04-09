#!/usr/bin/env python
from __future__ import division

import argparse
import itertools
import logging
import os
import re
import shutil
import subprocess
import sys

import Bio.PDB
import numpy

from bdb.bdb_utils import (get_pdb_header_and_trailer, init_bdb_logger,
                           write_whynot)

# Configure logging
_log = logging.getLogger("bdb")

UF = 8*3.14159265359**2

def check_beq(pdb_file_path, pdb_id=None, verbose=False):
    """Determine if Beq values are the same as the reported B-factors.

    The margin is 0.015 Angstrom**2

    Return a dictionary with values that are None if ANISOU records are absent.
    beq_identical: float that indicates the percentage of B-factors in the
                   ATOM records that could be reproduced within the margin by
                   calculating the Beq values from the ANISOU records.
    correct_uij  : False if a non-standard combination of the Uij values in the
                   ANISOU records was necessary to reproduce the B-factors.
    """
    pdb_id = pdb_file_path if pdb_id is None else pdb_id
    _log.info(("{0:4s} | Checking Beq values in ANISOU records...").format(pdb_id))
    margin = 0.015
    p = Bio.PDB.PDBParser(QUIET=not verbose)
    structure = p.get_structure(pdb_id, pdb_file_path)
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
                if numpy.isclose(b, beq, atol=margin):
                    eq = eq + 1
                elif check_combinations(anisou, b, margin, pdb_id):
                    """ e.g. 2a83, 2p6e, 2qik, 3bik, 3d95, 3d96, 3g5t
                    """
                    eq = eq + 1
                    correct_uij = False
                    _log.debug(("{0:4s} | B-factor reproduced by "\
                                "non-standard combination of Uij values in " \
                                "the ANISOU record of ATOM: {1:s}").format(
                                    pdb_id, atom.get_full_id()))
                else:
                    """ e.g 1g8t, 1kr7, 1llr, 1mgr, 1o9g, 1pm1, 1q7l, 1qjp,
                    1s2p, 1si6, 1sxu, 1sxy, 1sy0, 1sy2, 1ug6, 1x9q, 2a83, 2acp,
                    2at5, 2bwi, 2ceu, 2fri, 2frj, 2hmn, 2htx, 2j73, 2p6e, 2p6f,
                    2p6g, 2qfn, 2qik, 2v0a, 2xgb, 2xl6, 2xle, 2xlw, 3bwo, 3dqy,
                    3fde, 3g5t, 3jql, 3nju, 3nna, 3oxp
                    """
                    ne = ne + 1
                    _log.debug(("{0:4s} | Beq not identical to B-factor in " \
                                "ATOM record: {1:s} {2:3.2f} {3:3.2f}").format(
                                    pdb_id, atom.get_full_id(), b,beq))
        else:
            # TODO: Why break if an ATOM is None? This *could* happen at the
            #       very first. If so, isn't this a special case?
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
        if numpy.isclose(b, beq, atol=margin):
            reproduced = True
            _log.debug(("{0:4s} | B-factor could only be reproduced by " \
                        "combining non-standard Uij values "\
                        "{1:d} {2:d} {3:d}.").format(pdb_id, c[0], c[1], c[2]))
            break
    return reproduced

def determine_b_group(pdb_file_path, pdb_id=None, verbose=False):
    """Determine the most likely B-factor parameterization.

    Return a dictionary with separated output for protein and nucleic acid and
    a Boolean that indicates if the structure is a calpha trace.

    output can be one of the strings
    overall           e.g. 1etu
    residue_1ADP      e.g. the protein in 1hlz
    residue_2ADP      e.g. the DNA in 1hlz
    individual        most PDB files
    no_b-factors      e.g. 1mcb, 3zxa, 2yhx

    or None if protein or nucleic acid are not present.

    (margin 0.01 Angstrom**2)

    Warning: currently only the first protein and/or nucleid acid chains
    encountered are taken into account. The same parameterization is assumed
    for other chains (if they exist).
    """
    group = {
            "protein_b": None,
            "nucleic_b": None,
            "calpha_only": False,
            }
    margin = 0.01
    pdb_id = pdb_file_path if pdb_id is None else pdb_id
    _log.info(("{0:4s} | Determining most likely B-factor group type").format(
        pdb_id))
    structure = None
    try:
        p = Bio.PDB.PDBParser(QUIET=not verbose)
        structure = p.get_structure(pdb_id, pdb_file_path)
    except (AttributeError, IndexError, ValueError, AssertionError,
            Bio.PDB.PDBExceptions.PDBConstructionException):
        _log.error(("{0:4s} | Biopython Error.").format(pdb_id))
    if structure is not None:
        chains = structure.get_chains()
        for c in chains:
            if is_protein_chain(c):
                if group["protein_b"] is None:
                    group["protein_b"] = determine_b_group_chain(c)
            elif is_nucleic_chain(c):
                if group["nucleic_b"] is None:
                    group["nucleic_b"] = determine_b_group_chain(c)
            elif is_calpha_trace(c):
                if group["protein_b"] is None:
                    group["calpha_only"] = True
                    _log.info(("{0:4s} | Calpha-only chain(s) present").format(
                        pdb_id))
                    group["protein_b"] = determine_b_group_chain(c)
            else:
                _log.error(("{0:4s} | Chain {1:s}: no protein or nucleic " \
                            "acid chain found.").format(pdb_id, c.get_id()))
        _log.info(("{0:4s} | Most likely B-factor group type protein: {1:s}" \
                   "| nucleic acid: {2:s}.").format(pdb_id,
                    group["protein_b"] if group["protein_b"] is not None else\
                            "not present",
                    group["nucleic_b"] if group["nucleic_b"] is not None else\
                            "not present",))
    return group

def determine_b_group_chain(chain):
    """Return the most likely B-factor group type for this chain.

    Return a string

    overall           e.g. 1etu
    residue_1ADP      e.g. the protein in 1hlz
    residue_2ADP      e.g. the DNA in 1hlz
    individual        most PDB files
    (margin 0.01 Angstrom**2)

    Warning: the current approach is rather greedy as only the first four
    residues of the chain are taken into account. A uniform parameterization
    accross the chain is assumed. If multiple domains with different overall
    B-factors are present in the same chain, the ouput will still be overall.

    Note: if only the first three residues would have been considered,
    the approach would have been too greedy for 1hlz chain B
    """
    margin = 0.01
    residues = chain.get_residues()
    group = "individual"
    b_res = list()
    i = 0
    max_res = 4
    # 4 useful residues should be sufficient to make a decision
    while (i < max_res):
        try:
            res = residues.next()
        except StopIteration:
            # e.g. 1c0q
            _log.warn(("{0:4s} | Chain {1:s} has less than {2:d} useful " \
                       "residues composed of ATOMs.").format(
                           chain.get_full_id()[0], chain.get_id(), max_res))
            break;
        if res.get_id()[0] == " ": # Exclude HETATM and waters
            b_atom = list()
            for atom in res:
                # Exclude hydrogens and zero occupancy (many in e.g. 1etu)
                if not re.match("H", atom.get_name())\
                        and atom.get_occupancy() > 0:
                    b = atom.get_bfactor()
                    _log.debug(("{0:4s} | {1:s} - B-factor: {2:3.2f}").format(
                                    chain.get_full_id()[0], atom.get_full_id(),
                                    b,))
                    b_atom.append(b)
            # Useful atoms in this residue
            if len(b_atom) > 0:
                b_res.append(b_atom)
                i = i + 1
            # Determine the B-factor type for this residue if it is not CA-only
            if len(b_atom) > 1:
                b_atom = sorted(b_atom)
                if numpy.isclose(b_atom[-1], b_atom[0], atol=margin):
                    group = "residue_1ADP"
                elif len(b_atom) > 3 and \
                        numpy.allclose(
                                [b_atom[-1], b_atom[1]],
                                [b_atom[-2], b_atom[0]],
                                atol=margin,
                                ) and \
                        not numpy.isclose(
                                b_atom[-2],
                                b_atom[1],
                                atol=margin,
                                ):
                    group = "residue_2ADP"
                else:
                    group = "individual"
    if len(b_res) > max_res - 1 and numpy.isclose(
            b_res[0][0],
            b_res[-1][0],
            atol=margin):
        if numpy.isclose(b_res[-1][-1], 0):
            group = "no_b-factors"
        else:
            group = "overall"
    return group

def determine_b_group_chain_greedy(chain):
    """Return the most likely B-factor group type for this chain.

    Return a string

    Warning: the current approach is rather greedy as only the first four
    residues of the chain are taken into account. A uniform parameterization
    accross the chain is assumed.

    Note: if only the first three residues would have been considered,
    the approach would have been too greedy for 1hlz chain B

    Warning: more atoms need to be compared. E.g. from 1kkl:
    ATOM      1  N   GLU A 135      49.349  63.473  12.155  1.00 80.44
    ATOM      2  CA  GLU A 135      50.791  63.168  12.411  1.00 80.45
    ATOM      3  C   GLU A 135      50.998  61.654  12.472  1.00 79.50
    ATOM      4  O   GLU A 135      50.236  60.900  11.876  1.00 79.31
    ATOM      5  CB  GLU A 135      51.257  63.833  13.722  1.00 49.05
    Individual B-factors, not 2 per residue
    """
    residues = chain.get_residues()
    group = "individual"
    b_res = list()
    i = 0
    max_res = 4
    # 4 useful residues should be sufficient to make a decision
    while (i < max_res):
        res = residues.next()
        if res.get_id()[0] == " ": # Exclude HETATM and waters
            b_back = list()
            b_side = list()
            for atom in res:
                # Exclude hydrogens and zero occupancy (many in e.g. 1etu)
                if not re.match("H", atom.get_name())\
                        and atom.get_occupancy() > 0:
                    b = atom.get_bfactor()
                    _log.debug(("{0:4s} | {1:s} - B-factor: {2:3.2f}").format(
                                    chain.get_full_id()[0], atom.get_full_id(),
                                    b,))
                    if is_heavy_backbone(atom):
                        b_back.append(atom.get_bfactor())
                        if len(b_back) > 1: # Whithin backbone
                            if not numpy.isclose(
                                    b_back[0],
                                    b_back[1],
                                    atol=margin
                                    ):
                                #group = "individual"
                                i = max_res # stop comparing residues...
                                break; # or atoms
                    else:
                        b_side.append(atom.get_bfactor())
                        if len(b_side) > 1: # Whithin side-chain
                            if not numpy.isclose(
                                    b_side[0],
                                    b_side[1],
                                    atol=margin
                                    ):
                                #group = "individual"
                                i = max_res
                                break;
                    if len(b_back) > 0 and len(b_side) > 0:
                        # Between backbone and side-chain
                        if not numpy.isclose(
                                b_back[0],
                                b_side[0],
                                atol=margin
                                ):
                            group = "residue_2ADP"
                            i = max_res
                            break;
                        else:
                            group = "residue_1ADP"
                            break;
            #b_back.extend(b_side)
            b_res.append(b_back)
            i = i + 1 # useful atoms in this residue
    if len(b_res) > max_res - 1 and numpy.isclose(
            b_res[0][0],
            b_res[-1][0],
            atol=margin
            ):
        group = "overall"
    return group

def has_amino_acid_backbone(residue):
    """Return True if the residue's backbone looks like protein."""
    for atom in ("N", "CA", "C", "O"):
        if not atom in residue.child_dict:
            return False
    return True

def has_sugar_phosphate_backbone(residue):
    """Return True if the residue's backbone looks like nucleic acid."""
    for atom in ("P", "OP1", "OP2", "O5'", "C5'", "C4'",
                 "O4'", "C3'","O3'", "C2'", "C1'"):
        if not atom in residue.child_dict:
            return False
    return True

def is_heavy_backbone(atom):
    """Return True if the atom looks like a backbone atom."""
    return atom.get_name() in [
            "N", "CA", "C", "O", # Protein
            "P", "OP1", "OP2", "O5'", "C5'", "C4'",
            "O4'", "C3'","O3'", "C2'", "O2'", "C1'", # DNA/RNA
            ]

def is_calpha_trace(chain):
    """Return True if more than 75% of the atoms in the chain are ca atoms.

    The function accounts for unexpected residues and atoms (such as UNK and
    hetatms listed as atms) by calculating the percentage of ca atoms.

    Example: 1efg chain A contains 6 protein domains (each with a different
    overall B-factor) and GDP, chain B and C are composed of UNK residues.
    """
    ca = []
    for atom in chain.get_atoms():
        if atom.get_name() == "CA":
            ca.append(1)
        else:
            ca.append(0)
    ca_ratio = numpy.count_nonzero(ca) / len(ca)
    return ca_ratio >= 0.75

def is_nucleic_chain(chain):
    """Return True if the first 10 residues of the chain look like nucleotides.

    It is assumed mixed protein and nucleic acid chains don't exist.
    Therefore this approach is rather greedy.
    """
    residues = chain.get_residues()
    check_max = 10
    residues_checked = 0
    residues.next() # The first residue does not contain the phosphate,
                    # we rather start checking from the second residue
    for res in residues:
        if residues_checked < check_max \
                and res.get_id()[0] == " ": # Exclude HETATM and waters
            if not has_sugar_phosphate_backbone(res):
                return False
            residues_checked = residues_checked + 1
    return True

def is_protein_chain(chain):
    """Return True if the first 10 residues of the chain look like amino acids.

    It is assumed mixed protein and nucleic acid chains don't exist.
    Therefore this approach is rather greedy.
    """
    residues = chain.get_residues()
    check_max = 10
    residues_checked = 0
    for res in residues:
        if residues_checked < check_max \
                and res.get_id()[0] == " ": # Exclude HETATM and waters
            if not has_amino_acid_backbone(res):
                return False
            residues_checked = residues_checked + 1
    return True

def multiply(structure):
    """Multiply B-factors with 8*pi**2."""
    for atom in structure.get_atoms():
        atom.set_bfactor(UF * atom.get_bfactor())
    return structure

def transfer_header_and_trailer(pdb_file_path, xyzout):
    """Transfer header and trailer from pdb_file_path to xyzout."""
    transferred = False
    h, t = get_pdb_header_and_trailer(pdb_file_path)
    records = list()
    # Start with the header...
    records.extend(h)
    end = "END"
    try:
        with open(xyzout, "r") as pdb_out:
            # ... then copy coordinates
            for coord in pdb_out:
                # ... remember but skip END now
                if re.search(r"^END\s*$", coord):
                    end = coord.rstrip("\n")
                else:
                    records.append(coord.rstrip("\n"))
        # ... then copy the trailer
        records.extend(t)
        # ... and finally END
        records.append(end)
        # write a new file
        with open(xyzout + "2", "w") as pdb_out:
            for record in records:
                pdb_out.write("{0:s}\n".format(record))
        # replace the old file with the new file
        shutil.move(xyzout + "2", xyzout)
        transferred = True
    except IOError as ex:
        _log.error(ex)
    return transferred

def write_multiplied(pdb_file_path, xyzout, pdb_id=None, verbose=False):
    """Multiply the B-factors in the input PDB file with 8*pi^2."""
    pdb_id = pdb_id if pdb_id else "usio"
    _log.info(("{0:4s} | Calculating B-factors from Uiso values...").format(
        pdb_id))
    p = Bio.PDB.PDBParser(QUIET=not verbose)
    structure = p.get_structure(pdb_id, pdb_file_path)
    structure = multiply(structure)
    io = Bio.PDB.PDBIO()
    io.set_structure(structure)
    # Header and trailer records not present in this output file
    io.save(xyzout)
    return transfer_header_and_trailer(pdb_file_path, xyzout)

def report_beq(pdb_id, reproduced):
    """Report if Beqs are identical to B-factors."""
    if reproduced["beq_identical"] is None:
        _log.debug(("{0:4s} | No ANISOU records").format(pdb_id))
        return
    if not reproduced["correct_uij"]:
        _log.warn(("{0:4s} | One or more B-factors could only be reproduced " \
                   "by a non-standard combination of Uij values in the " \
                   "corresponding ANISOU record.").format(pdb_id))
    if reproduced["beq_identical"] == 1:
        _log.info(("{0:4s} | The B-factors in the ATOM records could all be " \
                   "reproduced within 0.015 A**2 by calculating Beq from " \
                   "the corresponding ANISOU records.").format(
                       pdb_id, 100 * (1 - reproduced["beq_identical"])))
    elif reproduced["beq_identical"] < 1:
        _log.warn(("{0:4s} | {1:3.2f}% of the B-factors in the ATOM records " \
                "could not be reproduced within 0.015 A**2 by calculating " \
                "Beq from the corresponding ANISOU records.").format(
                    pdb_id, 100 * (1 - reproduced["beq_identical"])))
    else:
        _log.info(("{0:4s} | No ANISOU records.").format(pdb_id))

if __name__ == "__main__":
    """Run Beq check or multiply Uiso with 8*pi^2."""
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
    calc = sub.add_parser(
            "calc",
            help="Multiply B-factor column with 8*pi**2"
            )
    calc.add_argument(
            "pdb_file_path",
            help="PDB file location."
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
            "pdb_file_path",
            help="PDB file location."
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
    args = parser.parse_args()
    pdb_id = args.pdbid if args.pdbid is not None else args.pdb_file_path
    _log = init_bdb_logger(pdb_id, global_log=True)
    if args.verbose:
        _log.setLevel(logging.DEBUG)
    if args.beq:
        # Check Beq mode
        report_beq(pdb_id, check_beq(args.pdb_file_path, pdb_id))
    elif args.group:
        # Check group mode
        determine_b_group(args.pdb_file_path, pdb_id, args.verbose)
    else:
        # Calc mode
        write_multiplied(args.pdb_file_path, args.xyzout, pdb_id, args.verbose)
