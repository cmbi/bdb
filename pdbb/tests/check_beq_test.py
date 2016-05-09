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
from nose.tools import eq_, ok_, raises

from pdbb.check_beq import (check_beq, check_combinations, determine_b_group,
                            get_structure, is_calpha_trace, is_phos_trace,
                            has_amino_acid_backbone,
                            has_sugar_phosphate_backbone,
                            is_heavy_backbone, is_nucleic_chain,
                            is_protein_chain, multiply_bfactors_8pipi)

import numpy as np


def test_check_beq_identical():
    """Tests check_beq.

    e.g. 3zzw, 3zit, 2yjl, etc. etc. etc.
    """
    pdb_file_path = "pdbb/tests/pdb/files/3zzw.pdb"
    pdb_id = "3zzw"
    structure = get_structure(pdb_file_path, pdb_id)
    result = check_beq(structure)
    eq_(result["beq_identical"], 1.0)
    eq_(result["correct_uij"], True)


def test_check_beq_incorrect_uij():
    """Tests check_beq.

    e.g. 2a83, 2p6e, 2qik, 3bik, 3d95, 3d96, 3g5t, etc.
    """
    pdb_file_path = "pdbb/tests/pdb/files/2a83.pdb"
    pdb_id = "2a83"
    structure = get_structure(pdb_file_path, pdb_id)
    result = check_beq(structure)
    eq_(result["beq_identical"], 0.9419321685508736)
    eq_(result["correct_uij"], False)


def test_check_beq_not_identical():
    """Tests check_beq.

    e.g 1g8t, 1kr7, 1llr, 1mgr, 1o9g, 1pm1, 1q7l, 1qjp,
    1s2p, 1si6, 1sxu, 1sxy, 1sy0, 1sy2, 1ug6, 1x9q, 2a83, 2acp,
    2at5, 2bwi, 2ceu, 2fri, 2frj, 2hmn, 2htx, 2j73, 2p6e, 2p6f,
    2p6g, 2qfn, 2qik, 2v0a, 2xgb, 2xl6, 2xle, 2xlw, 3bwo, 3dqy,
    3fde, 3g5t, 3jql, 3nju, 3nna, 3oxp, etc.
    """
    pdb_file_path = "pdbb/tests/pdb/files/1g8t.pdb"
    pdb_id = "1g8t"
    structure = get_structure(pdb_file_path, pdb_id)
    result = check_beq(structure)
    eq_(result["beq_identical"], 0.999124343257443)
    eq_(result["correct_uij"], True)


@raises(ValueError)
def test_check_beq_structure_none():
    """Tests that check_beq raises a value error if structure is None."""
    check_beq(None)


@raises(AssertionError)
def test_check_combinations_error():
    """Tests check_combinations of U values."""
    anisou = [1, 2, 3, 4, 5]
    check_combinations(anisou, 0, 0)


def test_check_combinations_first():
    """Tests check_combinations with standard combination of U values."""
    uij = [4, 5, 6, 1, 0, 0]
    beq = 8*np.pi**2 * (uij[0] + uij[1] + uij[2])/3
    tol = 1e-08
    reproduced = check_combinations(uij, beq, tol, check_first=True)
    eq_(reproduced, True)


def test_check_combinations():
    """Tests check_combinations with non-standard combination of U values."""
    uij = [1, 5, 6, 4, 0, 0]
    beq = 8*np.pi**2 * (uij[3] + uij[1] + uij[2])/3
    tol = 1e-08
    reproduced = check_combinations(uij, beq, tol, check_first=True)
    eq_(reproduced, True)

    reproduced = check_combinations(uij, beq, tol)
    eq_(reproduced, True)


def test_determine_b_group_protein_overall():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/1etu.pdb"
    pdb_id = "1etu"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_determine_b_group_protein_overall_2():
    """Tests that b_group is correctly determined.

        1az2 is has 16.67 for all N CA and C atoms
        and 20.54 for O and side-chain atoms.
        overall is probably the best term.
    """
    pdb_file_path = "pdbb/tests/pdb/files/1az2.pdb"
    pdb_id = "1az2"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_determine_b_group_protein_overall_3():
    """Tests that b_group is correctly determined.

        1c2y has the same B-factor everywhere,
        except for some THR and LEU atoms end-of-side-chain
        overall is probably the best term.
    """
    pdb_file_path = "pdbb/tests/pdb/files/1c2y.pdb"
    pdb_id = "1c2y"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_determine_b_group_protein_overall_4():
    """Tests that b_group is correctly determined.

        1av1 has the same B-factor for large stretches
        of residues; overall is probably the best term.
    """
    pdb_file_path = "pdbb/tests/pdb/files/1av1.pdb"
    pdb_id = "1av1"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], "overall")


def test_determine_b_group_protein_1ADP():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/3b8a.pdb"
    pdb_id = "3b8a"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], "residue_1ADP")


def test_determine_b_group_protein_2ADP():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/3zzt.pdb"
    pdb_id = "3zzt"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], "residue_2ADP")


def test_determine_b_group_nucleic_None():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/1av1.pdb"
    pdb_id = "1av1"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["nucleic_b"], None)


def test_determine_b_group_nucleic_2ADP():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/1hlz.pdb"
    pdb_id = "1hlz"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], "residue_1ADP")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], "residue_2ADP")


def test_determine_b_group_protein_individual():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/1crn.pdb"
    pdb_id = "1crn"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], "individual")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_determine_b_group_protein_individual2():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/2er0.pdb"
    pdb_id = "2er0"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], "individual")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_determine_b_group_nucleic_individual():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/100d.pdb"
    pdb_id = "100d"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], None)
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], "individual")


def test_determine_b_group_protein_no_b_factors():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/1mcb.pdb"
    pdb_id = "1mcb"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], "no_b-factors")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_determine_b_group_calpha_false():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/1av1.pdb"
    pdb_id = "1av1"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["calpha_only"], False)


def test_determine_b_group_calpha_true():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/1a1q.pdb"
    pdb_id = "1a1q"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["calpha_only"], True)


def test_determine_b_group_calpha_true_2():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/1efg.pdb"
    pdb_id = "1efg"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], True)
    eq_(result["nucleic_b"], None)
    eq_(result["phos_only"], False)


def test_determine_b_group_calpha_true_3():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/1efg.pdb"
    pdb_id = "1efg"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], True)


def test_determine_b_group_calpha_true_phos_true():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/3cw1.pdb"
    pdb_id = "3cw1"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], True)
    eq_(result["nucleic_b"], "overall")
    eq_(result["phos_only"], True)


def test_determine_b_group_too_short():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/438d.pdb"
    pdb_id = "438d"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], None)
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)
    eq_(result["phos_only"], False)


def test_determine_b_group_nousefulres():
    """Tests that b_group is correctly determined."""
    pdb_file_path = "pdbb/tests/pdb/files/1c0q.pdb"
    pdb_id = "1c0q"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure)
    eq_(result["protein_b"], "individual")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)
    eq_(result["phos_only"], False)


@raises(IOError)
def test_get_structure_invalid_path():
    """Tests get_structure."""
    pdb_file_path = "1crn.pdb"
    pdb_id = "1crn"
    get_structure(pdb_file_path, pdb_id, verbose=True)


def test_get_structure_pdbid_none():
    """Tests get_structure."""
    pdb_file_path = "pdbb/tests/pdb/files/1crn.pdb"
    pdb_id = None
    ok_(get_structure(pdb_file_path, pdb_id, verbose=True))


def test_is_calpha_trace_true():
    """Tests is_calpha_trace."""
    pdb_file_path = "pdbb/tests/pdb/files/1efg.pdb"
    pdb_id = "1efg"
    structure = get_structure(pdb_file_path, pdb_id)
    chain = structure.get_chains().next()
    result = is_calpha_trace(chain)
    eq_(result, True)


def test_is_calpha_trace_true_unk():
    """Tests is_calpha_trace."""
    pdb_file_path = "pdbb/tests/pdb/files/1efg.pdb"
    pdb_id = "1efg"
    structure = get_structure(pdb_file_path, pdb_id)
    chains = list(structure.get_chains())
    result = is_calpha_trace(chains[1])
    eq_(result, True)
    result = is_calpha_trace(chains[2])
    eq_(result, True)


def test_is_calpha_trace_true_2():
    """Tests is_calpha_trace."""
    pdb_file_path = "pdbb/tests/pdb/files/3cw1.pdb"
    pdb_id = "3cw1"
    structure = get_structure(pdb_file_path, pdb_id)
    chains = structure.get_chains()
    for c in chains:
        cid = c.get_id()
        if cid in ("D", "S", "T", "U", "A", "H", "I", "J", "B", "M", "N", "O",
                   "C", "P", "Q", "R", "1", "2", "F", "Z", "E", "W", "X", "Y",
                   "3", "5", "G", "6", "7", "8", "K", "0", "9", "L", "I"):
            result = is_calpha_trace(c)
            eq_(result, True)


def test_is_calpha_trace_false_prot():
    """Tests is_calpha_trace."""
    pdb_file_path = "pdbb/tests/pdb/files/1crn.pdb"
    pdb_id = "1crn"
    structure = get_structure(pdb_file_path, pdb_id)
    chain = structure.get_chains().next()
    result = is_calpha_trace(chain)
    eq_(result, False)


def test_is_calpha_trace_false_nuc():
    """Tests is_calpha_trace."""
    pdb_file_path = "pdbb/tests/pdb/files/100d.pdb"
    pdb_id = "100d"
    structure = get_structure(pdb_file_path, pdb_id)
    chain = structure.get_chains().next()
    result = is_calpha_trace(chain)
    eq_(result, False)


@raises(AttributeError)
def test_is_calpha_trace_none():
    """Tests is_calpha_trace."""
    is_calpha_trace(None)


def test_is_phos_trace_true():
    """Tests is_phos_trace."""
    pdb_file_path = "pdbb/tests/pdb/files/3cw1.pdb"
    pdb_id = "3cw1"
    structure = get_structure(pdb_file_path, pdb_id)
    chains = structure.get_chains()
    for c in chains:
        cid = c.get_id()
        if cid in ("V", "v", "w", "x"):
            result = is_phos_trace(c)
            eq_(result, True)


def test_is_phos_trace_false_nuc_chain():
    """Tests is_phos_trace."""
    pdb_file_path = "pdbb/tests/pdb/files/100d.pdb"
    pdb_id = "100d"
    structure = get_structure(pdb_file_path, pdb_id)
    chain = structure.get_chains().next()
    result = is_phos_trace(chain)
    eq_(result, False)


def test_is_phos_trace_false_prot_chain():
    """Tests is_phos_trace."""
    pdb_file_path = "pdbb/tests/pdb/files/1crn.pdb"
    pdb_id = "1crn"
    structure = get_structure(pdb_file_path, pdb_id)
    chain = structure.get_chains().next()
    result = is_phos_trace(chain)
    eq_(result, False)


@raises(AttributeError)
def test_is_phos_trace_none():
    """Tests is_phos_trace."""
    is_phos_trace(None)


def test_has_amino_acid_backbone_true():
    """Tests has_amino_acid_backbone with protein."""
    pdb_file_path = "pdbb/tests/pdb/files/1crn.pdb"
    pdb_id = "1crn"
    structure = get_structure(pdb_file_path, pdb_id)
    result = has_amino_acid_backbone(structure.get_residues().next())
    eq_(result, True)


def test_has_amino_acid_backbone_false():
    """Tests has_amino_acid_backbone with dna/rna."""
    pdb_file_path = "pdbb/tests/pdb/files/100d.pdb"
    pdb_id = "100d"
    structure = get_structure(pdb_file_path, pdb_id)
    residues = structure.get_residues()
    residue = residues.next()
    result = has_amino_acid_backbone(residue)
    eq_(result, False)

    residue = residues.next()
    result = has_amino_acid_backbone(residue)
    eq_(result, False)


def test_has_sugar_phosphate_backbone_false():
    """Tests has_amino_acid_backbone with protein."""
    pdb_file_path = "pdbb/tests/pdb/files/1crn.pdb"
    pdb_id = "1crn"
    structure = get_structure(pdb_file_path, pdb_id)
    result = has_sugar_phosphate_backbone(structure.get_residues().next())
    eq_(result, False)


def test_has_sugar_phosphate_backbone_true():
    """Tests has_amino_acid_backbone with dna/rna."""
    pdb_file_path = "pdbb/tests/pdb/files/100d.pdb"
    pdb_id = "100d"
    structure = get_structure(pdb_file_path, pdb_id)
    residues = structure.get_residues()
    residue = residues.next()
    result = has_sugar_phosphate_backbone(residue)
    eq_(result, False)  # DNA, misses P

    residue = residues.next()
    result = has_sugar_phosphate_backbone(residue)
    eq_(result, True)

    residues = list(structure.get_chains().next().get_residues())
    residue = residues[9]  # Residue 10 has the phosphate..
    result = has_sugar_phosphate_backbone(residue)
    eq_(result, True)


def test_is_heavy_backbone():
    """Tests is_heavy_backbone with protein, dna and rna."""
    pdb_file_path = "pdbb/tests/pdb/files/1crn.pdb"
    pdb_id = "1crn"
    structure = get_structure(pdb_file_path, pdb_id)
    result = is_heavy_backbone(structure.get_atoms().next())
    eq_(result, True)

    pdb_file_path = "pdbb/tests/pdb/files/100d.pdb"
    pdb_id = "100d"
    structure = get_structure(pdb_file_path, pdb_id)
    residues = structure.get_residues()
    residue = residues.next()
    atom = residue.get_list()[0]
    result = is_heavy_backbone(atom)
    eq_(result, True)

    residue = residues.next()
    atom = residue.get_list()[0]
    result = is_heavy_backbone(atom)
    eq_(result, True)


def test_is_nucleic_chain_false():
    """Tests is_nucleic_chain with protein."""
    pdb_file_path = "pdbb/tests/pdb/files/1crn.pdb"
    pdb_id = "1crn"
    structure = get_structure(pdb_file_path, pdb_id)
    result = is_nucleic_chain(structure.get_chains().next())
    eq_(result, False)


def test_is_nucleic_chain_true():
    """Tests is_nucleic_chain with dna/rna."""
    pdb_file_path = "pdbb/tests/pdb/files/100d.pdb"
    pdb_id = "100d"
    structure = get_structure(pdb_file_path, pdb_id)
    result = is_nucleic_chain(structure.get_chains().next())
    eq_(result, True)


def test_is_protein_chain_true():
    """Tests is_protein_chain with protein."""
    pdb_file_path = "pdbb/tests/pdb/files/1crn.pdb"
    pdb_id = "1crn"
    structure = get_structure(pdb_file_path, pdb_id)
    result = is_protein_chain(structure.get_chains().next())
    eq_(result, True)


def test_is_protein_chain_false():
    """Tests is_protein_chain with dna/rna."""
    pdb_file_path = "pdbb/tests/pdb/files/100d.pdb"
    pdb_id = "100d"
    structure = get_structure(pdb_file_path, pdb_id)
    result = is_protein_chain(structure.get_chains().next())
    eq_(result, False)


def test_multiply_bfactors_8pipi():
    """Tests that bfactors are correctly multiplied by 8*pi^2."""
    pdb_file_path = "pdbb/tests/pdb/files/1crn.pdb"
    pdb_id = "1crn"
    structure = get_structure(pdb_file_path, pdb_id)
    expected = []
    for atom in structure.get_atoms():
        expected.append(8*np.pi**2 * atom.get_bfactor())

    s_mult = multiply_bfactors_8pipi(structure)
    bvalues = []
    for atom in s_mult.get_atoms():
        bvalues.append(atom.get_bfactor())
    eq_(bvalues, expected)
