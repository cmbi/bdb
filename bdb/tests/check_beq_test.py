from nose.tools import eq_, ok_, raises

from bdb.check_beq import (determine_b_group, get_structure, is_calpha_trace,
                           is_phos_trace)

from Bio.PDB.PDBExceptions import PDBConstructionException



def test_determine_b_group_protein_overall():
    pdb_file_path = "bdb/tests/pdb/files/1etu.pdb"
    pdb_id = "1etu"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure, pdb_id)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_determine_b_group_protein_overall_2():
    """ 1az2 is has 16.67 for all N CA and C atoms
        and 20.54 for O and side-chain atoms.
        overall is probably the best term.
    """
    pdb_file_path = "bdb/tests/pdb/files/1az2.pdb"
    pdb_id = "1az2"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure, pdb_id)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_determine_b_group_protein_overall_3():
    """ 1c2y has the same B-factor everywhere,
        except for some THR and LEU atoms end-of-side-chain
        overall is probably the best term.
    """
    pdb_file_path = "bdb/tests/pdb/files/1c2y.pdb"
    pdb_id = "1c2y"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure, pdb_id)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_determine_b_group_protein_1ADP():
    pdb_file_path = "bdb/tests/pdb/files/1av1.pdb"
    pdb_id = "1av1"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure, pdb_id)
    eq_(result["protein_b"], "residue_1ADP")


def test_determine_b_group_nucleic_None():
    pdb_file_path = "bdb/tests/pdb/files/1av1.pdb"
    pdb_id = "1av1"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure, pdb_id)
    eq_(result["nucleic_b"], None)


def test_determine_b_group_nucleic_2ADP():
    pdb_file_path = "bdb/tests/pdb/files/1hlz.pdb"
    pdb_id = "1hlz"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure, pdb_id)
    eq_(result["protein_b"], "residue_1ADP")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], "residue_2ADP")


def test_determine_b_group_protein_individual():
    pdb_file_path = "bdb/tests/pdb/files/1crn.pdb"
    pdb_id = "1crn"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure, pdb_id)
    eq_(result["protein_b"], "individual")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_determine_b_group_nucleic_individual():
    pdb_file_path = "bdb/tests/pdb/files/100d.pdb"
    pdb_id = "100d"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure, pdb_id)
    eq_(result["protein_b"], None)
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], "individual")


def test_determine_b_group_protein_no_b_factors():
    pdb_file_path = "bdb/tests/pdb/files/1mcb.pdb"
    pdb_id = "1mcb"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure, pdb_id)
    eq_(result["protein_b"], "no_b-factors")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_determine_b_group_calpha_false():
    pdb_file_path = "bdb/tests/pdb/files/1av1.pdb"
    pdb_id = "1av1"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure, pdb_id)
    eq_(result["calpha_only"], False)


def test_determine_b_group_calpha_true():
    pdb_file_path = "bdb/tests/pdb/files/1a1q.pdb"
    pdb_id = "1a1q"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure, pdb_id)
    eq_(result["calpha_only"], True)


def test_determine_b_group_calpha_true_2():
    pdb_file_path = "bdb/tests/pdb/files/1efg.pdb"
    pdb_id = "1efg"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure, pdb_id)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], True)
    eq_(result["nucleic_b"], None)
    eq_(result["phos_only"], False)

def test_determine_b_group_calpha_true_3():
    pdb_file_path = "bdb/tests/pdb/files/1efg.pdb"
    pdb_id = "1efg"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure, pdb_id)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], True)


def test_determine_b_group_calpha_true_phos_true():
    pdb_file_path = "bdb/tests/pdb/files/3cw1.pdb"
    pdb_id = "3cw1"
    structure = get_structure(pdb_file_path, pdb_id)
    result = determine_b_group(structure, pdb_id)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], True)
    eq_(result["nucleic_b"], "overall")
    eq_(result["phos_only"], True)


@raises(IOError)
def test_get_structure_invalid_path():
    pdb_file_path = "1crn.pdb"
    pdb_id = "1crn"
    structure = get_structure(pdb_file_path, pdb_id, verbose=True)


def test_get_structure_pdbid_none():
    pdb_file_path = "bdb/tests/pdb/files/1crn.pdb"
    pdb_id = None
    structure = get_structure(pdb_file_path, pdb_id, verbose=True)


def test_get_structure_parse_error():
    pdb_file_path = "bdb/tests/pdb/files/4aph.pdb"
    pdb_id = "4aph"
    ok_(get_structure(pdb_file_path, pdb_id, verbose=False))


def test_is_calpha_trace_true():
    pdb_file_path = "bdb/tests/pdb/files/1efg.pdb"
    pdb_id = "1efg"
    structure = get_structure(pdb_file_path, pdb_id)
    chain = structure.get_chains().next()
    result = is_calpha_trace(chain)
    eq_(result, True)


def test_is_calpha_trace_true_unk():
    pdb_file_path = "bdb/tests/pdb/files/1efg.pdb"
    pdb_id = "1efg"
    structure = get_structure(pdb_file_path, pdb_id)
    chains = list(structure.get_chains())
    result = is_calpha_trace(chains[1])
    eq_(result, True)
    result = is_calpha_trace(chains[2])
    eq_(result, True)


def test_is_calpha_trace_true_2():
    pdb_file_path = "bdb/tests/pdb/files/3cw1.pdb"
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
    pdb_file_path = "bdb/tests/pdb/files/1crn.pdb"
    pdb_id = "1crn"
    structure = get_structure(pdb_file_path, pdb_id)
    chain = structure.get_chains().next()
    result = is_calpha_trace(chain)
    eq_(result, False)


def test_is_calpha_trace_false_nuc():
    pdb_file_path = "bdb/tests/pdb/files/100d.pdb"
    pdb_id = "100d"
    structure = get_structure(pdb_file_path, pdb_id)
    chain = structure.get_chains().next()
    result = is_calpha_trace(chain)
    eq_(result, False)


@raises(AttributeError)
def test_is_calpha_trace_none():
    result = is_calpha_trace(None)


def test_is_phos_trace_true():
    pdb_file_path = "bdb/tests/pdb/files/3cw1.pdb"
    pdb_id = "3cw1"
    structure = get_structure(pdb_file_path, pdb_id)
    chains = structure.get_chains()
    for c in chains:
        cid = c.get_id()
        if cid in ("V", "v", "w", "x"):
            result = is_phos_trace(c)
            eq_(result, True)


def test_is_phos_trace_false_nuc_chain():
    pdb_file_path = "bdb/tests/pdb/files/100d.pdb"
    pdb_id = "100d"
    structure = get_structure(pdb_file_path, pdb_id)
    chain = structure.get_chains().next()
    result = is_phos_trace(chain)
    eq_(result, False)


def test_is_phos_trace_false_prot_chain():
    pdb_file_path = "bdb/tests/pdb/files/1crn.pdb"
    pdb_id = "1crn"
    structure = get_structure(pdb_file_path, pdb_id)
    chain = structure.get_chains().next()
    result = is_phos_trace(chain)
    eq_(result, False)


@raises(AttributeError)
def test_is_phos_trace_none():
    result = is_phos_trace(None)


