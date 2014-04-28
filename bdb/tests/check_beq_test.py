from nose.tools import eq_

from bdb.check_beq import determine_b_group


def test_check_determine_b_group_protein_overall():
    pdb_file_path = "bdb/tests/pdb/files/1etu.pdb"
    pdb_id = "1etu"
    result = determine_b_group(pdb_file_path, pdb_id)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_check_determine_b_group_protein_overall_2():
    """ 1az2 is has 16.67 for all N CA and C atoms
        and 20.54 for O and side-chain atoms.
        overall is probably the best term.
    """
    pdb_file_path = "bdb/tests/pdb/files/1az2.pdb"
    pdb_id = "1az2"
    result = determine_b_group(pdb_file_path, pdb_id)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_check_determine_b_group_protein_overall_3():
    """ 1c2y has the same B-factor everywhere,
        except for some THR and LEU atoms end-of-side-chain
        overall is probably the best term.
    """
    pdb_file_path = "bdb/tests/pdb/files/1c2y.pdb"
    pdb_id = "1c2y"
    result = determine_b_group(pdb_file_path, pdb_id)
    eq_(result["protein_b"], "overall")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_check_determine_b_group_protein_1ADP():
    pdb_file_path = "bdb/tests/pdb/files/1av1.pdb"
    pdb_id = "1av1"
    result = determine_b_group(pdb_file_path, pdb_id)
    eq_(result["protein_b"], "residue_1ADP")


def test_check_determine_b_group_nucleic_None():
    pdb_file_path = "bdb/tests/pdb/files/1av1.pdb"
    pdb_id = "1av1"
    result = determine_b_group(pdb_file_path, pdb_id)
    eq_(result["nucleic_b"], None)


def test_check_determine_b_group_nucleic_2ADP():
    pdb_file_path = "bdb/tests/pdb/files/1hlz.pdb"
    pdb_id = "1hlz"
    result = determine_b_group(pdb_file_path, pdb_id)
    eq_(result["protein_b"], "residue_1ADP")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], "residue_2ADP")


def test_check_determine_b_group_protein_individual():
    pdb_file_path = "bdb/tests/pdb/files/1crn.pdb"
    pdb_id = "1crn"
    result = determine_b_group(pdb_file_path, pdb_id)
    eq_(result["protein_b"], "individual")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_check_determine_b_group_nucleic_individual():
    pdb_file_path = "bdb/tests/pdb/files/100d.pdb"
    pdb_id = "100d"
    result = determine_b_group(pdb_file_path, pdb_id)
    eq_(result["protein_b"], None)
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], "individual")


def test_check_determine_b_group_protein_no_b_factors():
    pdb_file_path = "bdb/tests/pdb/files/1mcb.pdb"
    pdb_id = "1mcb"
    result = determine_b_group(pdb_file_path, pdb_id)
    eq_(result["protein_b"], "no_b-factors")
    eq_(result["calpha_only"], False)
    eq_(result["nucleic_b"], None)


def test_check_determine_b_group_calpha_false():
    pdb_file_path = "bdb/tests/pdb/files/1av1.pdb"
    pdb_id = "1av1"
    result = determine_b_group(pdb_file_path, pdb_id)
    eq_(result["calpha_only"], False)


def test_check_determine_b_group_calpha_true():
    pdb_file_path = "bdb/tests/pdb/files/1a1q.pdb"
    pdb_id = "1a1q"
    result = determine_b_group(pdb_file_path, pdb_id)
    eq_(result["calpha_only"], True)


