from nose.tools import eq_, raises

from bdb.pdb.parser import parse_pdb_file, parse_exp_methods


@raises(ValueError)
def test_parser_invalid_file():
    p = parse_pdb_file("1crn.pdb")


def test_parser():
    pdb_records = parse_pdb_file("bdb/tests/pdb/1crn.pdb")
    eq_(len(pdb_records), 27)
    eq_(len(pdb_records["HEADER"]), 1)
    eq_(len(pdb_records["EXPDTA"]), 1)
    eq_(len(pdb_records["ATOM  "]), 327)
    eq_(len(pdb_records["REMARK"]), 224)


def test_exp_method_single_line_multiple():
    records = {"EXPDTA": ["    NEUTRON DIFFRACTION; X-RAY DIFFRACTION",]}
    exp_methods = parse_exp_methods(records)
    eq_(len(exp_methods), 2)
    eq_(exp_methods[0], "NEUTRON DIFFRACTION")
    eq_(exp_methods[1], "X-RAY DIFFRACTION")


def test_exp_method_single_line_single():
    records = {"EXPDTA": ["    X-RAY DIFFRACTION",]}
    exp_methods = parse_exp_methods(records)
    eq_(len(exp_methods), 1)
    eq_(exp_methods[0], "X-RAY DIFFRACTION")


@raises(ValueError)
def test_exp_method_none():
    records = {"REMARK": ["TEST",]}
    exp_methods = parse_exp_methods(records)
