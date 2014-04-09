from nose.tools import eq_, raises

from bdb.pdb.parser import (parse_pdb_file, parse_exp_methods, parse_btype,
                            parse_other_ref_remarks, is_bmsqav,
                            parse_format_date_version, parse_num_tls_groups,
                            parse_ref_prog)


@raises(ValueError)
def test_parser_invalid_file():
    parse_pdb_file("1crn.pdb")


def test_parser():
    pdb_records = parse_pdb_file("bdb/tests/pdb/1crn.pdb")
    eq_(len(pdb_records), 27)
    eq_(len(pdb_records["HEADER"]), 1)
    eq_(len(pdb_records["EXPDTA"]), 1)
    eq_(len(pdb_records["ATOM  "]), 327)
    eq_(len(pdb_records["REMARK"]), 224)


def test_parse_exp_method_single_line_multiple():
    records = {"EXPDTA": ["    NEUTRON DIFFRACTION; X-RAY DIFFRACTION", ]}
    exp_methods = parse_exp_methods(records)
    eq_(len(exp_methods), 2)
    eq_(exp_methods[0], "NEUTRON DIFFRACTION")
    eq_(exp_methods[1], "X-RAY DIFFRACTION")


def test_parse_exp_method_single_line_single():
    records = {"EXPDTA": ["    X-RAY DIFFRACTION", ]}
    exp_methods = parse_exp_methods(records)
    eq_(len(exp_methods), 1)
    eq_(exp_methods[0], "X-RAY DIFFRACTION")


@raises(ValueError)
def test_parse_exp_method_none():
    records = {"REMARK": ["TEST", ]}
    parse_exp_methods(records)


def test_parse_btype():
    """
    Tests that the btype value is correctly parsed from the pdb file remark
    records.

    This example is taken from 1ean.
    """
    records = {"REMARK": ["  3   B VALUE TYPE : LIKELY RESIDUAL", ]}
    btype = parse_btype(records)
    eq_(btype, "LIKELY RESIDUAL")


def test_parse_other_ref_remarks():
    records = {"REMARK": [
        "  3  OTHER REFINEMENT REMARKS: ANISOTROPIC B-FACTOR REFINEMENT",
        "  3  PERFORMED WITH SHELXL-97 GIVES AN R-FACTOR OF 13.4% AND AN R-",
        "  3  FREE OF 0.154.", ]}
    other_ref_remarks = parse_other_ref_remarks(records)
    expected = reduce(
        lambda x, y: x.replace("  3  OTHER REFINEMENT REMARKS: ", "")
        + '\n' + y[5:], records["REMARK"])
    eq_(other_ref_remarks, expected)


def test_is_bmsqav_true():
    ref_remarks = "THE QUANTITY GIVEN IN THE TEMPERATURE\n" \
        "FACTOR FIELD OF THE *ATOM* AND *HETATM* RECORDS BELOW IS U**2," \
        "WHICH IS THE MEAN-SQUARE AMPLITUDE OF ATOMIC VIBRATION. THE" \
        "TEMPERATURE FACTOR, B, CAN BE DERIVED BY THE FOLLOWING RELATION -" \
        "B = 8 * (PI)**2 * U**2. IT IS AN INDICATION OF POSSIBLE ERRORS IN" \
        "THE REFINEMENT THAT SOME ARE SLIGHTLY NEGATIVE."
    eq_(is_bmsqav(ref_remarks), True)


def test_is_bmsqav_false():
    ref_remarks = ""
    eq_(is_bmsqav(ref_remarks), False)


def test_parse_format_date_version():
    records = {"REMARK":
               ["  4 12GS COMPLIES WITH FORMAT V. 3.30, 13-JUL-11", ]}
    (f_vers, f_date) = parse_format_date_version(records)
    eq_(f_vers, 3.3)
    eq_(f_date, "13-JUL-11")


def test_parse_format_date_version_none():
    records = {"REMARK": ["", ]}
    (f_vers, f_date) = parse_format_date_version(records)
    eq_(f_vers, None)
    eq_(f_date, None)


def test_parse_num_tls_groups():
    records = {"REMARK": ["  3   NUMBER OF TLS GROUPS  : 6", ]}
    num_tls_groups = parse_num_tls_groups(records)
    eq_(num_tls_groups, 6)


def test_parse_num_tls_groups_null():
    records = {"REMARK": ["  3   NUMBER OF TLS GROUPS  : NULL", ]}
    num_tls_groups = parse_num_tls_groups(records)
    eq_(num_tls_groups, None)


def test_parse_num_tls_groups_non_existant():
    records = {"REMARK": ["", ]}
    num_tls_groups = parse_num_tls_groups(records)
    eq_(num_tls_groups, None)


def test_parse_ref_prog():
    records = {"REMARK": ["  3   PROGRAM     : X-PLOR 3.1", ]}
    ref_prog = parse_ref_prog(records)
    eq_(ref_prog, "X-PLOR 3.1")


def test_parse_ref_prog_null():
    records = {"REMARK": ["  3   PROGRAM     : NULL", ]}
    ref_prog = parse_ref_prog(records)
    eq_(ref_prog, None)


def test_parse_ref_prog_no_refinement():
    records = {"REMARK": ["  3   PROGRAM     : NO REFINEMENT", ]}
    ref_prog = parse_ref_prog(records)
    eq_(ref_prog, None)


def test_parse_ref_prog_none():
    records = {"REMARK": ["  3   PROGRAM     : NONE", ]}
    ref_prog = parse_ref_prog(records)
    eq_(ref_prog, None)
