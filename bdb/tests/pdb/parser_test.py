from nose.tools import eq_, raises

from bdb.pdb.parser import (parse_pdb_file, parse_exp_methods, parse_btype,
                            parse_other_ref_remarks, is_bmsqav,
                            parse_format_date_version, parse_num_tls_groups,
                            parse_ref_prog)


@raises(ValueError)
def test_parser_invalid_file():
    parse_pdb_file("1crn.pdb")


def test_parser():
    pdb_records = parse_pdb_file("bdb/tests/pdb/files/1crn.pdb")
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


def test_parse_btype_residual():
    """
    Tests that the btype value is correctly parsed from the pdb file remark
    records.

    This example is taken from e.g.
    1zup, 1zuz, 1zv4, 1zv7, 1zv8, 1zva, 1zvb, 1zvn, 1zvp, 1zvq, 1zvt, 1zvu,
    2zvk, 2zvl, 2zvm, 1zw5, 1zwz, 2zwn, 1zx1, 1zx2, 1zx5, 1zx6, 1zx8, 1zxk,
    1zxt, 2zxe, 2zxj, 2zxt, 1zyb, 1zzi, 1zzj, 2zz1, 2zz2, 2zz5, 2zz7, 2zz8,
    etc. etc.
    """
    records = {"REMARK": ["  3   B VALUE TYPE : LIKELY RESIDUAL", ]}
    btype = parse_btype(records)
    eq_(btype, "residual")


def test_parse_btype_residual_trailing_whitespace():
    """
    Tests that the btype value is correctly parsed from the pdb file remark
    records.

    This example is taken from e.g. 1ean.
    """
    records = {"REMARK": ["  3   B VALUE TYPE : LIKELY RESIDUAL           ", ]}
    btype = parse_btype(records)
    eq_(btype, "residual")


def test_parse_btype_unverified():
    """
    Tests that the btype value is correctly parsed from the pdb file remark
    records.

    This example is taken from e.g.
    2a2t, 3b4v, 3b9k, 3be8, 2bnx, 2bwg, 2byc, 2byo, 2byp, 2byr, 2c1j, 2c2g,
    2c2n, 2c43, 2c8t, 2c95, 2cdg, 2cfo, 2cfz, 2ch6, 2ci9, 2ck3, 2ckd, 2ckq,
    2clq, 3eit, 3ekt, 3ekv, 3eky, 3ffl, 3gg8, 3hvg, 2io4, 3io2, 2ixp, 2iyg,
    2izs, 2izt, 2j1n, 2j24, 2j7x, 2j8x, 2jc6, 2jdi, 2je0, 2jf6, 2jfd, 2jfk,
    2jg3, 2jgs, 2jhp, 2jjy, 2oyu, 2p28, 2uv0, 1uwx, 2uz2, 2v0t, 2v24, 2v3v,
    2v79, 2v7l, 2v8d, 2v8z, 2v9s, 2va0, 2vag, 2vc2, 2vdk, 2ve7, 2vef, 2vfc,
    2vg3, 2vgf, 2vgi, 2vka, 2vlf, 2vnu, 2vob, 2vsz, 2vu0, 2vu1, 2vv8, 2vvc,
    2vx2, 2vxk, 2vze, 2w0d, 2w3d, 2w6t, 2w75, 2w77, 1w89, 2w8f, 2wnc, 2wvo,
    2x23, 2x3w
    """
    records = {"REMARK": ["  3   B VALUE TYPE : UNVERIFIED", ]}
    btype = parse_btype(records)
    eq_(btype, "unverified")


def test_parse_btype_unverified_trailing_whitespace():
    """
    Tests that the btype value is correctly parsed from the pdb file remark
    records.

    This example is taken from e.g. 2x3w.
    """
    records = {"REMARK": ["  3   B VALUE TYPE : UNVERIFIED                ", ]}
    btype = parse_btype(records)
    eq_(btype, "unverified")


def test_parse_other_ref_remarks():
    """
    Tests that other refinement remarks are concatenated correctly from the
    pdb file remark records.

    This example is taken from 1a0m.
    """
    records = {"REMARK": [
        "  3                                                                      ",
        "  3  OTHER REFINEMENT REMARKS: ANISOTROPIC B-FACTOR REFINEMENT           ",
        "  3  PERFORMED WITH SHELXL-97 GIVES AN R-FACTOR OF 13.4% AND AN R-       ",
        "  3  FREE OF 0.154.                                                      ",
        "  4                                                                      ",
        "  4 1A0M COMPLIES WITH FORMAT V. 3.15, 01-DEC-08                         ",
        ]}
    other_ref_remarks = parse_other_ref_remarks(records)
    expected = "ANISOTROPIC B-FACTOR REFINEMENT PERFORMED WITH SHELXL-97 " \
               "GIVES AN R-FACTOR OF 13.4% AND AN R- FREE OF 0.154."
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


def test_parse_ref_prog_null_trailing_whitespace():
    records = {"REMARK": ["  3   PROGRAM     : NULL                       ", ]}
    ref_prog = parse_ref_prog(records)
    eq_(ref_prog, None)


def test_parse_ref_prog_no_refinement():
    records = {"REMARK": ["  3   PROGRAM     : NO REFINEMENT", ]}
    ref_prog = parse_ref_prog(records)
    eq_(ref_prog, None)


def test_parse_ref_prog_no_refinement_trailing_whitespace():
    records = {"REMARK": ["  3   PROGRAM     : NO REFINEMENT              ", ]}
    ref_prog = parse_ref_prog(records)
    eq_(ref_prog, None)


def test_parse_ref_prog_none():
    records = {"REMARK": ["  3   PROGRAM     : NONE", ]}
    ref_prog = parse_ref_prog(records)
    eq_(ref_prog, None)


def test_parse_ref_prog_none_trailing_whitespace():
    records = {"REMARK": ["  3   PROGRAM     : NONE                       ", ]}
    ref_prog = parse_ref_prog(records)
    eq_(ref_prog, None)
