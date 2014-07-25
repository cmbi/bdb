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
from datetime import datetime
from nose.tools import eq_, raises

from pdbb.pdb.parser import (parse_pdb_file, parse_dep_date, parse_exp_methods,
                             parse_btype, parse_other_ref_remarks, is_bmsqav,
                             parse_format_date_version, parse_num_tls_groups,
                             parse_ref_prog, is_tls_residual, is_tls_sum,
                             get_pdb_header_and_trailer)


@raises(ValueError)
def test_parser_invalid_file():
    parse_pdb_file("1crn.pdb")


def test_parser():
    pdb_records = parse_pdb_file("pdbb/tests/pdb/files/1crn.pdb")
    eq_(len(pdb_records), 27)
    eq_(len(pdb_records["HEADER"]), 1)
    eq_(len(pdb_records["EXPDTA"]), 1)
    eq_(len(pdb_records["ATOM  "]), 327)
    eq_(len(pdb_records["REMARK"]), 224)


def test_parse_dep_date():
    records = parse_pdb_file("pdbb/tests/pdb/files/1crn.pdb")
    dep_date = parse_dep_date(records)
    eq_(dep_date, datetime(1981, 4, 30))

    records = {"HEADER": ["   PLANT PROTEIN                         "
                          "  30-APR-81   1CRN"]}
    parse_dep_date(records)
    eq_(dep_date, datetime(1981, 4, 30))


@raises(ValueError)
def test_parse_dep_date_none():
    records = {}
    parse_dep_date(records)


@raises(ValueError)
def test_parse_dep_date_multiple():
    records = {"HEADER": [" HEADER1", "HEADER2"]}
    parse_dep_date(records)


@raises(ValueError)
def test_parse_dep_date_format():
    records = {"HEADER": ["   PLANT PROTEIN                         "
                          "  30-AAA-81   1CRN"]}
    parse_dep_date(records)
    records = {"HEADER": ["   PLANT PROTEIN                         "
                          "  30-04-81    1CRN"]}
    parse_dep_date(records)
    records = {"HEADER": ["   PLANT PROTEIN                         "
                          " 30-APR-81    1CRN"]}
    parse_dep_date(records)
    records = {"HEADER": ["   PLANT PROTEIN                         "
                          "   30-AAA-81  1CRN"]}
    parse_dep_date(records)
    records = {"HEADER": ["   PLANT PROTEIN                         "
                          "  30-04-1981  1CRN"]}
    parse_dep_date(records)
    records = {"HEADER": ["   PLANT PROTEIN                         "
                          "  30- 04-81  1CRN"]}
    parse_dep_date(records)
    records = {"HEADER": ["   PLANT PROTEIN                         "
                          "  81-APR-30  1CRN"]}
    parse_dep_date(records)
    records = {"HEADER": ["   PLANT PROTEIN                         "
                          "  1981-04-30 1CRN"]}
    parse_dep_date(records)


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


@raises(ValueError)
def test_parse_btype_unexpected():
    """
    Tests that the btype value is correctly parsed from the pdb file remark
    records.
    """
    records = {"REMARK": ["  3   B VALUE TYPE : SOMETHINGELSE", ]}
    parse_btype(records)


def test_parse_btype_unverified_trailing_whitespace():
    """
    Tests that the btype value is correctly parsed from the pdb file remark
    records.

    This example is taken from e.g. 2x3w.
    """
    records = {"REMARK": ["  3   B VALUE TYPE : UNVERIFIED                ", ]}
    btype = parse_btype(records)
    eq_(btype, "unverified")


def test_parse_btype_unverified_trailing_whitespace_none():
    """Trivially test that the btype value is None if not present."""
    records = {"REMARK": ["  3               ", ]}
    btype = parse_btype(records)
    eq_(btype, None)


def test_parse_other_ref_remarks():
    """
    Tests that other refinement remarks are concatenated correctly from the
    pdb file remark records.

    This example is taken from 1a0m.
    """
    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: ANISOTROPIC B-FACTOR REFINEMENT        "
        "   ",
        "  3  PERFORMED WITH SHELXL-97 GIVES AN R-FACTOR OF 13.4% AND AN R-    "
        "   ",
        "  3  FREE OF 0.154.                                                   "
        "   ",
        "  4                                                                   "
        "   ",
        "  4 1A0M COMPLIES WITH FORMAT V. 3.15, 01-DEC-08                      "
        "   ",
        ]}
    other_ref_remarks = parse_other_ref_remarks(records)
    expected = "ANISOTROPIC B-FACTOR REFINEMENT PERFORMED WITH SHELXL-97 " \
               "GIVES AN R-FACTOR OF 13.4% AND AN R- FREE OF 0.154."
    eq_(other_ref_remarks, expected)


def test_parse_other_ref_remarks_none():
    """
    Tests that other refinement remarks are None if not present in the
    pdb file remark records.

    This example is taken from 1crn.
    """
    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: NULL                                   "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    other_ref_remarks = parse_other_ref_remarks(records)
    eq_(other_ref_remarks, None)


def test_is_bmsqav_true_1():
    """Tests that b_msqav is correctly parsed from REMARK 3 records.

    e.g. 1czi, 1ent, 1smr, 4ape, etc.
    """
    ref_remarks = "THE QUANTITY GIVEN IN THE TEMPERATURE\n" \
        "FACTOR FIELD OF THE *ATOM* AND *HETATM* RECORDS BELOW IS U**2," \
        "WHICH IS THE MEAN-SQUARE AMPLITUDE OF ATOMIC VIBRATION. THE" \
        "TEMPERATURE FACTOR, B, CAN BE DERIVED BY THE FOLLOWING RELATION -" \
        "B = 8 * (PI)**2 * U**2. IT IS AN INDICATION OF POSSIBLE ERRORS IN" \
        "THE REFINEMENT THAT SOME ARE SLIGHTLY NEGATIVE."
    eq_(is_bmsqav(ref_remarks), True)

    ref_remarks = "MEAN-SQUARE AMPLITUDE OF ATOMIC VIBRATION"
    eq_(is_bmsqav(ref_remarks), True)


def test_is_bmsqav_true_2():
    """Tests that b_msqav is correctly parsed from REMARK 3 records.

    e.g. 1eed and 5pep
    """
    ref_remarks = "THE ATOMIC TEMPERATURE FACTORS IN THIS ENTRY ARE GIVEN AS" \
        "U VALUES NOT B VALUES.  B VALUES MAY BE CALCULATED BY THE" \
        "THE FOLLOWING: BISO (BISO = 8 PI==2== UISO)."
    eq_(is_bmsqav(ref_remarks), True)

    ref_remarks = "ISOTROPIC UISO VALUES ARE PROVIDED IN THE FIELD THAT" \
        "U VALUES NOT B VALUES.  B VALUES MAY BE CALCULATED BY THE" \
        "USUALLY CONTAINS B VALUES."
    eq_(is_bmsqav(ref_remarks), True)

    ref_remarks = "U**2"
    eq_(is_bmsqav(ref_remarks), True)

    ref_remarks = "UISO"
    eq_(is_bmsqav(ref_remarks), True)


def test_is_bmsqav_false():
    ref_remarks = ""
    eq_(is_bmsqav(ref_remarks), False)

    ref_remarks = "AMPLITUDE OF ATOMIC VIBRATION"
    eq_(is_bmsqav(ref_remarks), False)

    ref_remarks = "U*2"
    eq_(is_bmsqav(ref_remarks), False)

    ref_remarks = "U ISO"
    eq_(is_bmsqav(ref_remarks), False)


def test_is_bmsqav_na():
    ref_remarks = None
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


def test_parse_format_date_version_string():
    records = {"REMARK":
               ["  4 1ABC COMPLIES WITH FORMAT V. 3.30.1, 13-JUL-11", ]}
    (f_vers, f_date) = parse_format_date_version(records)
    eq_(f_vers, None)
    eq_(f_date, "13-JUL-11")


def test_parse_num_tls_groups():
    records = {"REMARK": ["  3   NUMBER OF TLS GROUPS  : 6 ", ]}
    num_tls_groups = parse_num_tls_groups(records)
    eq_(num_tls_groups, 6)

    records = {"REMARK": ["  3   NUMBER OF TLS GROUPS  :   10", ]}
    num_tls_groups = parse_num_tls_groups(records)
    eq_(num_tls_groups, 10)


def test_parse_num_tls_groups_null():
    records = {"REMARK": ["  3   NUMBER OF TLS GROUPS  : NULL ", ]}
    num_tls_groups = parse_num_tls_groups(records)
    eq_(num_tls_groups, None)


def test_parse_num_tls_groups_non_existant():
    records = {"REMARK": ["", ]}
    num_tls_groups = parse_num_tls_groups(records)
    eq_(num_tls_groups, None)


def test_parse_no_ref_prog():
    records = {"REMARK": ["  3                  ", ]}
    ref_prog = parse_ref_prog(records)
    eq_(ref_prog, None)


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


def test_is_tls_residual_remark3():
    """Tests that tls_residual is correctly parsed from REMARK 3 records."""
    records = {"REMARK": ["  3   ATOM RECORD CONTAINS RESIDUAL B FACTORS ONLY",
                          ]}
    tls_residual = is_tls_residual(records)
    eq_(tls_residual, True)


def test_is_tls_residual_false():
    """Tests that tls_residual is correctly parsed from REMARK 3 records."""
    records = {"REMARK": ["  3                             ", ]}
    tls_residual = is_tls_residual(records)
    eq_(tls_residual, False)

    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: NULL                                   "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_residual = is_tls_residual(records)
    eq_(tls_residual, False)


def test_is_tls_residual_other1():
    """Tests that tls_residual is correctly parsed from REMARK 3 records."""
    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: RESIDUAL B FACTORS ONLY                "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_residual = is_tls_residual(records)
    eq_(tls_residual, True)

    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: RESIDUAL   U-  VALUE   ONLY            "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_residual = is_tls_residual(records)
    eq_(tls_residual, True)

    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: RESIDUAL B FACTORS ON LY               "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_residual = is_tls_residual(records)
    eq_(tls_residual, False)


def test_is_tls_residual_other2():
    """Tests that tls_residual is correctly parsed from REMARK 3 records.

    e.g. 2qnr
    """
    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: 1. HYDROGENS HAVE BEEN ADDED IN THE    "
        "   ",
        "  3  RIDING POSITIONS. 2. COOT, O, MOLPROBITY PROGRAMS HAVE ALSO BEEN "
        "   ",
        "  3  USED IN THE REFINEMENT. 3. ATOMIC B-FACTORS SHOWN ARE RESIDUALS  "
        "   ",
        "  3  FROM TLS REFINEMENT. 4. CAVEAT: RESIDUES 84-89 IN CHAIN A AND 74-"
        "   ",
        "  3  90 IN CHAIN B COULD NOT BE RELIABLY ASSIGNED IN THE AMINO ACID   "
        "   ",
        "  3  SEQUENCE AND WERE MODELED AS ALANINES.                           "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_residual = is_tls_residual(records)
    eq_(tls_residual, True)

    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS:   ATOMIC B-VALUE   ARE  RESIDUALS FROM "
        "   ",
        "  3  TLS REFINEMENT                                                   "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_residual = is_tls_residual(records)
    eq_(tls_residual, True)

    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: ATOMIC B VALUES ARE RESIDUALS FROM TLS "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_residual = is_tls_residual(records)
    eq_(tls_residual, True)

    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: ATOMIC U-VALUES ARE RESIDUALS          "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_residual = is_tls_residual(records)
    eq_(tls_residual, True)


def test_is_tls_residual_other3():
    """Tests that tls_residual is correctly parsed from REMARK 3 records.

    The second example is taken from 1uwv.
    """
    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: B FACTORS ARE RESIDUAL B-FACTORS       "
        "   ",
        ]}
    tls_residual = is_tls_residual(records)
    eq_(tls_residual, True)

    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: HYDROGENS HAVE BEEN ADDED IN THE       "
        "   ",
        "  3  RIDING POSITIONS. B-FACTORS ARE RESIDUAL B-FACTORS, (WHICH DO    "
        "   ",
        "  3  NOT INCLUDE THE CONTRIBUTION FROM THE TLS PARAMETERS). USE       "
        "   ",
        "  3  TLSANL (DISTRIBUTED WITH CCP4) TO OBTAIN THE FULL B-FACTOR.      "
        "   ",
        ]}
    tls_residual = is_tls_residual(records)
    eq_(tls_residual, True)

    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: HYDROGENS HAVE BEEN ADDED IN THE       "
        "   ",
        "  3  RIDING POSITIONS. B-FACTORS ARE RESIDUAL B-FACTORS, (WHICH DO    "
        "   ",
        "  3  NOT INCLUDE THE CONTRIBUTION FROM THE TLS PARAMETERS).           "
        "   ",
        ]}
    tls_residual = is_tls_residual(records)
    eq_(tls_residual, True)

    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: HYDROGENS HAVE BEEN ADDED IN THE       "
        "   ",
        "  3  RIDING POSITIONS. B-FACTORS ARE RESIDUAL B-FACTORS,              "
        "   ",
        "  3  USE TLSANL (DISTRIBUTED WITH CCP4) TO OBTAIN THE FULL B-FACTOR.  "
        "   ",
        ]}
    tls_residual = is_tls_residual(records)
    eq_(tls_residual, True)


def test_is_tls_sum_remark3():
    """Tests that tls_sum is correctly parsed from REMARK 3 records."""
    records = {"REMARK": ["  3   ATOM RECORD CONTAINS SUM OF TLS AND "
                          "RESIDUAL B FACTORS", ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, True)


def test_is_tls_sum_false():
    """Tests that tls_sum is correctly parsed from REMARK 3 records."""
    records = {"REMARK": ["  3                             ", ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, False)

    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: NULL                                   "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, False)

    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: NO HINTS ABOUT B FACTORS               "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, False)


def test_is_tls_sum_other1():
    """Tests that tls_sum is correctly parsed from REMARK 3 records.

    e.g. 2r99
    """
    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: HYDROGENS HAVE BEEN ADDED IN THE       "
        "   ",
        "  3  RIDING POSITIONS, ATOM RECORD CONTAINS SUM OF TLS AND RESIDUAL   "
        "   ",
        "  3  B FACTORS, ANISOU RECORD CONTAINS SUM OF TLS AND RESIDUAL U      "
        "   ",
        "  3  FACTORS                                                          "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, True)

    records = {"REMARK": [
        "  3  OTHER REFINEMENT REMARKS:                                        "
        "   ",
        "  3  SUM OF   TLS    AND RESIDUAL   U  VALUE                          "
        "   ",
        ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, True)


def test_is_tls_sum_other2():
    """Tests that tls_sum is correctly parsed from REMARK 3 records."""
    records = {"REMARK": [
        "  3  OTHER REFINEMENT REMARKS:                                        "
        "   ",
        "  3  B FACTORS WITH  TLS      ADDED                                   "
        "   ",
        ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, True)


def test_is_tls_sum_other3():
    """Tests that tls_sum is correctly parsed from REMARK 3 records.

    e.g. 2wnj
    """
    records = {"REMARK": [
        "  3                                                                   "
        "   ",
        "  3  OTHER REFINEMENT REMARKS: HYDROGENS HAVE BEEN ADDED IN THE       "
        "   ",
        "  3   RIDING POSITIONS. GLOBAL B-FACTORS, CONTAINING RESIDUAL         "
        "   ",
        "  3   AND TLS COMPONENT HAVE BEEN DEPOSITED                           "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, True)

    records = {"REMARK": [
        "  3  OTHER REFINEMENT REMARKS:                                        "
        "   ",
        "  3  B-FACTORS RESIDUAL + TLS COMPONENTS                              "
        "   ",
        ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, True)


def test_is_tls_sum_other4():
    """Tests that tls_sum is correctly parsed from REMARK 3 records.

    The example is taken from 2ix9.
    """
    records = {"REMARK": [
        "  3  OTHER REFINEMENT REMARKS: HYDROGENS HAVE BEEN ADDED IN THE       "
        "   ",
        "  3  RIDING POSITIONS.B FACTORS CORRESPOND TO THE OVERALL B FACTORS   "
        "   ",
        "  3  EQUAL TO THE RESIDUAL PLUS THE TLS COMPONENT.                    "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, True)

    records = {"REMARK": [
        "  3  OTHER REFINEMENT REMARKS: RESIDUAL PLUS TLS                      "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, True)


def test_is_tls_sum_other5():
    """Tests that tls_sum is correctly parsed from REMARK 3 records.

    The example is taken from 3msq.
    """
    records = {"REMARK": [
        "  3  OTHER REFINEMENT REMARKS:                                        "
        "   ",
        "  3  INCORPORATION. 2. B-FACTORS CONTAIN BOTH TLS AND RESIDUAL        "
        "   ",
        "  3  COMPONENTS.                                                      "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, True)

    records = {"REMARK": [
        "  3  OTHER REFINEMENT REMARKS:                                        "
        "   ",
        "  3  B FACTOR CONTAINS TLS AND RESIDUAL COMPONENTS                    "
        "   ",
        ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, True)


def test_is_tls_sum_other6():
    """Tests that tls_sum is correctly parsed from REMARK 3 records.

    The example is taken from 1ocx.
    """
    records = {"REMARK": [
        "  3  OTHER REFINEMENT REMARKS: HYDROGENS HAVE BEEN ADDED IN THE       "
        "   ",
        "  3  RIDING POSITIONS. INDIVIDUAL B-FACTOR REFINEMENT WAS PRECEEDED   "
        "   ",
        "  3  BY TLS REFINEMENT IN WHICH EACH MONOMER WAS DIVIDED INTO THREE   "
        "   ",
        "  3  TLS GROUPS. THE DEPOSITED STRUCTURE SHOWS ANISOTROPIC B          "
        "   ",
        "  3  FACTORS THAT RESULT FROM THE COMBINATION OF THE TLS COMPONENTS   "
        "   ",
        "  3  WITH THE RESIDUAL INDIVIDUAL B FACTORS. THE CCP4 PROGRAM         "
        "   ",
        "  3  TLSANAL WAS USED TO COMBINE THE TWO.                             "
        "   ",
        "  4                                                                   "
        "   ",
        ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, True)

    records = {"REMARK": [
        "  3  OTHER REFINEMENT REMARKS:                                        "
        "   ",
        "  3  COMBINATION OF TLS COMPONENT WITH RESIDUAL B FACTOR              "
        "   ",
        ]}
    tls_sum = is_tls_sum(records)
    eq_(tls_sum, True)


def test_get_header_and_trailer():
    """Tests that the header and trailer are correctly parsed from a PDB file.

    Test entry ht.pdb:

HEADER    TEST PROTEIN                            00-JAN-00   1ABC
TITLE     TEST
MODEL        1
ATOM      1
ANISOU
SIGUIJ
ATOM      2
ANISOU
SIGUIJ
HETATM
HETATM
ATOM      3
ANISOU
SIGUIJ
TER       4
ENDMDL
MODEL        2
ATOM      1
ANISOU
SIGUIJ
ATOM      2
ANISOU
SIGUIJ
HETATM
HETATM
ATOM      3
ANISOU
SIGUIJ
TER       4
ENDMDL
CONECT
CONECT
MASTER
END

    """
    header, trailer = get_pdb_header_and_trailer("pdbb/tests/pdb/files/ht.pdb")
    head_exp = [
        "HEADER    TEST PROTEIN                            "
        "00-JAN-00   1ABC              ",
        "TITLE     TEST                                    "
        "                              "]
    trai_exp = [
        "CONECT                                            "
        "                              ",
        "CONECT                                            "
        "                              ",
        "MASTER                                            "
        "                              "]
    eq_(header, head_exp)
    eq_(trailer, trai_exp)


@raises(IOError)
def test_get_header_and_trailer_ioerror():
    """Tests that the header and trailer are correctly parsed from a PDB file.

    Test file not found:
    """
    header, trailer = get_pdb_header_and_trailer("ht.pdb")
    eq_(header, [])
    eq_(trailer, [])
