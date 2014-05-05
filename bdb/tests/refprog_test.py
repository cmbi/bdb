from mock import patch
from nose.tools import eq_, ok_, raises

from bdb.check_beq import get_structure
from bdb.refprog import (decide_refprog, decide_refprog_restrain,
                         except_refprog_warn, filter_progs, last_used,
                         is_bdb_includable_refprog, one_of_the_two,
                         parse_refprog, get_refi_data)
from bdb.pdb.parser import parse_pdb_file


BNEQ_MSG = "Not enough B-factors could be reproduced from ANISOU records"

def test_is_bdb_includable_refprog_none():
    refprog = None
    result = is_bdb_includable_refprog(refprog)
    eq_(result, True)


def test_is_bdb_includable_refprog_cns():
    refprog = "CNS"
    result = is_bdb_includable_refprog(refprog)
    eq_(result, True)


def test_is_bdb_includable_refprog_xplor():
    refprog = "X-PLOR"
    result = is_bdb_includable_refprog(refprog)
    eq_(result, True)


def test_is_bdb_includable_refprog_buster():
    refprog = "BUSTER"
    result = is_bdb_includable_refprog(refprog)
    eq_(result, True)


def test_is_bdb_includable_refprog_refmac():
    refprog = "REFMAC"
    result = is_bdb_includable_refprog(refprog)
    eq_(result, True)


def test_is_bdb_includable_refprog_phenix():
    refprog = "PHENIX"
    result = is_bdb_includable_refprog(refprog)
    eq_(result, True)


# REFMAC

def test_decide_refprog_refmac_remediation_residual():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "has_anisou": False,
                "format_vers": 3.30, "b_type": "residual",
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: residual B-factors (wwPDB remediation 2011)"
    expected = (True, False, True, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_remediation_unverified():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "has_anisou": False,
                "format_vers": 3.30, "b_type": "unverified",
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: B-factor type could not be determined "\
          "(wwPDB remediation 2011)"
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_remediation_full():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "has_anisou": False,
                "format_vers": 3.30, "b_type": None,
                "tls_residual": False, "tls_sum": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: full B-factors (wwPDB remediation 2011)"
    expected = (True, True, False, msg)
    eq_(result, expected)

    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "has_anisou": False,
                "format_vers": 3.30, "b_type": "something",
                "tls_residual": False, "tls_sum": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "beq_identical": 0.9,
                "format_vers": 3.30, "b_type": None, "has_anisou": True,
                "tls_residual": False, "tls_sum": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_refmac_residual():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "has_anisou": False,
                "format_vers": "not3.30", "b_type": "residual",
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: residual B-factors (wwPDB remediation)"
    expected = (True, False, True, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_unverified():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "has_anisou": False,
                "format_vers": "not3.30", "b_type": "unverified",
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: B-factor type could not be determined "\
          "(wwPDB remediation)"
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_residual_and_full():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "has_anisou": False,
                "format_vers": "not3.30", "b_type": None,
                "tls_residual": True, "tls_sum": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: residual and full B-factors (REMARK 3)"
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_tls_residual_anisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": True, "tls_sum": False, "has_anisou": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: TLS group(s), residual B-factors (REMARK 3) "\
          "and ANISOU records. {}".format(BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_tls_residual_noanisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": True, "tls_sum": False, "has_anisou": False,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: TLS group(s), residual B-factors (REMARK 3) "\
          "without ANISOU records"
    expected = (True, False, True, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_tls_sum_anisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": False, "tls_sum": True, "has_anisou": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: TLS group(s), full B-factors (REMARK 3) "\
          "and ANISOU records. {}".format(BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_tls_sum_noanisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": False, "tls_sum": True, "has_anisou": False,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: TLS group(s), full B-factors (REMARK 3) "\
          "without ANISOU records"
    expected = (True, True, False, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_tls_bexcept_noanisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "has_anisou": False,
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": False, "tls_sum": False,
                "other_refinement_remarks": "RESIDUAL U VALUES"}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: TLS group(s) and, possibly, residual or full B-factors "\
          "(REMARK 3, unrecognized format)"
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "B-FACTORS ARE RESIDUALS "
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "B-FACTORS SUM  "
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_refmac_tls_bexcept_anisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "has_anisou": True,
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": False, "tls_sum": False, "beq_identical": 0.9,
                "other_refinement_remarks": "RESIDUAL U VALUES"}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: TLS group(s) and, possibly, residual or full B-factors "\
          "(REMARK 3, unrecognized format)"
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "B-FACTORS ARE RESIDUALS "
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "B-FACTORS SUM  "
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_refmac_tls_nohints_anisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "has_anisou": True,
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": False, "tls_sum": False, "beq_identical": 0.9,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: TLS group(s), no B-factor type details (REMARK 3). {}".\
            format(BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "MAXIMUM LIKELIHOOD RESIDUALS"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL FEATURES"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL ELECTRONS"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL ELECTRON DENSITY"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL DENSITY "
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_refmac_tls_nohints_noanisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "has_anisou": False,
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": False, "tls_sum": False,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: TLS group(s), no B-factor type details (REMARK 3) and no "\
          "ANISOU records"
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "MAXIMUM LIKELIHOOD RESIDUALS"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL FEATURES"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL ELECTRONS"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL ELECTRON DENSITY"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL DENSITY "
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_refmac_notls_tlsremark_noanisou_nohints():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": False, "has_anisou": False,
                "other_refinement_remarks": "SOMETHING MENTIONING TLS..."}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: TLS remark without TLS group(s)"
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_notls_tlsremark_anisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": False, "has_anisou": True,
                "other_refinement_remarks": "SOMETHING MENTIONING TLS..."}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: TLS remark without TLS group(s). {}".format(BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_notls_tlsremark_noanisou_sum():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": True, "has_anisou": False,
                "other_refinement_remarks": "SOMETHING MENTIONING TLS..."}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: TLS remark without TLS group(s). Full B-factors (REMARK 3)"
    expected = (True, True, False, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_notls_tlsremark_anisou_sum():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": True, "has_anisou": True,
                "other_refinement_remarks": "SOMETHING MENTIONING TLS..."}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: TLS remark without TLS group(s). {}. Full B-factors "\
          "(REMARK 3)".format(BNEQ_MSG)
    expected = (True, True, False, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_notls_tlsremark_noanisou_residual():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": True, "tls_sum": False, "has_anisou": False,
                "other_refinement_remarks": "SOMETHING MENTIONING TLS..."}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: TLS remark without TLS group(s). Residual B-factors "\
          "(REMARK 3)"
    expected = (True, False, True, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_notls_tlsremark_anisou_residual():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": True, "tls_sum": False, "has_anisou": True,
                "other_refinement_remarks": "SOMETHING MENTIONING TLS..."}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: TLS remark without TLS group(s). {}. Residual B-factors "\
          "(REMARK 3)".format(BNEQ_MSG)
    expected = (True, False, True, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_notls_residual_anisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": True, "tls_sum": False, "has_anisou": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: residual B-factors without TLS group(s) (REMARK 3). {}".\
            format(BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_notls_residual_noanisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": True, "tls_sum": False, "has_anisou": False,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: residual B-factors without TLS group(s) (REMARK 3)"
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_notls_sum_anisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": True, "has_anisou": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: full B-factors without TLS group(s) (REMARK 3). {}".format(
            BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_notls_sum_noanisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": True, "has_anisou": False,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: full B-factors without TLS group(s) (REMARK 3)"
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_refmac_notls_bexcept_anisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": False, "has_anisou": True,
                "other_refinement_remarks": "RESIDUAL B-FACTOR"}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: possibly, residual or full B-factors (REMARK 3, "\
          "unrecognized format). No TLS groups. {}".format(BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUALS U-VALUES"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = " SUM  AND U-VALUES"
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_refmac_notls_bexcept_noanisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "has_anisou": False,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": False,
                "other_refinement_remarks": "RESIDUAL B-FACTOR"}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: possibly, residual or full B-factors (REMARK 3, "\
          "unrecognized format). No TLS groups"
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUALS U-VALUES"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = " SUM  AND U-VALUES "
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_refmac_notls_nohints_anisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test", "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": False, "has_anisou": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: probably full/mixed anisotropic refinement. {}".\
            format(BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUALS"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUALS"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = " SUM  "
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_refmac_notls_nohints_noanisou():
    pdb_info = {"prog_last": ["REFMAC"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": False, "has_anisou": False,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "REFMAC: probably full B-factors"
    expected = (True, True, False, msg)
    eq_(result, expected)


# OTHER

def test_decide_refprog_OTHER_remediation_residual():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "has_anisou": False,
                "format_vers": 3.30, "b_type": "residual",
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: residual B-factors (wwPDB remediation 2011)"
    expected = (False, False, True, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_remediation_unverified():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "has_anisou": False,
                "format_vers": 3.30, "b_type": "unverified",
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: B-factor type could not be determined "\
          "(wwPDB remediation 2011)"
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_remediation_full():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "has_anisou": False,
                "format_vers": 3.30, "b_type": "somethingunexpected",
                "tls_residual": False, "tls_sum": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: unexpected B-factor type annotation (wwPDB remediation)"
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "beq_identical": 0.9,
                "format_vers": 3.30, "b_type": "something", "has_anisou": True,
                "tls_residual": False, "tls_sum": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "has_anisou": False,
                "format_vers": 3.31, "b_type": "something",
                "tls_residual": False, "tls_sum": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_OTHER_residual():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "has_anisou": False,
                "format_vers": "not3.30", "b_type": "residual",
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: residual B-factors (wwPDB remediation)"
    expected = (False, False, True, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_unverified():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "has_anisou": False,
                "format_vers": "not3.30", "b_type": "unverified",
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: B-factor type could not be determined "\
          "(wwPDB remediation)"
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_residual_and_full():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "has_anisou": False,
                "format_vers": "not3.30", "b_type": None,
                "tls_residual": True, "tls_sum": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: residual and full B-factors (REMARK 3)"
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_tls_residual_anisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": True, "tls_sum": False, "has_anisou": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: TLS group(s), residual B-factors (REMARK 3) "\
          "and ANISOU records. {}".format(BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_tls_residual_noanisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": True, "tls_sum": False, "has_anisou": False,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: TLS group(s), residual B-factors (REMARK 3) "\
          "without ANISOU records"
    expected = (False, False, True, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_tls_sum_anisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": False, "tls_sum": True, "has_anisou": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: TLS group(s), full B-factors (REMARK 3) "\
          "and ANISOU records. {}".format(BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_tls_sum_noanisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": False, "tls_sum": True, "has_anisou": False,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: TLS group(s), full B-factors (REMARK 3) "\
          "without ANISOU records"
    expected = (False, True, False, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_tls_bexcept_noanisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "has_anisou": False,
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": False, "tls_sum": False,
                "other_refinement_remarks": "RESIDUAL U VALUES"}
    result = decide_refprog(pdb_info)
    msg = "OTHER: TLS group(s) and, possibly, residual or full B-factors "\
          "(REMARK 3, unrecognized format)"
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "B-FACTORS ARE RESIDUALS "
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "B-FACTORS SUM  "
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_OTHER_tls_bexcept_anisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "has_anisou": True,
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": False, "tls_sum": False, "beq_identical": 0.9,
                "other_refinement_remarks": "RESIDUAL U VALUES"}
    result = decide_refprog(pdb_info)
    msg = "OTHER: TLS group(s) and, possibly, residual or full B-factors "\
          "(REMARK 3, unrecognized format)"
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "B-FACTORS ARE RESIDUALS "
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "B-FACTORS SUM  "
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_OTHER_tls_nohints_anisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "has_anisou": True,
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": False, "tls_sum": False, "beq_identical": 0.9,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: TLS group(s), no B-factor type details (REMARK 3). {}".\
            format(BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "MAXIMUM LIKELIHOOD RESIDUALS"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL FEATURES"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL ELECTRONS"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL ELECTRON DENSITY"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL DENSITY "
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_OTHER_tls_nohints_noanisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "has_anisou": False,
                "format_vers": "not3.30", "b_type": None, "tls_groups": 1,
                "tls_residual": False, "tls_sum": False,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: TLS group(s), no B-factor type details (REMARK 3) and no "\
          "ANISOU records"
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "MAXIMUM LIKELIHOOD RESIDUALS"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL FEATURES"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL ELECTRONS"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL ELECTRON DENSITY"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUAL DENSITY "
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_OTHER_notls_tlsremark_noanisou_nohints():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": False, "has_anisou": False,
                "other_refinement_remarks": "SOMETHING MENTIONING TLS..."}
    result = decide_refprog(pdb_info)
    msg = "OTHER: TLS remark without TLS group(s)"
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_notls_tlsremark_anisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": False, "has_anisou": True,
                "other_refinement_remarks": "SOMETHING MENTIONING TLS..."}
    result = decide_refprog(pdb_info)
    msg = "OTHER: TLS remark without TLS group(s). {}".format(BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_notls_tlsremark_noanisou_sum():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": True, "has_anisou": False,
                "other_refinement_remarks": "SOMETHING MENTIONING TLS..."}
    result = decide_refprog(pdb_info)
    msg = "OTHER: TLS remark without TLS group(s). Full B-factors (REMARK 3)"
    expected = (False, True, False, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_notls_tlsremark_anisou_sum():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": True, "has_anisou": True,
                "other_refinement_remarks": "SOMETHING MENTIONING TLS..."}
    result = decide_refprog(pdb_info)
    msg = "OTHER: TLS remark without TLS group(s). {}. Full B-factors "\
          "(REMARK 3)".format(BNEQ_MSG)
    expected = (False, True, False, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_notls_tlsremark_noanisou_residual():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": True, "tls_sum": False, "has_anisou": False,
                "other_refinement_remarks": "SOMETHING MENTIONING TLS..."}
    result = decide_refprog(pdb_info)
    msg = "OTHER: TLS remark without TLS group(s). Residual B-factors "\
          "(REMARK 3)"
    expected = (False, False, True, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_notls_tlsremark_anisou_residual():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": True, "tls_sum": False, "has_anisou": True,
                "other_refinement_remarks": "SOMETHING MENTIONING TLS..."}
    result = decide_refprog(pdb_info)
    msg = "OTHER: TLS remark without TLS group(s). {}. Residual B-factors "\
          "(REMARK 3)".format(BNEQ_MSG)
    expected = (False, False, True, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_notls_residual_anisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": True, "tls_sum": False, "has_anisou": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: residual B-factors without TLS group(s) (REMARK 3). {}".\
            format(BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_notls_residual_noanisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": True, "tls_sum": False, "has_anisou": False,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: residual B-factors without TLS group(s) (REMARK 3)"
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_notls_sum_anisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": True, "has_anisou": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: full B-factors without TLS group(s) (REMARK 3). {}".format(
            BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_notls_sum_noanisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": True, "has_anisou": False,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: full B-factors without TLS group(s) (REMARK 3)"
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_OTHER_notls_bexcept_anisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": False, "has_anisou": True,
                "other_refinement_remarks": "RESIDUAL B-FACTOR"}
    result = decide_refprog(pdb_info)
    msg = "OTHER: possibly, residual or full B-factors (REMARK 3, "\
          "unrecognized format). No TLS groups. {}".format(BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUALS U-VALUES"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = " SUM  AND U-VALUES"
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_OTHER_notls_bexcept_noanisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "has_anisou": False,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": False,
                "other_refinement_remarks": "RESIDUAL B-FACTOR"}
    result = decide_refprog(pdb_info)
    msg = "OTHER: possibly, residual or full B-factors (REMARK 3, "\
          "unrecognized format). No TLS groups"
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUALS U-VALUES"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = " SUM  AND U-VALUES "
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_OTHER_notls_nohints_anisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test", "beq_identical": 0.9,
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": False, "has_anisou": True,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: probably full/mixed anisotropic refinement. {}".\
            format(BNEQ_MSG)
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUALS"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = "RESIDUALS"
    result = decide_refprog(pdb_info)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = " SUM  "
    result = decide_refprog(pdb_info)
    eq_(result, expected)


def test_decide_refprog_OTHER_notls_nohints_noanisou():
    pdb_info = {"prog_last": ["OTHER"], "pdb_id": "test",
                "format_vers": "not3.30", "b_type": None, "tls_groups": None,
                "tls_residual": False, "tls_sum": False, "has_anisou": False,
                "other_refinement_remarks": ""}
    result = decide_refprog(pdb_info)
    msg = "OTHER: probably full B-factors"
    expected = (True, True, False, msg)
    eq_(result, expected)

    pdb_info["other_refinement_remarks"] = None
    result = decide_refprog(pdb_info)
    eq_(result, expected)


# RESTRAIN

def test_decide_refprog_restrain_u():
    pdb_info = {"prog_last": ["RESTRAIN"],
            "pdb_id": "test",
            "other_refinement_remarks": " THE QUANTITY PRESENTED IN THE "
            "TEMPERATURE FACTOR FIELD IS U.",
            "b_msqav": False,
            "b_type": None,
            "has_anisou": False,
            "tls_groups": None,
            "tls_residual": False,
            "tls_sum": False}
    msg = "RESTRAIN: B-factor field contains \"U\" (REMARK 3)"
    result = decide_refprog_restrain(pdb_info)
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_restrain_uu():
    pdb_info = {"prog_last": ["RESTRAIN"], "pdb_id": "test", "b_msqav": True,
                "other_refinement_remarks": "", "has_anisou": False}
    msg = "RESTRAIN: B-factor field contains mean square atomic displacement "\
          "(REMARK 3)"
    result = decide_refprog_restrain(pdb_info)
    expected = (True, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_restrain_uu_2():
    pdb_info = {"prog_last": ["RESTRAIN"], "pdb_id": "test", "b_msqav": True,
                "other_refinement_remarks": "", "has_anisou": False}
    msg = "RESTRAIN: B-factor field contains mean square atomic displacement "\
          "(REMARK 3)"
    result = decide_refprog(pdb_info)
    expected = (True, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_restrain_nothandled():
    pdb_info_def = {"prog_last": ["RESTRAIN"],
            "pdb_id": "test",
            "beq_identical": 0.9,
            "other_refinement_remarks": "",
            "b_msqav": False,
            "b_type": None,
            "has_anisou": False,
            "tls_groups": None,
            "tls_residual": False,
            "tls_sum": False}

    pdb_info = pdb_info_def.copy()
    pdb_info["b_type"] = "residual"
    result = decide_refprog_restrain(pdb_info)
    msg = "RESTRAIN: unexpected content cannot (yet) be handled"
    expected = (False, False, False, msg)
    eq_(result, expected)

    pdb_info = pdb_info_def.copy()
    pdb_info["has_anisou"] = True
    result = decide_refprog_restrain(pdb_info)
    eq_(result, expected)

    pdb_info = pdb_info_def.copy()
    pdb_info["tls_groups"] = 1
    result = decide_refprog_restrain(pdb_info)
    eq_(result, expected)

    pdb_info = pdb_info_def.copy()
    pdb_info["tls_residual"] = True
    result = decide_refprog_restrain(pdb_info)
    eq_(result, expected)

    pdb_info = pdb_info_def.copy()
    pdb_info["tls_sum"] = True
    result = decide_refprog_restrain(pdb_info)
    eq_(result, expected)


def test_decide_refprog_restrain_useful():
    pdb_info_def = {"prog_last": ["RESTRAIN"],
            "pdb_id": "test",
            "other_refinement_remarks": "",
            "b_msqav": False,
            "b_type": None,
            "has_anisou": False,
            "tls_groups": None,
            "tls_residual": False,
            "tls_sum": False}
    result = decide_refprog_restrain(pdb_info_def)
    msg = "RESTRAIN: probably full B-factors"
    expected = (True, True, False, msg)
    eq_(result, expected)


@raises(AssertionError)
def test_decide_refprog_type_error_list():
    pdb_info = {"prog_last": "notalist", "has_anisou": False, "pdb_id": "test"}
    decide_refprog(pdb_info)

    pdb_info = {"prog_last": 1}
    decide_refprog(pdb_info)


@raises(AssertionError)
def test_decide_refprog_type_error_string():
    pdb_info = {"prog_last": [1], "has_anisou": False, "pdb_id": "test"}
    decide_refprog(pdb_info)


def test_decide_refprog_too_many():
    pdb_info = {"prog_last": ["1", "2"], "has_anisou": False, "pdb_id": "test"}
    result = decide_refprog(pdb_info)
    msg = "Combination of refinement programs cannot (yet) be included "\
          "in the bdb: 1 and 2"
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_decide_refprog_too_few():
    pdb_info = {"prog_last": [], "has_anisou": False, "pdb_id": "test"}
    result = decide_refprog(pdb_info)
    msg = "Program(s) in REMARK 3 not interpreted as refinement program(s)"
    expected = (False, False, False, msg)
    eq_(result, expected)


@patch("bdb.refprog.is_bdb_includable_refprog", return_value=False)
def test_decide_refprog_nonincludable_refprog(*args):
    pdb_info = {"prog_last": ["A"], "has_anisou": False, "pdb_id": "test"}
    result = decide_refprog(pdb_info)
    msg = "A: this program cannot (yet) be included"
    expected = (False, False, False, msg)
    eq_(result, expected)


def test_except_refprog_warn():
    eq_(except_refprog_warn(), None)


def test_filter_progs_empty():
    pin = []
    pv = []
    result = filter_progs(pin, pv)
    eq_(result, ([],[]))


@raises(AssertionError)
def test_filter_progs_type_error():
    pin = "program"
    pv = []
    filter_progs(pin, pv)

    pin = []
    pv = "version"
    filter_progs(pin, pv)


@raises(AssertionError)
def test_filter_progs_different_length():
    pin = []
    pv = ["version"]
    filter_progs(pin, pv)

    pin = ["program"]
    pv = []
    filter_progs(pin, pv)


def test_filter_progs_no_filter():
    pin = ["REFMAC", "SHELX"]
    pv = ["5.6", "L"]
    result = filter_progs(pin, pv)
    eq_(result, (pin, pv))


def test_filter_progs_filtered():
    pin = ["REFMAC", "COOT"]
    pv = ["5.7", "-"]
    result = filter_progs(pin, pv)
    expected = (["REFMAC"], ["5.7"])
    eq_(result, expected)

    pin = ["SOLVE/RESOLVE", "TOM/FRODO", "XTALVIEW"]
    pv = ["-", "-", "-"]
    result = filter_progs(pin, pv)
    expected = ([], [])
    eq_(result, expected)

    pop_progs = ["ARP", "ARP/WARP", "CHAIN", "COOT", "DM", "FRODO", "HKL-3000",
                 "LAFIRE", "MOLPROBITY", "O", "OOPS", "PIKSOL", "PRODRG",
                 "PROTEIN", "QUANTA", "SCWRL", "SFALL", "SOLVE/RESOLVE", "TOM",
                 "TOM/FRODO", "XFIT", "XPLEO", "XTALVIEW"]
    pv = ["-"]
    for p in pop_progs:
        result = filter_progs([p], pv)
        eq_(result, ([], []))


def test_get_refi_data_1crn():
    pdb_file_path = "bdb/tests/pdb/files/1crn.pdb"
    pdb_id = "1crn"
    structure = get_structure(pdb_file_path, pdb_id)
    records = parse_pdb_file(pdb_file_path)
    pdb_info = get_refi_data(records, structure, pdb_id)
    eq_(pdb_info["assume_iso"], True)
    eq_(pdb_info["b_msqav"], False)
    eq_(pdb_info["b_type"], None)
    eq_(pdb_info["beq_identical"], None)
    eq_(pdb_info["correct_uij"], None)
    eq_(pdb_info["decision"], "PROLSQ: probably full B-factors")
    eq_(pdb_info["format_date"], "13-JUL-11")
    eq_(pdb_info["format_vers"], 3.3)
    eq_(pdb_info["has_anisou"], False)
    eq_(pdb_info["is_bdb_includable"], True)
    eq_(pdb_info["other_refinement_remarks"], "")
    eq_(pdb_info["pdb_id"], "1crn")
    eq_(pdb_info["prog_inter"], ["PROLSQ"])
    eq_(pdb_info["prog_last"], ["PROLSQ"])
    eq_(pdb_info["prog_vers"], ["PROLSQ"])
    eq_(pdb_info["ref_prog"], ["PROLSQ"])
    eq_(pdb_info["refprog"], "PROLSQ")
    eq_(pdb_info["req_tlsanl"], False)
    eq_(pdb_info["tls_groups"], None)
    eq_(pdb_info["tls_residual"], False)
    eq_(pdb_info["tls_sum"], False)


def test_get_refi_data_3zzw():
    pdb_file_path = "bdb/tests/pdb/files/3zzw.pdb"
    pdb_id = "3zzw"
    structure = get_structure(pdb_file_path, pdb_id)
    records = parse_pdb_file(pdb_file_path)
    pdb_info = get_refi_data(records, structure, pdb_id)
    eq_(pdb_info["assume_iso"], True)
    eq_(pdb_info["b_msqav"], False)
    eq_(pdb_info["b_type"], None)
    eq_(pdb_info["beq_identical"], 1.0)
    eq_(pdb_info["correct_uij"], True)
    eq_(pdb_info["decision"], "Assuming full isotropic B-factors because "\
            "enough B-factors can be reproduced from the ANISOU records")
    eq_(pdb_info["format_date"], "13-JUL-11")
    eq_(pdb_info["format_vers"], 3.3)
    eq_(pdb_info["has_anisou"], True)
    eq_(pdb_info["is_bdb_includable"], True)
    eq_(pdb_info["other_refinement_remarks"], "IDEAL-DIST CONTACT TERM "\
            "CONTACT SETUP.  ALL ATOMS HAVE CCP4 ATOM TYPE FROM LIBRARY.")
    eq_(pdb_info["pdb_id"], "3zzw")
    eq_(pdb_info["prog_inter"], ["BUSTER"])
    eq_(pdb_info["prog_last"], ["BUSTER"])
    eq_(pdb_info["prog_vers"], ["2.11.1"])
    eq_(pdb_info["ref_prog"], ["BUSTER 2.11.1"])
    eq_(pdb_info["refprog"], "BUSTER 2.11.1")
    eq_(pdb_info["req_tlsanl"], False)
    eq_(pdb_info["tls_groups"], 8)
    eq_(pdb_info["tls_residual"], False)
    eq_(pdb_info["tls_sum"], False)


def test_get_refi_data_2wnl():
    pdb_file_path = "bdb/tests/pdb/files/2wnl.pdb"
    pdb_id = "2wnl"
    structure = get_structure(pdb_file_path, pdb_id)
    records = parse_pdb_file(pdb_file_path)
    pdb_info = get_refi_data(records, structure, pdb_id)
    eq_(pdb_info["assume_iso"], True)
    eq_(pdb_info["b_msqav"], False)
    eq_(pdb_info["b_type"], None)
    eq_(pdb_info["beq_identical"], None)
    eq_(pdb_info["correct_uij"], None)
    eq_(pdb_info["decision"], "REFMAC: TLS group(s), full B-factors "\
                              "(REMARK 3) without ANISOU records")
    eq_(pdb_info["format_date"], "01-DEC-08")
    eq_(pdb_info["format_vers"], 3.2)
    eq_(pdb_info["has_anisou"], False)
    eq_(pdb_info["is_bdb_includable"], True)
    eq_(pdb_info["other_refinement_remarks"], "HYDROGENS HAVE BEEN ADDED IN "\
            "THE  RIDING POSITIONS. GLOBAL B-FACTORS, CONTAINING RESIDUAL  "\
            "AND TLS COMPONENT HAVE BEEN DEPOSITED.")
    eq_(pdb_info["pdb_id"], "2wnl")
    eq_(pdb_info["prog_inter"], ["REFMAC"])
    eq_(pdb_info["prog_last"], ["REFMAC"])
    eq_(pdb_info["prog_vers"], ["5.5.0102"])
    eq_(pdb_info["ref_prog"], ["REFMAC 5.5.0102"])
    eq_(pdb_info["refprog"], "REFMAC 5.5.0102")
    eq_(pdb_info["req_tlsanl"], False)
    eq_(pdb_info["tls_groups"], 10)
    eq_(pdb_info["tls_residual"], False)
    eq_(pdb_info["tls_sum"], True)


# TODO fix patch
#@patch("bdb.pdb.parser.parse_ref_prog", return_value=None)
#def test_get_refprog_no_refprog(*args):
#    pdb_file_path = "bdb/tests/pdb/files/1crn.pdb"
#    pdb_id = "1crn"
#    structure = get_structure(pdb_file_path, pdb_id)
#    records = parse_pdb_file(pdb_file_path)
#    pdb_info = get_refi_data(records, structure, pdb_id)
#    msg = "Refinement program parse error"
#    eq_(pdb_info["refprog"], None)
#    eq_(pdb_info["decision"], msg)


def test_last_used_empty():
    pin = []
    pv = []
    result = last_used(pin, pv)
    eq_(result, [])


@raises(AssertionError)
def test_last_used_type_error():
    pin = "program"
    pv = []
    last_used(pin, pv)

    pin = []
    pv = "version"
    last_used(pin, pv)


@raises(AssertionError)
def test_last_used_different_length():
    pin = []
    pv = ["version"]
    last_used(pin, pv)

    pin = ["program"]
    pv = []
    last_used(pin, pv)


def test_last_used_filtered():
    pin = ["REFMAC", "COOT"]
    pv = ["5.7", "-"]
    result = last_used(pin, pv)
    expected = ["REFMAC"]
    eq_(result, expected)

    pin = ["SOLVE/RESOLVE", "TOM/FRODO", "XTALVIEW"]
    pv = ["-", "-", "-"]
    result = last_used(pin, pv)
    expected = []
    eq_(result, expected)


def test_last_used_shelx():
    pin = ["REFMAC", "SHELX"]
    pv = ["5.8", "L"]
    result = last_used(pin, pv)
    expected = ["SHELX"]
    eq_(result, expected)

    pin = ["REFMAC", "SHELX"]
    pv = ["5.7", "H"]
    result = last_used(pin, pv)
    expected = ["SHELX"]
    eq_(result, expected)

    pin = ["REFMAC", "SHELX"]
    pv = ["5.6", "-"]
    result = last_used(pin, pv)
    expected = ["REFMAC"]
    eq_(result, expected)

    pin = ["REFMAC", "SHELX"]
    pv = ["5.5", "np"]
    result = last_used(pin, pv)
    expected = ["REFMAC"]
    eq_(result, expected)

    pin = ["SHELX"]
    pv = ["-"]
    result = last_used(pin, pv)
    eq_(result, pin)

    pin = ["SHELX"]
    pv = ["L"]
    result = last_used(pin, pv)
    eq_(result, pin)


def test_last_used_same_refprogs():
    # e.g. 1c4k
    pin = ["X-PLOR", "X-PLOR"]
    pv = ["3.1", "3.8"]
    result = last_used(pin, pv)
    expected = ["X-PLOR"]
    eq_(result, expected)


def test_last_used_corels():
    pin = ["REFMAC", "CORELS"]
    pv = ["5.8", "-"]
    result = last_used(pin, pv)
    expected = ["REFMAC"]
    eq_(result, expected)

    pin = ["CORELS"]
    pv = ["-"]
    result = last_used(pin, pv)
    eq_(result, pin)


def test_last_used_known_combis():
    # order doesn't matter
    pin = ["CNS", "CCP4"]
    pv = ["-", "-"]
    result = last_used(pin, pv)
    expected = ["CNS"]
    eq_(result, expected)

    known_combis = [("CCP4", "CNS"), ("CCP4", "X-PLOR"), ("CNS", "BUSTER"),
                    ("CNS", "PHENIX.REFINE"), ("CNS", "REFMAC"),
                    ("CNS", "TNT"), ("CNX", "BUSTER"),
                    ("CNX", "PHENIX.REFINE"), ("CNX", "REFMAC"),
                    ("CNX", "TNT"), ("EREF", "PROLSQ"), ("EREF", "X-PLOR"),
                    ("PROLSQ", "TNT"), ("PROLSQ", "REFMAC"),
                    ("REFMAC", "MAIN"), ("RESTRAIN", "X-PLOR"),
                    ("TWIN_LSQ", "CNS"), ("TNT", "REFMAC"),
                    ("X-PLOR", "BUSTER"), ("X-PLOR", "CNS"),
                    ("X-PLOR", "PHENIX.REFINE"), ("X-PLOR", "PROLSQ"),
                    ("X-PLOR", "REFMAC"), ("X-PLOR", "TNT"), ]
    for c in known_combis:
        result = last_used([c[0], c[1]], pv)
        expected = [c[1]]
        eq_(result, expected)



def test_last_used_unknown_combis():
    pin = ["X-PLOR", "CNS", "REFMAC", "PHENIX"]
    pv = ["3.85", "1.3", "5.8", "-"]
    result = last_used(pin, pv)
    expected = ["PHENIX", "REFMAC"]
    eq_(result, expected)


def test_last_used_sort():
    pin = ["B", "A"]
    pv = ["-", "-"]
    result = last_used(pin, pv)
    expected = ["A", "B"]
    eq_(result, expected)


def test_one_of_the_two():
    lst = ["A", "B"]
    loser = "A"
    winner = "B"
    result = one_of_the_two(lst, loser, winner)
    expected = ["B"]
    eq_(result, expected)


@raises(AssertionError)
def test_one_of_the_two_type_error():
    lst = "A"
    loser = "A"
    winner = "B"
    one_of_the_two(lst, loser, winner)


def test_one_of_the_two_fail():
    lst = ["A", "B"]
    loser = "A"
    winner = "C"
    result = one_of_the_two(lst, loser, winner)
    expected = ["A", "B"]
    eq_(result, expected)

    lst = ["A", "B"]
    loser = "C"
    winner = "B"
    result = one_of_the_two(lst, loser, winner)
    expected = ["A", "B"]
    eq_(result, expected)

    lst = ["A", "B"]
    loser = "B"
    winner = "B"
    result = one_of_the_two(lst, loser, winner)
    expected = ["A", "B"]
    eq_(result, expected)


def test_parse_refprog_null():
    refprog = "NULL"
    result = parse_refprog(refprog)
    expected = ([None], [None], [None])
    eq_(result, expected)


def test_parse_refprog_predefined_exceptions():
    refprog = "X-PLOR 3.1 AND 3.85"
    result = parse_refprog(refprog)
    expected =  (["X-PLOR 3.85"], ["X-PLOR"], ["3.85"])
    eq_(result, expected)

    refprog = "X-PLOR 3.1, 3.816"
    result = parse_refprog(refprog)
    expected =  (["X-PLOR 3.816"], ["X-PLOR"], ["3.816"])
    eq_(result, expected)

    refprog = "CNS 1.1 & 1.3"
    result = parse_refprog(refprog)
    expected =  (["CNS 1.3"], ["CNS"], ["1.3"])
    eq_(result, expected)

    refprog = "CNS 0.4, O, OOPS"
    result = parse_refprog(refprog)
    expected =  (["CNS 0.4", "O", "OOPS"],
                 ["CNS", "O", "OOPS"],
                 ["0.4", "-", "-"])
    eq_(result, expected)

    refprog = "CNS 0.1-0.4"
    result = parse_refprog(refprog)
    expected =  (["CNS 0.4"], ["CNS"], ["0.4"])
    eq_(result, expected)

    refprog = "CNS 0.9,1.0,1.1"
    result = parse_refprog(refprog)
    expected =  (["CNS 1.1"], ["CNS"], ["1.1"])
    eq_(result, expected)

    refprog = "CNS 1.3 WITH DEN REFINEMENT"
    result = parse_refprog(refprog)
    expected =  (["CNS 1.3"], ["CNS"], ["1.3"])
    eq_(result, expected)

    refprog = "CNS 1.2 (USING XTAL_TWIN UTILITIES)"
    result = parse_refprog(refprog)
    expected =  (["CNS 1.2"], ["CNS"], ["1.2"])
    eq_(result, expected)

    refprog = "PHENIX.REFINE_REFMAC 5.5.0070"
    result = parse_refprog(refprog)
    expected =  (["PHENIX.REFINE", "REFMAC 5.5.0070"],
                 ["PHENIX.REFINE", "REFMAC"],
                 ["-", "5.5.0070"])
    eq_(result, expected)

    refprog = "PHENIX (CCI APPS 2007_04_06_1210)"
    result = parse_refprog(refprog)
    expected =  (["PHENIX (PHENIX.REFINE: 2007_04_06_1210)"],
                 ["PHENIX.REFINE"],
                 ["2007_04_06_1210"])
    eq_(result, expected)

    refprog = "PHENIX VERSION 1.8_1069 (PHENIX.REFINE)"
    result = parse_refprog(refprog)
    expected =  (["PHENIX (PHENIX.REFINE: 1.8_1069)"], ["PHENIX.REFINE"],
                 ["1.8_1069"])
    eq_(result, expected)

    refprog = "PHENIX 1.6.2_432 - REFINE"
    result = parse_refprog(refprog)
    expected =  (["PHENIX (PHENIX.REFINE: 1.6.2_432)"], ["PHENIX.REFINE"],
                 ["1.6.2_432"])
    eq_(result, expected)

    refprog = "PHENIX REFINE"
    result = parse_refprog(refprog)
    expected = (["PHENIX.REFINE"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX, REFINE"
    result = parse_refprog(refprog)
    expected = (["PHENIX.REFINE"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX AUTOREFINE"
    result = parse_refprog(refprog)
    expected = (["PHENIX.REFINE"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "REFMAC 5.1.24/TLS"
    result = parse_refprog(refprog)
    expected =  (["REFMAC 5.1.24"], ["REFMAC"], ["5.1.24"])
    eq_(result, expected)

    refprog = "REFMAC 5.2.0005 24/04/2001"
    result = parse_refprog(refprog)
    expected =  (["REFMAC 5.2.0005"], ["REFMAC"], ["5.2.0005"])
    eq_(result, expected)

    refprog = "REFMAC 5.2.0019 24/04/2001"
    result = parse_refprog(refprog)
    expected =  (["REFMAC 5.2.0019"], ["REFMAC"], ["5.2.0019"])
    eq_(result, expected)

    refprog = "REFMAC X-PLOR 3.843"
    result = parse_refprog(refprog)
    expected =  (["REFMAC", "X-PLOR 3.843"],
                 ["REFMAC", "X-PLOR"],
                 ["-", "3.843"])
    eq_(result, expected)

    refprog = "REFMAC5 5.2.0019"
    result = parse_refprog(refprog)
    expected =  (["REFMAC 5.2.0019"], ["REFMAC"], ["5.2.0019"])
    eq_(result, expected)

    refprog = "REFMAC 5.5.0109 (AND PHENIX)"
    result = parse_refprog(refprog)
    expected =  (["REFMAC 5.5.0109", "PHENIX.REFINE"],
                 ["REFMAC", "PHENIX.REFINE"],
                 ["5.5.0109", "-"])
    eq_(result, expected)

    refprog = "BUSTER, BETA VERSION"
    result = parse_refprog(refprog)
    expected =  (["BUSTER BETA"], ["BUSTER"], ["BETA"])
    eq_(result, expected)

    refprog = "TNT BUSTER/TNT"
    result = parse_refprog(refprog)
    expected =  (["TNT", "BUSTER/TNT"], ["TNT", "BUSTER"], ["-", "-"])
    eq_(result, expected)

    refprog = "O, VERSION 9.0.7"
    result = parse_refprog(refprog)
    expected =  (["O 9.0.7"], ["O"], ["9.0.7"])
    eq_(result, expected)


def test_parse_refprog_nosplit():
    # 4ow3
    refprog = "PHENIX (PHENIX.REFINE: DEV_1549+SVN)      "
    result = parse_refprog(refprog)
    expected =  (["PHENIX (PHENIX.REFINE: DEV_1549+SVN)"], ["PHENIX.REFINE"],
                 ["np"])
    eq_(result, expected)

    refprog = "BUSTER/TNT"
    result = parse_refprog(refprog)
    expected =  (["BUSTER/TNT"], ["BUSTER"], ["-"])
    eq_(result, expected)

    refprog = "BUSTER/RESOLVE"
    result = parse_refprog(refprog)
    expected =  (["BUSTER/RESOLVE"], ["BUSTER"], ["np"])
    eq_(result, expected)

    refprog = "ARP/WARP"
    result = parse_refprog(refprog)
    expected =  (["ARP/WARP"], ["ARP/WARP"], ["-"])
    eq_(result, expected)

    refprog = "WARP/ARP"
    result = parse_refprog(refprog)
    expected =  (["WARP/ARP"], ["OTHER"], ["np"])
    eq_(result, expected)


def test_parse_refprog_split():
    refprog = "TNT/BUSTER"
    result = parse_refprog(refprog)
    expected =  (["TNT", "BUSTER"], ["TNT", "BUSTER"], ["-", "-"])
    eq_(result, expected)

    refprog = "TNT, BUSTER"
    result = parse_refprog(refprog)
    expected =  (["TNT", "BUSTER"], ["TNT", "BUSTER"], ["-", "-"])
    eq_(result, expected)

    refprog = "TNT&BUSTER"
    result = parse_refprog(refprog)
    expected =  (["TNT", "BUSTER"], ["TNT", "BUSTER"], ["-", "-"])
    eq_(result, expected)

    refprog = "TNT ; BUSTER"
    result = parse_refprog(refprog)
    expected =  (["TNT", "BUSTER"], ["TNT", "BUSTER"], ["-", "-"])
    eq_(result, expected)

    refprog = "TNT +, BUSTER"
    result = parse_refprog(refprog)
    expected =  (["TNT", "BUSTER"], ["TNT", "BUSTER"], ["-", "-"])
    eq_(result, expected)

    refprog = "TNT AND BUSTER,"
    result = parse_refprog(refprog)
    expected =  (["TNT", "BUSTER"], ["TNT", "BUSTER"], ["-", "-"])
    eq_(result, expected)

    refprog = "/REFMAC5.7/////COOT//"
    result = parse_refprog(refprog)
    expected =  (["REFMAC5.7", "COOT"], ["REFMAC", "COOT"], ["5.7", "-"])
    eq_(result, expected)


def test_parse_refprog_no_versions():
    refprogs = ["ARP", "BILDER", "CCP4", "CEDAR", "CHAIN", "CORELS", "DM",
                "DYNAMIX", "EREF", "FFX", "FMLS/VP", "FRODO", "HIPHOP",
                "IMPLOR", "LAFIRE", "LALS", "MAIN", "MOLPROBITY", "MOLLY",
                "MOPRO", "NCNS", "NCNS-TINKER", "NMREF", "PIKSOL", "PHASER",
                "PMB", "POLYVISION", "PRIMEX", "PRODRG", "PROTEIN", "QUANTA",
                "RESTRAIN", "SCWRL", "SHARP", "SFALL", "SOLVE",
                "SOLVE/RESOLVE", "TIBBITTS", "TOM", "TOM/FRODO", "XFIT",
                "XPLEO", "XTALVIEW"]
    for refprog in refprogs:
        result = parse_refprog(refprog)
        expected = ([refprog], [refprog], ["-"])
        eq_(result, expected)


def test_parse_refprog_no_versions_with_version():
    refprog = "CCP4 VERSION"
    result = parse_refprog(refprog)
    expected = ([refprog], ["CCP4"], ["np"])
    eq_(result, expected)


def test_parse_refprog_refmac():
    refprog = "refamc"
    result = parse_refprog(refprog)
    expected = (["REFAMC"], ["REFMAC"], ["-"])
    eq_(result, expected)

    refprog = "REFMEC"
    result = parse_refprog(refprog)
    expected = (["REFMEC"], ["REFMAC"], ["-"])
    eq_(result, expected)

    # Theoretically possible...
    refprog = "CCCCCC"
    result = parse_refprog(refprog)
    expected = (["CCCCCC"], ["REFMAC"], ["-"])
    eq_(result, expected)

    refprog = "REFMAC"
    result = parse_refprog(refprog)
    expected = (["REFMAC"], ["REFMAC"], ["-"])
    eq_(result, expected)

    refprog = "REFMAC 5.7"
    result = parse_refprog(refprog)
    expected = (["REFMAC 5.7"], ["REFMAC"], ["5.7"])
    eq_(result, expected)

    refprog = "REFMAC5.7"
    result = parse_refprog(refprog)
    expected = (["REFMAC5.7"], ["REFMAC"], ["5.7"])
    eq_(result, expected)

    refprog = "REFMAC V 5.7"
    result = parse_refprog(refprog)
    expected = (["REFMAC V 5.7"], ["REFMAC"], ["5.7"])
    eq_(result, expected)

    refprog = "REFMAC V. 5.7"
    result = parse_refprog(refprog)
    expected = (["REFMAC V. 5.7"], ["REFMAC"], ["np"])
    eq_(result, expected)

    refprog = "REFMAC NONE"
    result = parse_refprog(refprog)
    expected = (["REFMAC NONE"], ["REFMAC"], ["NONE"])
    eq_(result, expected)

    refprog = "REFMAC                   5.8                  "
    result = parse_refprog(refprog)
    expected = (["REFMAC                   5.8"], ["REFMAC"], ["np"])
    eq_(result, expected)

    refprog = "REFMAC 5.8                  "
    result = parse_refprog(refprog)
    expected = (["REFMAC 5.8"], ["REFMAC"], ["5.8"])
    eq_(result, expected)


def test_parse_refprog_cns():
    refprog = "CNS"
    result = parse_refprog(refprog)
    expected = (["CNS"], ["CNS"], ["-"])
    eq_(result, expected)

    refprog = "CNS          "
    result = parse_refprog(refprog)
    expected = (["CNS"], ["CNS"], ["-"])
    eq_(result, expected)

    refprog = "CNS (2)"
    result = parse_refprog(refprog)
    expected = (["CNS (2)"], ["CNS"], ["2"])
    eq_(result, expected)

    refprog = "CNS ( 2 )"
    result = parse_refprog(refprog)
    expected = (["CNS ( 2 )"], ["CNS"], ["np"])
    eq_(result, expected)

    refprog = "CNS 2"
    result = parse_refprog(refprog)
    expected = (["CNS 2"], ["CNS"], ["2"])
    eq_(result, expected)

    refprog = "CNS2.0"
    result = parse_refprog(refprog)
    expected = (["CNS2.0"], ["CNS"], ["2.0"])
    eq_(result, expected)

    refprog = "CNS V. ABC123.456"
    result = parse_refprog(refprog)
    expected = (["CNS V. ABC123.456"], ["CNS"], ["ABC123.456"])
    eq_(result, expected)

    refprog = "CNS-3"
    result = parse_refprog(refprog)
    expected = (["CNS-3"], ["CNS"], ["3"])
    eq_(result, expected)

    # Theoretically possible...
    refprog = "CNS-V.3"
    result = parse_refprog(refprog)
    expected = (["CNS-V.3"], ["CNS"], ["V.3"])
    eq_(result, expected)

    # Theoretically possible...
    refprog = "CNS (V.3)"
    result = parse_refprog(refprog)
    expected = (["CNS (V.3)"], ["CNS"], ["V.3"])
    eq_(result, expected)

    # Theoretically possible...
    refprog = "CNS-"
    result = parse_refprog(refprog)
    expected = (["CNS-"], ["CNS"], ["np"])
    eq_(result, expected)


def test_parse_refprog_twinlsq():
    refprog = "TWIN_LSQ"
    result = parse_refprog(refprog)
    expected = (["TWIN_LSQ"], ["CNS"], ["TWIN_LSQ"])
    eq_(result, expected)

    refprog = "TWIN_LSQ 4"
    result = parse_refprog(refprog)
    expected = (["TWIN_LSQ 4"], ["CNS"], ["np"])
    eq_(result, expected)


def test_parse_refprog_cnx():
    refprog = "CNX"
    result = parse_refprog(refprog)
    expected = (["CNX"], ["CNX"], ["-"])
    eq_(result, expected)

    refprog = "CNX          "
    result = parse_refprog(refprog)
    expected = (["CNX"], ["CNX"], ["-"])
    eq_(result, expected)

    refprog = "CNX (ACCELRYS)"
    result = parse_refprog(refprog)
    expected = (["CNX (ACCELRYS)"], ["CNX"], ["-"])
    eq_(result, expected)

    refprog = "CNX V0"
    result = parse_refprog(refprog)
    expected = (["CNX V0"], ["CNX"], ["np"])
    eq_(result, expected)

    refprog = "CNX 0"
    result = parse_refprog(refprog)
    expected = (["CNX 0"], ["CNX"], ["0"])
    eq_(result, expected)

    refprog = "CNX 0.9-9"
    result = parse_refprog(refprog)
    expected = (["CNX 0.9-9"], ["CNX"], ["0.9-9"])
    eq_(result, expected)

    refprog = "CNX0.9-9"
    result = parse_refprog(refprog)
    expected = (["CNX0.9-9"], ["CNX"], ["0.9-9"])
    eq_(result, expected)


def test_parse_refprog_xplor():
    refprog = "X-PLOR"
    result = parse_refprog(refprog)
    expected = (["X-PLOR"], ["X-PLOR"], ["-"])
    eq_(result, expected)

    refprog = "XPLOR"
    result = parse_refprog(refprog)
    expected = (["XPLOR"], ["X-PLOR"], ["-"])
    eq_(result, expected)

    refprog = "X-PLOR (ONLINE)"
    result = parse_refprog(refprog)
    expected = (["X-PLOR (ONLINE)"], ["X-PLOR"], ["-"])
    eq_(result, expected)

    refprog = "X-PLOR V5.4321"
    result = parse_refprog(refprog)
    expected = (["X-PLOR V5.4321"], ["X-PLOR"], ["V5.4321"])
    eq_(result, expected)

    refprog = "X-PLOR V 5"
    result = parse_refprog(refprog)
    expected = (["X-PLOR V 5"], ["X-PLOR"], ["np"])
    eq_(result, expected)


def test_parse_refprog_phenix_ensemble_refinement():
    refprog = "PHENIX.ENSEMBLE_REFINEMENT"
    result = parse_refprog(refprog)
    expected = (["PHENIX.ENSEMBLE_REFINEMENT"], ["PHENIX.ENSEMBLE_REFINEMENT"],
                ["np"])
    eq_(result, expected)

    refprog = "PHENIX (PHENIX.ENSEMBLE_REFINEMENT: dev-9999)"
    result = parse_refprog(refprog)
    expected = (["PHENIX (PHENIX.ENSEMBLE_REFINEMENT: DEV-9999)"],
                ["PHENIX.ENSEMBLE_REFINEMENT"], ["DEV-9999"])
    eq_(result, expected)


def test_parse_refprog_phenix():
    refprog = "PHENIX"
    result = parse_refprog(refprog)
    expected = (["PHENIX"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX.REFINE"
    result = parse_refprog(refprog)
    expected = (["PHENIX.REFINE"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX.REFINEMENT"
    result = parse_refprog(refprog)
    expected = (["PHENIX.REFINEMENT"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX (PHENIX.REFINE)"
    result = parse_refprog(refprog)
    expected = (["PHENIX (PHENIX.REFINE)"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX (PHENIX)"
    result = parse_refprog(refprog)
    expected = (["PHENIX (PHENIX)"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX.REFINE (PHENIX.REFINE)"
    result = parse_refprog(refprog)
    expected = (["PHENIX.REFINE (PHENIX.REFINE)"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX.REFINEMENT (PHENIX.REFINE)"
    result = parse_refprog(refprog)
    expected = (["PHENIX.REFINEMENT (PHENIX.REFINE)"], ["PHENIX.REFINE"],
                ["-"])
    eq_(result, expected)

    refprog = "PHENIX.REFINEMENT (PHENIX)"
    result = parse_refprog(refprog)
    expected = (["PHENIX.REFINEMENT (PHENIX)"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX 1.8"
    result = parse_refprog(refprog)
    expected = (["PHENIX 1.8"], ["PHENIX.REFINE"], ["1.8"])
    eq_(result, expected)

    refprog = "PHENIX.REFINE 1.8AB"
    result = parse_refprog(refprog)
    expected = (["PHENIX.REFINE 1.8AB"], ["PHENIX.REFINE"], ["np"])
    eq_(result, expected)

    refprog = "PHENIX.REFINE: 1.8AB"
    result = parse_refprog(refprog)
    expected = (["PHENIX.REFINE: 1.8AB"], ["PHENIX.REFINE"], ["1.8AB"])
    eq_(result, expected)

    refprog = "PHENIX (1.8_1060)"
    result = parse_refprog(refprog)
    expected = (["PHENIX (1.8_1060)"], ["PHENIX.REFINE"], ["1.8_1060"])
    eq_(result, expected)

    refprog = "PHENIX (PHENIX.REFINE: 1.8)"
    result = parse_refprog(refprog)
    expected = (["PHENIX (PHENIX.REFINE: 1.8)"], ["PHENIX.REFINE"], ["1.8"])
    eq_(result, expected)

    refprog = "AND PHENIX.REFINE"
    result = parse_refprog(refprog)
    expected = (["PHENIX.REFINE"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)


def test_parse_refprog_buster():
    refprog = "BUSTER"
    result = parse_refprog(refprog)
    expected = (["BUSTER"], ["BUSTER"], ["-"])
    eq_(result, expected)

    refprog = "AUTOBUSTER 4"
    result = parse_refprog(refprog)
    expected = (["AUTOBUSTER 4"], ["BUSTER"], ["4"])
    eq_(result, expected)

    refprog = "BUSTER 2.10.0"
    result = parse_refprog(refprog)
    expected = (["BUSTER 2.10.0"], ["BUSTER"], ["2.10.0"])
    eq_(result, expected)

    refprog = "BUSTER-TNT 2.11.2"
    result = parse_refprog(refprog)
    expected = (["BUSTER-TNT 2.11.2"], ["BUSTER"], ["2.11.2"])
    eq_(result, expected)

    refprog = "BUSTER-TNT 2.X"
    result = parse_refprog(refprog)
    expected = (["BUSTER-TNT 2.X"], ["BUSTER"], ["2.X"])
    eq_(result, expected)

    refprog = "BUSTER-TNT V. 1.1.1"
    result = parse_refprog(refprog)
    expected = (["BUSTER-TNT V. 1.1.1"], ["BUSTER"], ["1.1.1"])
    eq_(result, expected)

    refprog = "BUSTER-TNT BUSTER 2.8.0"
    result = parse_refprog(refprog)
    expected = (["BUSTER-TNT BUSTER 2.8.0"], ["BUSTER"], ["2.8.0"])
    eq_(result, expected)

    refprog = "BUSTER TNT"
    result = parse_refprog(refprog)
    expected = (["BUSTER TNT"], ["BUSTER"], ["-"])
    eq_(result, expected)

    refprog = "BUSTER/TNT"
    result = parse_refprog(refprog)
    expected = (["BUSTER/TNT"], ["BUSTER"], ["-"])
    eq_(result, expected)

    refprog = "AUTOBUSTER"
    result = parse_refprog(refprog)
    expected = (["AUTOBUSTER"], ["BUSTER"], ["-"])
    eq_(result, expected)

    refprog = "BUSTER/TNT a b c 1 2 3"
    result = parse_refprog(refprog)
    expected = (["BUSTER/TNT A B C 1 2 3"], ["BUSTER"], ["np"])
    eq_(result, expected)



def test_parse_refprog_tnt():
    refprog = "TNT"
    result = parse_refprog(refprog)
    expected = (["TNT"], ["TNT"], ["-"])
    eq_(result, expected)

    refprog = "TNT V. 5-E"
    result = parse_refprog(refprog)
    expected = (["TNT V. 5-E"], ["TNT"], ["5-E"])
    eq_(result, expected)

    refprog = "TNT 5.6.1"
    result = parse_refprog(refprog)
    expected = (["TNT 5.6.1"], ["TNT"], ["5.6.1"])
    eq_(result, expected)

    refprog = "TNT 5F"
    result = parse_refprog(refprog)
    expected = (["TNT 5F"], ["TNT"], ["5F"])
    eq_(result, expected)

    refprog = "TNT V. 5-F PRERELEASE"
    result = parse_refprog(refprog)
    expected = (["TNT V. 5-F PRERELEASE"], ["TNT"], ["5-F PRERELEASE"])
    eq_(result, expected)

    refprog = "TNT 5-F PRERELEASE"
    result = parse_refprog(refprog)
    expected = (["TNT 5-F PRERELEASE"], ["TNT"], ["5-F PRERELEASE"])
    eq_(result, expected)

    refprog = "TNT 5F 5G"
    result = parse_refprog(refprog)
    expected = (["TNT 5F 5G"], ["TNT"], ["np"])
    eq_(result, expected)


def test_parse_refprog_shelx():
    refprog = "SHELX"
    result = parse_refprog(refprog)
    expected = (["SHELX"], ["SHELX"], ["-"])
    eq_(result, expected)

    refprog = "SHELX97"
    result = parse_refprog(refprog)
    expected = (["SHELX97"], ["SHELX"], ["97"])
    eq_(result, expected)

    refprog = "SHELXH"
    result = parse_refprog(refprog)
    expected = (["SHELXH"], ["SHELX"], ["H"])
    eq_(result, expected)

    refprog = "SHELXS"
    result = parse_refprog(refprog)
    expected = (["SHELXS"], ["SHELX"], ["S"])
    eq_(result, expected)

    refprog = "SHELXL"
    result = parse_refprog(refprog)
    expected = (["SHELXL"], ["SHELX"], ["L"])
    eq_(result, expected)

    refprog = "SHELX-97"
    result = parse_refprog(refprog)
    expected = (["SHELX-97"], ["SHELX"], ["97"])
    eq_(result, expected)

    refprog = "SHELXL-97"
    result = parse_refprog(refprog)
    expected = (["SHELXL-97"], ["SHELX"], ["L-97"])
    eq_(result, expected)

    refprog = "SHELXH-97"
    result = parse_refprog(refprog)
    expected = (["SHELXH-97"], ["SHELX"], ["H-97"])
    eq_(result, expected)

    refprog = "SHELX-97-1"
    result = parse_refprog(refprog)
    expected = (["SHELX-97-1"], ["SHELX"], ["97-1"])
    eq_(result, expected)

    refprog = "SHELX 97-1"
    result = parse_refprog(refprog)
    expected = (["SHELX 97-1"], ["SHELX"], ["97-1"])
    eq_(result, expected)

    refprog = "SHELXL 97"
    result = parse_refprog(refprog)
    expected = (["SHELXL 97"], ["SHELX"], ["L 97"])
    eq_(result, expected)

    refprog = "SHELXL97"
    result = parse_refprog(refprog)
    expected = (["SHELXL97"], ["SHELX"], ["L97"])
    eq_(result, expected)

    refprog = "SHELX-L"
    result = parse_refprog(refprog)
    expected = (["SHELX-L"], ["SHELX"], ["L"])
    eq_(result, expected)

    refprog = "SHELXL-L"
    result = parse_refprog(refprog)
    expected = (["SHELXL-L"], ["SHELX"], ["np"])
    eq_(result, expected)


def test_parse_refprog_arpwarp():
    refprog = "ARP/WARP"
    result = parse_refprog(refprog)
    expected = (["ARP/WARP"], ["ARP/WARP"], ["-"])
    eq_(result, expected)

    refprog = "ARP/WARP V. 1.0"
    result = parse_refprog(refprog)
    expected = (["ARP/WARP V. 1.0"], ["ARP/WARP"], ["1.0"])
    eq_(result, expected)

    refprog = "ARP/WARP 1.0"
    result = parse_refprog(refprog)
    expected = (["ARP/WARP 1.0"], ["ARP/WARP"], ["np"])
    eq_(result, expected)


def test_parse_refprog_coot():
    refprog = "COOT"
    result = parse_refprog(refprog)
    expected = (["COOT"], ["COOT"], ["-"])
    eq_(result, expected)

    refprog = "COOT 0.1-3-PRE-1"
    result = parse_refprog(refprog)
    expected = (["COOT 0.1-3-PRE-1"], ["COOT"], ["0.1-3-PRE-1"])
    eq_(result, expected)

    refprog = "COOT 0.0.33"
    result = parse_refprog(refprog)
    expected = (["COOT 0.0.33"], ["COOT"], ["0.0.33"])
    eq_(result, expected)

    refprog = "COOT V. 0.0.33"
    result = parse_refprog(refprog)
    expected = (["COOT V. 0.0.33"], ["COOT"], ["0.0.33"])
    eq_(result, expected)

    refprog = "COOT $"
    result = parse_refprog(refprog)
    expected = (["COOT $"], ["COOT"], ["np"])
    eq_(result, expected)


def test_parse_refprog_o():
    refprog = "O"
    result = parse_refprog(refprog)
    expected = (["O"], ["O"], ["-"])
    eq_(result, expected)

    refprog = "O 9.0.7"
    result = parse_refprog(refprog)
    expected = (["O 9.0.7"], ["O"], ["np"])
    eq_(result, expected)

    refprog = "O9.0.7"
    result = parse_refprog(refprog)
    expected = (["O9.0.7"], ["O"], ["9.0.7"])
    eq_(result, expected)


def test_parse_refprog_profft():
    refprog = "PROFFT"
    result = parse_refprog(refprog)
    expected = (["PROFFT"], ["PROLSQ"], ["PROFFT"])
    eq_(result, expected)

    refprog = "PROFFT (MODIFIED BY Z.OTWINOWSKI)"
    result = parse_refprog(refprog)
    expected = ([refprog], ["PROLSQ"], ["PROFFT MODIFIED BY Z.OTWINOWSKI"])
    eq_(result, expected)

    refprog = "PROFFT anything else"
    result = parse_refprog(refprog)
    expected = (["PROFFT ANYTHING ELSE"], ["PROLSQ"], ["np"])
    eq_(result, expected)


def test_parse_refprog_prolsq():
    refprog = "PROLSQ"
    result = parse_refprog(refprog)
    expected = (["PROLSQ"], ["PROLSQ"], ["PROLSQ"])
    eq_(result, expected)

    refprog = "PROLSQ (MODIFIED BY G.J.QUIGLEY)"
    result = parse_refprog(refprog)
    expected = ([refprog], ["PROLSQ"], ["PROLSQ MODIFIED BY G.J.QUIGLEY"])
    eq_(result, expected)

    refprog = "PROLSQ anything else"
    result = parse_refprog(refprog)
    expected = (["PROLSQ ANYTHING ELSE"], ["PROLSQ"], ["np"])
    eq_(result, expected)


def test_parse_refprog_protin():
    refprog = "PROTIN"
    result = parse_refprog(refprog)
    expected = (["PROTIN"], ["PROLSQ"], ["PROLSQ"])
    eq_(result, expected)

    refprog = "PROTIN anything else"
    result = parse_refprog(refprog)
    expected = (["PROTIN ANYTHING ELSE"], ["PROLSQ"], ["np"])
    eq_(result, expected)


def test_parse_refprog_nuclsq():
    refprog = "NUCLSQ"
    result = parse_refprog(refprog)
    expected = (["NUCLSQ"], ["PROLSQ"], ["NUCLSQ"])
    eq_(result, expected)

    refprog = "NUCLSQ (MODIFIED BY G.J.QUIGLEY)"
    result = parse_refprog(refprog)
    expected = ([refprog], ["PROLSQ"], ["NUCLSQ MODIFIED BY G.J.QUIGLEY"])
    eq_(result, expected)

    refprog = "NUCLSQ anything else"
    result = parse_refprog(refprog)
    expected = (["NUCLSQ ANYTHING ELSE"], ["PROLSQ"], ["np"])
    eq_(result, expected)


def test_parse_refprog_nuclin():
    refprog = "NUCLIN"
    result = parse_refprog(refprog)
    expected = (["NUCLIN"], ["PROLSQ"], ["NUCLSQ"])
    eq_(result, expected)

    refprog = "PROTIN anything else"
    result = parse_refprog(refprog)
    expected = (["PROTIN ANYTHING ELSE"], ["PROLSQ"], ["np"])
    eq_(result, expected)


def test_parse_refprog_gprlsa():
    refprog = "GPRLSA"
    result = parse_refprog(refprog)
    expected = (["GPRLSA"], ["PROLSQ"], ["GPRLSA"])
    eq_(result, expected)

    refprog = "GPRLSA anything else"
    result = parse_refprog(refprog)
    expected = (["GPRLSA ANYTHING ELSE"], ["PROLSQ"], ["np"])
    eq_(result, expected)


def test_parse_refprog_deriv():
    refprog = "DERIV"
    result = parse_refprog(refprog)
    expected = (["DERIV"], ["PROLSQ"], ["DERIV"])
    eq_(result, expected)

    refprog = "DERIV anything else"
    result = parse_refprog(refprog)
    expected = (["DERIV ANYTHING ELSE"], ["PROLSQ"], ["np"])
    eq_(result, expected)


def test_parse_refprog_cerius():
    refprog = "CERIUS"
    result = parse_refprog(refprog)
    expected = (["CERIUS"], ["CERIUS"], ["-"])
    eq_(result, expected)

    refprog = "CERIUS 2"
    result = parse_refprog(refprog)
    expected = (["CERIUS 2"], ["CERIUS"], ["2"])
    eq_(result, expected)

    refprog = "CERIUS2"
    result = parse_refprog(refprog)
    expected = (["CERIUS2"], ["CERIUS"], ["2"])
    eq_(result, expected)

    refprog = "CERIUS anything else"
    result = parse_refprog(refprog)
    expected = (["CERIUS ANYTHING ELSE"], ["CERIUS"], ["np"])
    eq_(result, expected)


def test_parse_refprog_hkl3000():
    refprog = "HKL-3000"
    result = parse_refprog(refprog)
    expected = (["HKL-3000"], ["HKL-3000"], ["-"])
    eq_(result, expected)

    refprog = "HKL3000"
    result = parse_refprog(refprog)
    expected = (["HKL3000"], ["HKL-3000"], ["-"])
    eq_(result, expected)

    refprog = "HKL-3000 anything else"
    result = parse_refprog(refprog)
    expected = (["HKL-3000 ANYTHING ELSE"], ["HKL-3000"], ["np"])
    eq_(result, expected)


def test_parse_refprog_gromos():
    refprog = "GROMOS"
    result = parse_refprog(refprog)
    expected = (["GROMOS"], ["GROMOS"], ["-"])
    eq_(result, expected)

    refprog = "GROMOS87"
    result = parse_refprog(refprog)
    expected = (["GROMOS87"], ["GROMOS"], ["87"])
    eq_(result, expected)

    refprog = "GROMOS anything else"
    result = parse_refprog(refprog)
    expected = (["GROMOS ANYTHING ELSE"], ["GROMOS"], ["np"])
    eq_(result, expected)


def test_parse_refprog_other():
    refprogs = ["AGARWAL FAST-FOURIER TRANSFORM LEAST-SQUARES",
                "HENDRICKSON-KONNERT LEAST-SQUARES REFINEMENT",
                "* OF T. A. JONES.  THE R VALUE IS 0.180.",
                "RESTRAINED RECIPROCAL-SPACE LEAST-SQUARES",
                "SIMULATED ANNEALING METHOD",
                "CONSTRAINED RECIPROCAL-SPACE LEAST-SQUARES",
                "FAST-FOURIER LEAST-SQUARES REFINEMENT",
                "JACK-LEVITT",
                "REAL-SPACE REFINEMENT",
                "SOME OTHER REFPROG",]
    for p in refprogs:
        result = parse_refprog(p)
        expected = ([p], ["OTHER"], ["np"])
