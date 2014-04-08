from nose.tools import eq_

from bdb.expdta import check_exp_methods


def test_check_exp_methods_too_many():
    records = {"EXPDTA": ["    NEUTRON DIFFRACTION; X-RAY DIFFRACTION",]}
    result = check_exp_methods(records, "1crn")
    eq_(result["expdta_useful"], False)
    eq_(result["expdta"], ["NEUTRON DIFFRACTION", "X-RAY DIFFRACTION"])


def test_check_exp_methods_none():
    records = {"REMARK": ["    TEST",]}
    result = check_exp_methods(records, "1crn")
    eq_(result["expdta_useful"], False)
    eq_(result["expdta"], [])


def test_check_exp_methods_not_suitable():
    records = {"EXPDTA": ["    NEUTRON DIFFRACTION",]}
    result = check_exp_methods(records, "1crn")
    eq_(result["expdta_useful"], False)
    eq_(result["expdta"], ["NEUTRON DIFFRACTION"])


def test_check_exp_methods():
    records = {"EXPDTA": ["    X-RAY DIFFRACTION",]}
    result = check_exp_methods(records, "1crn")
    eq_(result["expdta_useful"], True)
    eq_(result["expdta"], ["X-RAY DIFFRACTION"])
