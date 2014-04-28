from nose.tools import eq_

from bdb.refprog import (is_bdb_includable_refprog)


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

