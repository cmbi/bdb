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
from nose.tools import eq_

from pdbb.expdta import check_exp_methods


def test_check_exp_methods_too_many():
    records = {"EXPDTA": ["    NEUTRON DIFFRACTION; X-RAY DIFFRACTION", ]}
    result = check_exp_methods(records, "test")
    eq_(result["expdta_useful"], False)
    eq_(result["expdta"], ["NEUTRON DIFFRACTION", "X-RAY DIFFRACTION"])


def test_check_exp_methods_none():
    records = {"REMARK": ["    TEST", ]}
    result = check_exp_methods(records, "test")
    eq_(result["expdta_useful"], False)
    eq_(result["expdta"], [])


def test_check_exp_methods_not_suitable():
    records = {"EXPDTA": ["    NEUTRON DIFFRACTION                        ", ]}
    result = check_exp_methods(records, "test")
    eq_(result["expdta_useful"], False)
    eq_(result["expdta"], ["NEUTRON DIFFRACTION"])


def test_check_exp_methods():
    records = {"EXPDTA": ["    X-RAY DIFFRACTION", ]}
    result = check_exp_methods(records, "test")
    eq_(result["expdta_useful"], True)
    eq_(result["expdta"], ["X-RAY DIFFRACTION"])
