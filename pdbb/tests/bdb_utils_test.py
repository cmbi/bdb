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
from nose.tools import eq_, raises

from pdbb.bdb_utils import (is_valid_directory, is_valid_file, is_valid_pdbid)

import argparse


def test_is_valid_directory():
    """Trivial test to check if directory exists."""
    parser = argparse.ArgumentParser()
    direct = "pdbb/tests/pdb/files"
    eq_(is_valid_directory(parser, direct), direct)


@raises(SystemExit)
def test_is_valid_directory_error():
    """Trivial test to check if directory exists."""
    parser = argparse.ArgumentParser()
    direct = "foo"
    is_valid_directory(parser, direct)


def test_is_valid_file():
    """Trivial test to check if file exists."""
    parser = argparse.ArgumentParser()
    file_path = "pdbb/tests/pdb/files/1crn.pdb"
    eq_(is_valid_file(parser, file_path), file_path)


@raises(SystemExit)
def test_is_valid_file_error():
    """Trivial test to check if file exists."""
    parser = argparse.ArgumentParser()
    file_path = "foo"
    is_valid_file(parser, file_path)


@raises(SystemExit)
def test_is_valid_file_empty():
    """Trivial test to check if file is empty."""
    parser = argparse.ArgumentParser()
    file_path = "pdbb/tests/pdb/files/empty"
    is_valid_file(parser, file_path)


def test_is_valid_pdbid():
    """Trivial test to check if the pdbid is valid."""
    parser = argparse.ArgumentParser()
    pdb_id = "1crn"
    eq_(is_valid_pdbid(parser, pdb_id), pdb_id)

    pdb_id = "1CRN"
    eq_(is_valid_pdbid(parser, pdb_id), pdb_id)

    pdb_id = "CcRn"
    eq_(is_valid_pdbid(parser, pdb_id), pdb_id)

    pdb_id = "1111"
    eq_(is_valid_pdbid(parser, pdb_id), pdb_id)


@raises(TypeError)
def test_is_valid_pdbid_error():
    """Trivial test to check if the pdbid is valid."""
    parser = argparse.ArgumentParser()

    pdb_id = None
    is_valid_pdbid(parser, pdb_id)

    pdb_id = 1
    is_valid_pdbid(parser, pdb_id)


@raises(SystemExit)
def test_is_valid_pdbid_false():
    """Trivial test to check if the pdbid is valid."""
    parser = argparse.ArgumentParser()
    pdb_id = "11111"
    is_valid_pdbid(parser, pdb_id)

    pdb_id = "111"
    is_valid_pdbid(parser, pdb_id)

    pdb_id = "(crn"
    is_valid_pdbid(parser, pdb_id)
