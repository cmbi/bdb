from nose.tools import eq_, raises

from bdb.bdb_utils import (is_valid_directory, is_valid_file, is_valid_pdbid)

import argparse


def test_is_valid_directory():
    """Trivial test to check if directory exists."""
    parser = argparse.ArgumentParser()
    direct = "bdb/tests/pdb/files"
    eq_(is_valid_directory(parser, direct), direct)


def test_is_valid_file():
    """Trivial test to check if file exists."""
    parser = argparse.ArgumentParser()
    file_path = "bdb/tests/pdb/files/1crn.pdb"
    eq_(is_valid_file(parser, file_path), file_path)


@raises(SystemExit)
def test_is_valid_file_error():
    """Trivial test to check if file exists."""
    parser = argparse.ArgumentParser()
    file_path = "1crn.pdb"
    is_valid_file(parser, file_path)


@raises(SystemExit)
def test_is_valid_file_error():
    """Trivial test to check if file is empty."""
    parser = argparse.ArgumentParser()
    file_path = "bdb/tests/pdb/files/empty"
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
def test_is_valid_pdbid():
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


