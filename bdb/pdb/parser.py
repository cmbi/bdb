import logging
import os


_log = logging.getLogger(__name__)


def parse_pdb_file(pdb_file_path):
    """
    Parses the given pdb file, returning a dict where the key is the
    record name (e.g. 'ATOM   ') and the value is a list of all lines of that
    record name type.

    No validation is performed on the content of the pdb file.

    If the file at pdb_file_path doesn't exist, a ValueError is raised.
    """
    _log.info("Parsing pdb file {}".format(pdb_file_path))

    if not os.path.exists(pdb_file_path):
        _log.error("'{}' not found".format(pdb_file_path))
        raise ValueError("'{}' not found".format(pdb_file_path))

    with open(pdb_file_path) as pdb_file:
        records = {}
        for record in pdb_file:
            record_name = record[0:6]

            # If this is the first occurrence of a record name, initialise
            # the value with an empty list.
            if record_name not in records:
                records[record_name] = []
            records[record_name].append(record[7:])
        _log.debug("Parsed {0} records".format(len(records)))
        return records


def parse_exp_methods(pdb_records):
    """
    Parses the experiment methods from the pdb EXPDTA records, returning a
    list.

    The EXPDTA record is mandatory and can span multiple lines. If more than
    one experiment method is used, they are separated by a semi-colon.

    Although only a limited number of experiment methods are allowed, this
    function performs no validation of the values provided.

    If no EXPDTA records are found, a ValueError is raised.
    If no experimental method is found, an empty list is returned.
    """
    _log.info("Parsing experiment methods from EXPDTA records")
    if "EXPDTA" not in pdb_records:
        _log.error("No EXPDTA records found")
        raise ValueError("No EXPDTA records found")

    exp_methods = []
    for record in pdb_records["EXPDTA"]:
        for method in record.split(";"):
            exp_methods.append(method.strip())
    _log.debug("Found {} experiment methods".format(len(exp_methods)))
    return exp_methods
