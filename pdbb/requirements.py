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
import logging
_log = logging.getLogger(__name__)

import subprocess
import sys


# CCP4 dependencies
ccp4_software = [
    "tlsanl", ]

# Other dependencies
# None

"""
Test if CCP4 environment has been set up properly.

Strategy:
    Since CCP4 programs normally only except interactive input,
    we cannot test for a 0 exit code without having to use example files.
    Instead, we catch a subprocess.CalledProcessError (raised in case the
    program exists) before an OSError (raised if the program cannot be called),
    after our attempt to call the program with an invalid argument.
"""
for p in ccp4_software:
    try:
        p = subprocess.check_call(
            [p, "and_an_invalid_argument"],
            stdin=subprocess.PIPE,   # suppress output
            stdout=subprocess.PIPE,  # suppress output
            stderr=subprocess.PIPE)  # suppress output
    except subprocess.CalledProcessError:
        _log.debug("{0:s} set up properly".format(p))
    except OSError:
        _log.error("{0:s} could not be executed. Install {0:s} and set up "
                   "the CCP4 environment properly.".format(p))
        sys.exit(1)
