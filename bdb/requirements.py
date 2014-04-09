#!/usr/bin/env python
import logging
import subprocess
import sys

# Configure logging
_log = logging.getLogger("bdb")

# CCP4 dependencies
ccp4_software = [
    "tlsanl",
    ]

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
