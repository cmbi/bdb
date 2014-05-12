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
from distutils.core import setup


setup(
    name='bdb',
    version='0.6.1',
    description='A databank of PDB entries with full isotropic B-factors.',
    author='Wouter Touw',
    author_email='wouter.touw@radboudumc.nl',
    license = "GNU GPLv3",
    url='https://github.com/cmbi/bdb',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=[
        'pdbb',
        'pdbb.pdb',
        'pdbb.tests',
        'pdbb.tests.pdb',
    ],
    scripts=['scripts/mkbdb', ],
)
