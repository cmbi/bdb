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
