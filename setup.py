from distutils.core import setup


setup(
    name='bdb',
    version='0.5.3.5',
    description='A databank of PDB entries with full isotropic B-factors.',
    author='Wouter Touw',
    author_email='wouter.touw@radboudumc.nl',
    url='https://github.com/cmbi/bdb',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=[
        'bdb',
        'bdb.pdb',
        'bdb.tests',
        'bdb.tests.pdb',
    ],
    scripts=['scripts/bdb', ],
)
