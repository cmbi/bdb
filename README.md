# bdb
If protein engineers, homology modellers, biologists, or bioinformaticians
need B-factors from a PDB file, they normally want full isotropic B-factors.
Normally it is indeed the full B-factor that is stored in the B-factor field
in the ATOM records of a PDB file. However, sometimes the field contains
"residual" B-factors or atomic mean-square displacements instead of B-factors.
The BDB contains PDB files with full isotropic B-factors in the B-factor field
if the original PDB file contains enough information to determine the content
of the B-factor field and calculate the full B-factor if necessary.

BDB files can be viewed and downloaded [here][1].
More information on BDB entries, their construction, their web page, how to
download the entire BDB databank, and how to cite the BDB can be found
[here][2].
</br>

# License

GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

# Installation

## From source

Download the source code and run the python installer:

* `wget https://github.com/cmbi/bdb/releases/bdb-<version>.tar.gz`
* `tar -zxvf bdb-<version>.tar.gz`
* `cd bdb-<version>`
* `sudo python setup.py install`

# Development

If you'd like to contribute by adding features or fixing bugs, follow the steps
in this section to setup a development environment.

## Pre-requisites

The following pre-requisites are required by bdb and must be installed
manually:

* virtualenv
* virtualenvwrapper
* [ccp4 software suite][3]

## Setup

* Run `mkvirtualenv --no-site-packages bdb` to create a virtual environment,
  and switch to that environment with `workon bdb`.
* Run `git clone https://github.com/cmbi/bdb.git` to obtain the latest
  source code.
* Run `git checkout develop` to switch to the development branch. It's also
  recommended that you create a feature branch off develop.
* Run `pip install -r requirements` to install the module dependencies.

The last command installs the required module dependencies into the virtual
environment you've created, isolating the development environment from the rest
of your system.

## Tests

The following command will run all the unit tests for bdb:

    nosetests

To see the code coverage of the unit tests, run:

    nosetests --with-coverage --cover-package=pdbb

[1]: http://www.cmbi.umcn.nl/bdb/
[2]: http://www.cmbi.umcn.nl/bdb/about/
[3]: http://www.ccp4.ac.uk/
