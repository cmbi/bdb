# bdb

**TODO: Write a short paragraph introducing bdb**

# License

**TODO: Pick a license and write down here which one it is**

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
* [ccp4 software suite][1]

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

    nosetests --with-coverage --cover-package=bdb

[1]: http://www.ccp4.ac.uk/
