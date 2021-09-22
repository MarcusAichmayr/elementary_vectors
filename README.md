# Sign vectors

## Description

a Sage package to work with sign vectors, oriented matroids and elementary vectors

## License

Distributed under the terms of the GNU General Public License (GPL, see the
LICENSE file), either version 3 or (at your option) any later version

- http://www.gnu.org/licenses/

## Requirements

Sage 9.0 or later is recommended. Some features should work with older versions.

## Installation

### Local install from source

Download the source from the git repository::

    $ git clone https://github.com/MarcusAichmayr/sign_vectors.git

Change to the root directory and run::

    $ sage -pip install --upgrade --no-index -v .

For convenience this package contains a `makefile <makefile>`_ with this
and other often used commands. Should you wish too, you can use the
shorthand::

    $ make install

### Install from GitHub

To download and install the latest development version on a system where Sage
was built from source or installed from official packages, run

    $ sage -pip install git+https://github.com/MarcusAichmayr/sign_vectors.git

or

    $ sage -pip install --user git+https://github.com/MarcusAichmayr/sign_vectors.git

The optional --user flag causes the package to be installed in your .sage directory instead of the Sage installation tree.

Testing
-------

To run the test suite, install the package and run the command

    $ make test

at the root of the git checkout.
