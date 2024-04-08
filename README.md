# Sign vectors

## Description

a Sage package to work with elementary vectors, sign vectors and oriented matroids

This package is still in development. Please inform me if something doesn't work as expected.

## License

Distributed under the terms of the GNU General Public License (GPL, see the
LICENSE file), either version 3 or (at your option) any later version

- http://www.gnu.org/licenses/

## Requirements

Sage 9.0 or later is recommended. Some features should work with older versions.

## Installation

### Install from GitHub (recommended)

To download and install the latest development version on a system where Sage
was built from source or installed from official packages, run

    sage -pip install git+https://github.com/MarcusAichmayr/elementary_vectors.git

or

    sage -pip install --user git+https://github.com/MarcusAichmayr/elementary_vectors.git

The optional `--user` flag causes the package to be installed in your `.sage` directory instead of the Sage installation tree.

### Local install from source

Download the source from the git repository:

    git clone https://github.com/MarcusAichmayr/elementary_vectors.git

Change to the root directory of the repository and run:

    sage -pip install --upgrade --no-index -v .

You can also run instead the shorthand:

    make install

### Documentation

The documentation of this package can be found on GitHub:
https://marcusaichmayr.github.io/elementary_vectors/index.html

If you want to generate it yourself, run

    make doc

or

    make doc-pdf

at the root directory of the repository.

### Testing

To run the test suite, install the package and run the command

    make test

at the root directory of the repository.
