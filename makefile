.PHONY: install test

install:
	sage -pip install --upgrade --no-index .

test:
    # Todo test elementary_vectors
	sage tests/test_sign_vectors.py -v # run unittests for class Signvector
	sage tests/test_utility.py -v # run unittests for utility functions
	sage tests/test_functions/test_functions.py -v # run test if appropriate functions are available # Todo: Should be done for each module separately.
	sage tests/test_exists_vector.py -v # run tests

