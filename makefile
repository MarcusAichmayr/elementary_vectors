.PHONY: install test

install:
	sage -pip install --upgrade --no-index .

test_elementary_vectors:
	sage tests/elementary_vectors/test_elementary_vectors.py -v
	sage tests/elementary_vectors/test_utility.py -v
	sage tests/elementary_vectors/test_reductions.py -v
	sage tests/elementary_vectors/test_exists_vector.py -v

test_sign_vectors:
	sage tests/sign_vectors/test_sign_vectors.py -v # run unittests for class Signvector
	sage tests/sign_vectors/test_utility.py -v # run unittests for utility functions
	sage tests/sign_vectors/oriented_matroids/test_functions.py -v # run test if appropriate functions are available # Todo: Should be done for each module separately.

test:
	make test_elementary_vectors
	make test_sign_vectors

