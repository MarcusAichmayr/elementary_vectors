.PHONY: install test

install:
	sage -pip install --upgrade --no-index .

test_elementary_vectors:
	sage -t elementary_vectors/elementary_vectors.py
	sage -t elementary_vectors/utility.py
	sage tests/elementary_vectors/test_reductions.py -v
	sage tests/elementary_vectors/test_exists_vector.py -v

test_sign_vectors:
	sage -t sign_vectors/sign_vectors.py
	sage tests/sign_vectors/test_utility.py -v # run unittests for utility functions
	sage tests/sign_vectors/oriented_matroids/test_functions.py -v # run test if appropriate functions are available # Todo: Should be done for each module separately.

test:
	make test_elementary_vectors
	make test_sign_vectors

