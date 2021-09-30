.PHONY: install test

install:
	sage -pip install --upgrade --no-index .

test_elementary_vectors:
	sage -t elementary_vectors/

test_sign_vectors:
	sage -t sign_vectors/

test:
	sage -t elementary_vectors/ sign_vectors/

doc:
	cd docs && make html

doc-pdf:
	cd docs && make latexpdf
