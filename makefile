.PHONY: install test

install:
	sage -pip install --upgrade .

test:
	sage -t elementary_vectors/ sign_vectors/ vectors_in_intervals/

doc:
	cd docs && make html

doc-pdf:
	cd docs && make latexpdf
