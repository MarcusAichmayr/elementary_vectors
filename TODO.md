# TODO

## elementary_vectors

* [x] Add example with `return_minors=True` in `elementary_vectors` and `elementary_vectors_from_matrix`.
* [x] Add examples to `non_negative_vectors`.
* [x] Add examples to `positive_elementary_vectors`.
    * [ ] add examples with variables

### elementary_vectors.elementary_vectors

* Do we want more optional arguments for `elementary_vectors` to have more control on reductions?
    * optional argument: `cancel_common_factor`?
* [x] add examples for other rings: QQ, GF(7)

* [x] Change output if `return_minors=True` to `[evs, m]` instead of `[m, evs]`.


### elementary_vectors.positive_elementary_vectors
* In `positive_elementary_vector`, all minors could be 0.
  In this case, the result is wrong.
  We should recompute the elementary vectors in this case.
* Would it make sense to pass the elementary vectors to `positive_elementary_vector`?
  A problem would be that for a certain choice of parameters all elementary vectors can be 0.
* Wird vor `positive_elementary_vector` zum Beispiel `x > 0` angenommen,
  so kommt diese Bedingung nicht beim Output vor.
  Man könnte natürlich `assumptions()` dazugeben, anstatt eine Liste rekursiv zu befüllen.
  Jedoch werden dann alle assumptions ausgegeben. (Also auch `a is integer`.)
* Durch `reduce=True` könnte der Elementary Vector `(1,1,0,0)` aus `(x,x,0,0)` entstehen.
  Ist aber durch Annahme `x=0`, so haben wir den Vektor `(1,1,0,0)` der kein Elementary Vector ist.

CoCalc:
* Dateien zu Elementary Vectors zusammenfügen.


### elementary_vectors.utility

* `sign_determined`
    * [x] improve name (was `has_sign`)
    * [x] add docstring
    * [x] add tests

* `conformal_elimination`
    * Should this be here?
    * [ ] add examples
    * [ ] add tests

### elementary_vectors.reductions

* [x] improve docstrings
    * [x] add examples

### elementary_vectors.exists_vector

* [ ] If True and certificate is True, compute a corresponding vector of the vector space. (witness)
* [x] add examples


## sign_vectors

### sign_vectors.sign_vectors

* [x] add docstrings
* [x] use doctests instead of separate test file

### sign_vectors.utility

* [x] add examples
* [x] move closure, contraction, deletion into this file
* [x] use doctests instead of separate test file

#### closure

* [x] improve docstring
    * [x] add examples
    * [x] define closure
* [x] add tests
* [x] use doctests instead of separate test file

#### contraction, deletion

* [x] add examples
* [x] add tests
* [x] use doctests instead of separate test file

#### subvector

* [x] add examples

### sign_vectors.oriented_matroids

* [x] add examples
    * (optional) add examples for all functions
* [x] use only one file (called `oriented_matroids.py`) for all functions:
    * `cocircuits_from_matrix`
    * `covectors_from_cocircuits`
    * `topes_from_cocircuits`
    * `lower_faces`
    * `face_enumeration`
    
    * `topes_from_matrix`
    * `covectors_from_topes`
    * `cocircuits_from_topes`
    * `covectors_from_matrix`
    * [ ] test files in CoCalc.
      `from_sign_vectors.oriented_matroids import *` should still work.

## Sphinx documentation

I adapted the files from here: https://github.com/mkauers/ore_algebra

* [x] Make this work.

