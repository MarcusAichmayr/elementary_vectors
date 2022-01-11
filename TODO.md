# TODO

## elementary_vectors

* [ ] Add explanations and examples to demonstrate functions in this module. (in `__init__`?)
  * several reductions
  * `elementary_vectors`
  * `exists_vector`
* Every .py file should have a docstring.
  * [ ] `vectors_in_intervals`
    * [ ] check for existence of vectors, compute vectors

### elementary_vectors.functions

* [x] Add example with `return_minors=True` in `elementary_vectors` and `elementary_vectors_from_matrix`.
* [x] Add examples to `non_negative_vectors`.
* [x] Add examples to `positive_elementary_vectors`.
  * [ ] add examples with variables
* Do we want more optional arguments for `elementary_vectors` to have more control on reductions?
  * optional argument: `cancel_common_factor`?
* [x] add examples for other rings: `QQ`, `GF(7)`
* [x] Change output if `return_minors=True` to `[evs, m]` instead of `[m, evs]`.


#### elementary_vectors.functions.positive_elementary_vectors

* In `positive_elementary_vector`, all minors could be 0.
  In this case, the result is wrong.
  We should recompute the elementary vectors in this case.
* Would it make sense to pass the elementary vectors to `positive_elementary_vector`?
  A problem would be that for a certain choice of parameters all elementary vectors can be 0.
* If we assume for instance `x > 0` before we call `positive_elementary_vector`,
  then this condition is not in the output.
  To prevent this, we could add `assumptions()`, instead of filling a list recursively.
  But then, we would return all assumptions. (e.g. `a is integer`)
* With `reduce=True`, we could obtain the elementary vector `(1,1,0,0)` from `(x,x,0,0)`.
  For `x=0`, we have the wrong elementary vector `(1,1,0,0)`.

CoCalc:
* Put functions for elementary vectors in one file. For maintaining, it would be easier to work more with the documentation.

### elementary_vectors.vectors_in_intervals

* [x] implement `construct_vector`
  * [ ] add examples
    * [ ] several examples, where it works
    * [ ] easy examples
    * [x] examples where it does not work

* [x] implement `construct_normal_vector`
  * [x] Fix bug:
    If there exists no solution, `construct_normal_vector` might still find something.
  * [ ] add examples
    * [ ] several examples, where it works
    * [ ] easy examples
    * [ ] examples where it does not work


#### elementary_vectors.vectors_in_intervals.exists_vector

* [ ] If True and certificate is True, compute a corresponding vector of the vector space. (witness)
* [x] add examples
* [ ] Should we remove certificate part?
  * If `construct_vector` is implemented, do we need this?
* [ ] Do we need the argument `kernel`?


### elementary_vectors.utility

* `sign_determined`
  * [x] improve name (was `has_sign`)
  * [x] add docstring
  * [x] add tests

* `conformal_elimination`
  * [ ] Should this be here?
  * [ ] add examples
  * [ ] add tests

### elementary_vectors.reductions

* [x] improve docstrings
  * [x] add examples
* [ ] add examples to docstring on top


## sign_vectors

### sign_vectors.sign_vectors

* [x] Add explanations and examples to demonstrate functions in this module on top of this file.
* [x] add docstrings
* [x] use doctests instead of separate test file

### sign_vectors.oriented_matroids

* [x] add examples
  * [x] add examples to all important functions
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
* [ ] Is it possible to move references at the end of the documentation?

### sign_vectors.utility

* [x] add examples
* [x] move closure, contraction, deletion into this file
* [x] use doctests instead of separate test file
* [ ] Can we improve the implementation of `closure`?

## Sphinx documentation

I adapted the files from here: https://github.com/mkauers/ore_algebra

* [x] Make this work.
* Do we want to have tests in the documentation?
