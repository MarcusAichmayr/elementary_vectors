# TODO

## 1. elementary_vectors

### 1.1. elementary_vectors.functions

#### 1.1.1. elementary_vectors.functions.positive_elementary_vectors

* In `positive_elementary_vector`, all minors could be 0.
  In this case, the result is wrong.
  We should recompute the elementary vectors in this case.
* Would it make sense to pass the elementary vectors to `positive_elementary_vector`?
  A problem would be that for a certain choice of parameters all elementary vectors can be 0.
* If we assume for instance `x > 0` before we call `positive_elementary_vector`,
  then this condition is not in the output.
  To prevent this, we could add `assumptions()`, instead of filling a list recursively.
  But then, we would return all assumptions. (e.g. `a is integer`)
* With `cancel_factor=True`, we could obtain the elementary vector `(1,1,0,0)` from `(x,x,0,0)`.
  For `x=0`, we have a wrong elementary vector `(1,1,0,0)`.

### 1.2. elementary_vectors.functions_dd

* `determine_sign(X, a, M)`
  * Should this function work for real vectors (as `X` and `a`)? Does not make much sense.
  * [ ] This function works for cocircuits. For other sign vectors, it might fail.

### 1.4. elementary_vectors.utility

* `kernel_vector_support_given`
  * [ ] rename this function

* `conformal_elimination`
  * [ ] should this be here?
  * [ ] add examples
  * [ ] add tests

## 2. sign_vectors

### 2.1. sign_vectors.sign_vectors

* [ ] `is_harmonious_to` should work for a vector as input.
* [ ] improve new implementation of `SignVector`
* [ ] is list/tuple of Sign objects faster than sets for support and positive support?
* [ ] Cython might not be that useful if length of sign vectors is small (e.g. <100)

### 2.3. sign_vectors.utility

* [x] Improve the implementation of `closure`
  * [x] use sets
  * We could improve the implementation by using Combinations of support.
* [ ] implement lower, upper, total closure
  * We might be able to reuse code.

## 3. Sphinx documentation

NOTE:
I adapted the files from here:
https://github.com/mkauers/ore_algebra

* [ ] Do we want to have tests in the documentation?
  * Currently, tests are in the documentation
