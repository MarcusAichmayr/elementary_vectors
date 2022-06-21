# TODO

<!-- TOC -->
- [1. elementary_vectors](#1-elementary_vectors)
  - [1.1. elementary_vectors.functions](#11-elementary_vectorsfunctions)
    - [1.1.1. elementary_vectors.functions.positive_elementary_vectors](#111-elementary_vectorsfunctionspositive_elementary_vectors)
  - [1.2. elementary_vectors.functions_dd](#12-elementary_vectorsfunctions_dd)
  - [1.3. elementary_vectors.vectors_in_intervals](#13-elementary_vectorsvectors_in_intervals)
    - [1.3.1. elementary_vectors.vectors_in_intervals.exists_vector](#131-elementary_vectorsvectors_in_intervalsexists_vector)
  - [1.4. elementary_vectors.utility](#14-elementary_vectorsutility)
  - [1.5. elementary_vectors.reductions](#15-elementary_vectorsreductions)
- [2. sign_vectors](#2-sign_vectors)
  - [2.1. sign_vectors.sign_vectors](#21-sign_vectorssign_vectors)
  - [2.2. sign_vectors.oriented_matroids](#22-sign_vectorsoriented_matroids)
  - [2.3. sign_vectors.utility](#23-sign_vectorsutility)
- [3. Sphinx documentation](#3-sphinx-documentation)
<!-- /TOC -->

## 1. elementary_vectors

* [ ] Add explanations and examples to demonstrate functions in this module. (in `__init__`?)
  * several reductions
  * `elementary_vectors`
  * `exists_vector`
* [ ] Should we rename this file?
  * It is a bit difficult to get the documentation of this file.
  Compare:
  ```
  import elementary_vectors
  elementary_vectors?

  from elementary_vectors import *
  elementary_vectors?
  ```

### 1.1. elementary_vectors.functions

* [x] Add examples to `positive_elementary_vectors`.
  * [ ] add examples with variables
* Do we want more optional arguments for `elementary_vectors` to have more control on reductions?
  * optional argument: `cancel_common_factor`?

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

CoCalc:
* Put functions for elementary vectors in one file. For maintaining, it would be easier to work more with the documentation.

### 1.2. elementary_vectors.functions_dd

* `determine_sign(X, a, M)`
  * Should this function work for real vectors (as `X` and `a`)? Does not make much sense.
  * [ ] This function works for cocircuits. For other sign vectors, it might fail.

### 1.3. elementary_vectors.vectors_in_intervals

* [ ] Right now, this works only for matrices.
  What if elementary vectors are given? Can we use those to construct a vector in the corresponding vector space?

#### 1.3.1. elementary_vectors.vectors_in_intervals.exists_vector

* [x] Removed optional argument `kernel`
* [x] Removed optional argument `certificate`
  * In case, someone wants to find the certificate, they can iterate over all elementary vectors and return the vector that does not satisfy `exists_normal_vector`.

### 1.4. elementary_vectors.utility

* `vector_from_matrix`
  * [ ] rename this function

* `conformal_elimination`
  * [ ] Should this be here?
  * [ ] add examples
  * [ ] add tests

### 1.5. elementary_vectors.reductions

* [ ] add examples to docstring on top
* [ ] (optional) write functions `reduce_factor_of_vector` and `reduce_factor_of_list`
  * [ ] `reduce_factor` should use these functions
* [ ] improve `reduce_factor` by canceling denominators if applicable


## 2. sign_vectors

### 2.1. sign_vectors.sign_vectors

* [ ] `is_harmonious` should work for a vector as input.
* [ ] improve new implementation of `SignVector`

### 2.2. sign_vectors.oriented_matroids

* [ ] The default value of `kernel` should be true for all functions in this module.
* [ ] Is it useful to add a function `cocircuits_from_elementary_vectors`?
  It would be used in `cocircuits_from_matrix` and `adjacent` of `double_description`.
* [ ] Is it possible to move references to the end of the documentation?
* [ ] It might be useful to use a class for oriented matroids.
  * sign vectors could be stored in a better way

### 2.3. sign_vectors.utility

* [ ] Can we improve the implementation of `closure`?
* [ ] improve docstring of `adjacent`
  * define adjacency.

## 3. Sphinx documentation

Note: I adapted the files from here: https://github.com/mkauers/ore_algebra

* Do we want to have tests in the documentation?
