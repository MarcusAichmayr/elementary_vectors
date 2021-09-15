# TODO

## elementary_vectors

### elementary_vectors.elementary_vectors

* Vielleicht für `elementary_vectors` mehr optionale Argumente, um reductions besser zu kontrollieren?
    * optionales Argument: `cancel_common_factor` oder so
* Tests für verschiedene Ringe RR, QQ

* Output if ``return_minors=True`` is ``[m, evs]``.
  Should this be ``[evs, m]`` instead?

* Positive Elementary Vectors
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

#### Docstrings

* [ ] Add example with ``return_minors=True`` in `elementary_vectors` and `elementary_vectors_from_matrix`.

* [ ] Add examples to ``non_negative_elementary_vectors``.
* [ ] Add examples to ``positive_elementary_vectors``.


#### Tests

* [ ] Add tests to ``non_negative_elementary_vectors``.
* [ ] Add tests to ``positive_elementary_vectors``.

### elementary_vectors.utility

* `has_sign`
    * [ ] improve name of `has_sign`
    * [ ] add docstring
    * [ ] add tests

* `conformal_elimination`
    * Should this be here?
    * [ ] add examples
    * [ ] add tests

### elementary_vectors.reductions

* [ ] improve docstrings
    * [ ] add examples

### elementary_vectors.exists_vector

* [ ] If True and certificate is True, compute a corresponding vector of the vector space. (witness)
* [ ] add examples


## sign_vectors

### sign_vectors.sign_vectors

* [x] add docstrings

### sign_vectors.utility

* [ ] add examples
* [ ] move closure, contraction, deletion into this file

#### contraction, deletion

* [ ] add examples
* [ ] add tests

#### closure

* [ ] improve docstring
    * [ ] add examples
    * [ ] define closure
* [ ] add tests

### sign_vectors.oriented_matroids

* [ ] add examples
* [ ] add tests
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


