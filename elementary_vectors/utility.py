"""Utility functions."""

#############################################################################
#  Copyright (C) 2022                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.symbolic.ring import SR
from sign_vectors import sign_vector
from sage.modules.free_module_element import vector


def sign_determined(a):
    r"""
    Check whether the sign of a number or symbolic expression ``a`` is uniquely determined.

    EXAMPLES::

        sage: from elementary_vectors.utility import sign_determined

    Integers have always a unique sign::

        sage: sign_determined(2)
        True
        sage: sign_determined(-5)
        True

    Now, we consider a variable::

        sage: var('a')
        a
        sage: sign_determined(a)
        False
        sage: assume(a >= 0)
        sage: sign_determined(a)
        False
        sage: assume(a != 0)
        sage: sign_determined(a)
        True
        sage: sign_determined(a - 1)
        False
    """
    return bool(SR(a) > 0 or SR(a) < 0 or SR(a) == 0)


def vector_from_matrix(M, I):
    r"""
    Return a vector in the right kernel of ``M`` such that
    the support is a subset of ``I``.

    INPUT:

    - ``M`` -- a matrix

    - ``I`` -- a list of indices

    OUTPUT:
    a vector ``v`` in the right kernel of ``M`` such that
    the support is a subset of ``I``.

    EXAMPLES::

        sage: from elementary_vectors.utility import vector_from_matrix
        sage: M = matrix([[1, 2, 0, 0], [0, 1, -1, 0]])
        sage: vector_from_matrix(M, [0, 1, 2])
        (2, -1, -1, 0)
        sage: vector_from_matrix(M, [3])
        (0, 0, 0, 1)
        sage: vector_from_matrix(M, [0, 3])
        (0, 0, 0, 1)
    """
    n = M.ncols()
    M_I = M.matrix_from_columns(I)
    try:
        l = list(M_I.right_kernel_matrix()[0])
    except IndexError:
        raise(ValueError("Right kernel of ``M`` restricted to the columns ``I`` is empty."))
    for k in range(n):
        if not k in I:
            l.insert(k, 0)
    return vector(l)


def conformal_elimination(x, y, S=None):
    r"""
    Apply conformal elimination to two real vectors to find a new vector.

    INPUT:

    - ``x`` -- a real vector

    - ``y`` -- a real vector

    - ``S`` -- a list of indices (default: ``[]``)

    OUTPUT:

    Returns a new vector ``z = x + a y`` where ``a > 0``, such that ``z[e] == 0``
    for some ``e`` in ``S`` and ``Z_S <= X_S`` and ``Z_f = (X o Y)_f`` for ``f``
    not in ``D(X, Y)``. Here, ``X``, ``Y`` and ``Z`` are the sign vectors
    corresponding to ``x``, ``y`` and ``z``.

    .. NOTE::

        If ``S`` is the empty list ``[]``, the whole list of separating elements
        will be considered instead. (default)
    """
    if S is None:
        S = []
    if x.length() != y.length():
        raise ValueError('Vectors have different length.')
    X = sign_vector(x)
    Y = sign_vector(y)
    D = X.separating_elements(Y)
    if D == []:
        raise ValueError('List of separating elements is empty.')
    if S == []:
        S = D
    elif not all(s in D for s in S):
        raise ValueError('S is not a subset of D.')
    lam = max([x[e]/y[e] for e in S])  # x[e]/y[e] < 0 since e in D
    return x - lam*y
