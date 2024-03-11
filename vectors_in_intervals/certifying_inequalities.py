r"""
Certifying solvability of linear inequality systems using elementary vectors
============================================================================

Given is a system of linear inequalities.
We are interested whether the system has a solution.
Further, we want to certify the result using an elementary vector (support-minimal vector).

Many of the results are derived from theorems in [Ben00]_.

.. [Ben00] A. Ben-Israel.
    Motzkin's transposition theorem, and the related theorems of Farkas, Gordan and Stiemke.
    2000.

Inhomogeneous systems
~~~~~~~~~~~~~~~~~~~~~

We deal with linear inequality systems given in the form
:math:`A x \leq b, B x < c` with matrices :math:`A`, :math:`B` and vectors :math:`b`, :math:`c`.
We are interested whether the system  has a solution :math:`x`.
This can be checked using :func:`vectors_in_intervals.existence.exists_vector`.
If no solution exists, this function can be used to certify the result.

To demonstrate this, consider the following example::

    sage: from vectors_in_intervals import *
    sage: from vectors_in_intervals.certifying_inequalities import *
    sage: A = matrix([[1, 2], [0, 1]])
    sage: B = matrix([[2, 3]])
    sage: b = vector([1, -1])
    sage: c = vector([2])
    sage: exists_vector(matrix_inhomogeneous(A, B).T, bc_to_intervals(b, c))
    True

However, to certify existence of a solution, we need to consider the alternative system.
This system can be described by two matrices and two lists of intervals::

    sage: M1 = matrix_inhomogeneous1_alternative(A, B, b, c)
    sage: I1 = intervals_inhomogeneous1_alternative(A, B, b, c)
    sage: M1
    [ 1  0| 2]
    [ 2  1| 3]
    [-----+--]
    [ 1  0| 0]
    [ 0  1| 0]
    [-----+--]
    [ 0  0| 1]
    [-----+--]
    [-1  1|-2]
    sage: I1
    [{0}, {0}, [0, 1], [0, 1], [0, 1], (0, 1]]
    sage: M2 = matrix_inhomogeneous2_alternative(A, B, b, c)
    sage: I2 = intervals_inhomogeneous2_alternative(A, B, b, c)
    sage: M2
    [ 1  0| 2]
    [ 2  1| 3]
    [-----+--]
    [ 1  0| 0]
    [ 0  1| 0]
    [-----+--]
    [ 0  0| 1]
    [-----+--]
    [-1  1|-2]
    [-----+--]
    [ 0  0| 1]
    sage: I2
    [{0}, {0}, [0, 1], [0, 1], [0, 1], [0, 1], (0, 1]]
    sage: exists_vector(M1.T, I1)
    False
    sage: exists_vector(M2.T, I2)
    False
    sage: exists_vector(M1.T, I1, certify=True)
    (5, -2, 1, 0, 0, 2)
    sage: exists_vector(M2.T, I2, certify=True)
    (-1, 0, 1, 0, 0, 0, 2)

The package offers a single function that certifies existence of a solution::

    sage: certify_inhomogeneous(A, B, b, c)
    (True, [(5, -2, 1, 0, 0, 2), (-1, 0, 1, 0, 0, 0, 2)])

We consider another example::

    sage: A = matrix([[1, 0], [1, 1]])
    sage: B = matrix([[-1, -1]])
    sage: b = vector([1, 0])
    sage: c = vector([0])
    sage: certify_inhomogeneous(A, B, b, c)
    (False, [(0, 1, 1)])
    
Homogeneous systems
~~~~~~~~~~~~~~~~~~~

There is also a homogeneous version of Motzkin's transposition theorem.
Here, we have three matrices :math:`A`, :math:`B` and :math:`C` and deals with the system
:math:`A x > 0, B x \geq 0, C x = 0`::

    sage: A = matrix([[1, 2], [0, 1]])
    sage: B = matrix([[2, 3]])
    sage: C = matrix([[-1, 0]])
    sage: exists_vector(matrix_homogeneous(A, B, C).T, intervals_homogeneous(A, B, C))
    True

The certify the result, we consider the alternative system
which consists of a single matrix and intervals::

    sage: M = matrix_homogeneous_alternative(A, B, C)
    sage: I = intervals_homogeneous_alternative(A, B, C)
    sage: M
    [ 1  0| 2|-1]
    [ 2  1| 3| 0]
    [-----+--+--]
    [ 1  0| 0| 0]
    [ 0  1| 0| 0]
    [-----+--+--]
    [ 0  0| 1| 0]
    [-----+--+--]
    [ 1  1| 0| 0]
    sage: I
    [{0}, {0}, [0, 1], [0, 1], [0, 1], (0, 1]]
    sage: exists_vector(M.T, I)
    False
    sage: exists_vector(M.T, I, certify=True)
    (0, 1, -1, 0, -3, -1)
    sage: certify_homogeneous(A, B, C)
    (True, (0, 1, -1, 0, -3, -1))

We consider another example::

    sage: A = matrix([[1, 0], [0, 1]])
    sage: B = matrix([[2, -3]])
    sage: C = matrix([[-1, -1]])
    sage: certify_homogeneous(A, B, C)
    (False, (1, 1, 0, 1))
"""

#############################################################################
#  Copyright (C) 2024                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)        #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

import warnings
from sage.combinat.combination import Combinations
from sage.matrix.constructor import matrix, zero_matrix, identity_matrix, ones_matrix
from sage.modules.free_module_element import zero_vector
from sage.rings.infinity import Infinity

from elementary_vectors import elementary_vectors
from elementary_vectors.reductions import reduce_vectors_support
from .utility import interval_from_bounds


def matrix_inhomogeneous(A, B, b=None, c=None):
    r"""Return a block matrix for the system ``A x <= b, B x < c``."""
    return matrix.block([[A], [B]])


def matrix_homogeneous(A, B, C):
    r"""Return a block matrix for the system ``A x > 0, B x >= 0, C x = 0``."""
    return matrix.block([[A], [B], [C]])


def b_to_intervals(b):
    r"""Converts the right hand side of ``A x <= b`` to intervals."""
    return [interval_from_bounds(-Infinity, bi) for bi in b]


def c_to_intervals(c):
    r"""Converts the right hand side of ``B x < c`` to intervals."""
    return [interval_from_bounds(-Infinity, ci, False, False) for ci in c]


def bc_to_intervals(b, c):
    r"""Converts both right hand sides of ``A x <= b, B x < c`` to intervals."""
    return b_to_intervals(b) + c_to_intervals(c)


def intervals_homogeneous(A, B, C):
    r"""Compute the intervals corresponding to ``A x > 0, B x >= 0, C x = 0``."""
    return (
        A.nrows() * [interval_from_bounds(0, Infinity, False, False)] +
        B.nrows() * [interval_from_bounds(0, Infinity)] +
        C.nrows() * [interval_from_bounds(0, 0)]
    )


def exists_orthogonal_vector_inhomogeneous(v, b, c):
    r"""
    Return whether there exists an orthogonal vector ``z = (z_1, z_2)`` to ``v`` satisfying ``z_1 <= b, z_2 < c``.
    
    INPUT:
    
    - ``v`` -- a vector
    
    - ``b`` -- a vector
    
    - ``c`` -- a vector
    
    .. SEEALSO::
    
        :func:`vectors_in_intervals.existence.exists_orthogonal_vector`
    
    TESTS::
    
        sage: from vectors_in_intervals.certifying_inequalities import *
        sage: from vectors_in_intervals import *
        sage: A = matrix([[1, 2], [0, 1], [2, 0]])
        sage: B = matrix([[2, 3], [0, 2]])
        sage: b = vector([1, -1, 2])
        sage: c = vector([2, 3])
        sage: I = bc_to_intervals(b, c)
        sage: M = matrix.block([[A.T, B.T]])
        sage: for v in elementary_vectors(M, generator=True):
        ....:     assert(exists_orthogonal_vector(v, I) == exists_orthogonal_vector_inhomogeneous(v, b, c))
    """
    m1 = len(b)
    m2 = len(c)
    def condition(v):
        if all(vk >= 0 for vk in v):
            scalarproduct = (
                sum(v[k] * b[k] for k in range(m1))
                +
                sum(v[k + m1] * c[k] for k in range(m2))
            )
            if scalarproduct < 0:
                return True
            if scalarproduct <= 0 and any(v[k] for k in range(m1, m1 + m2)):
                return True
        return False
    return not condition(v) and not condition(-v)


def exists_orthogonal_vector_homogeneous(v, range_strict, range_non_strict):
    r"""
    Return whether there exists an orthogonal vector to ``v`` lying by homogeneous intervals.
    
    INPUT:
    
    - ``v`` -- a vector
    
    - ``range_strict`` -- range of strict inequalities
    
    - ``range_non_strict`` -- range of non-strict inequalities
    
    OUTPUT:
    Checks if there exists an orthogonal vector to ``v`` that lies in homogeneous intervals.
    The intervals are ``(0, oo)`` for components in ``range_strict``,
    ``[0, oo)`` for components in ``range_non_strict``
    and ``{0}`` for all remaining components.

    .. SEEALSO::
    
        :func:`vectors_in_intervals.existence.exists_orthogonal_vector`
    
    TESTS::
    
        sage: from vectors_in_intervals.certifying_inequalities import *
        sage: from vectors_in_intervals import *
        sage: A = matrix([[1, 2], [0, 1]])
        sage: B = matrix([[2, 3], [0, -1]])
        sage: C = matrix([[-1, 0], [1, 1]])
        sage: M = matrix.block([[A.T, B.T, C.T]])
        sage: I = intervals_homogeneous(A, B, C)
        sage: for v in elementary_vectors(M, generator=True):
        ....:     assert(exists_orthogonal_vector(v, I) == exists_orthogonal_vector_homogeneous(v, range(A.nrows()), range(A.nrows(), A.nrows() + B.nrows())))
    """
    return not (
        any(v[k] != 0 for k in range_strict)
        and
        (
            (
                all(v[k] >= 0 for k in range_strict)
                and
                all(v[k] >= 0 for k in range_non_strict)
            )
            or
            (
                all(v[k] <= 0 for k in range_strict)
                and
                all(v[k] <= 0 for k in range_non_strict)
            )
        )
    )


def elementary_vectors_generator_trailing_nonzero(M):
    r"""
    Return generator of elementary vectors that are non-zero at last component.
    
    INPUT:
    
    - ``M`` -- a matrix
    
    .. SEEALSO::
    
        :func:`elementary_vectors.functions.elementary_vectors`
    """
    try:
        ind = M.pivot_rows()
        M = M.matrix_from_rows(ind)
    except (ArithmeticError, NotImplementedError):
        warnings.warn('Could not determine rank of matrix. Expect wrong result!')

    m, n = M.dimensions()
    minors = {}
    ring = M.base_ring()

    def ev_from_support(indices):
        r"""
        Return the elementary vector corresponding to ``indices``.

        INPUT:

        - ``indices`` -- a list of indices
        """
        v = zero_vector(ring, n)
        for pos, k in enumerate(indices):
            Ik = tuple(i for i in indices if i != k)
            try:
                minor = minors[Ik]
            except KeyError:
                minors[Ik] = M.matrix_from_columns(Ik).det()
                minor = minors[Ik]
            v[k] = (-1)**pos * minor
        return v

    evs = (ev_from_support(indices + [n - 1]) for indices in Combinations(n - 1, m))
    evs = (v for v in evs if v)
    evs = reduce_vectors_support(evs, generator=True)

    return evs


def matrix_inhomogeneous1_alternative(A, B, b, c):
    r"""
    Matrix of first alternative system for ``A x <= b, B x < c``.
    
    Returns a matrix.
    """
    m_A = A.nrows()
    m_B = B.nrows()

    M = matrix.block([
        [A.T, B.T],
        [identity_matrix(m_A), zero_matrix(m_A, m_B)],
        [zero_matrix(m_B, m_A), identity_matrix(m_B)],
        [matrix(1, -b), matrix(1, -c)] # number of rows is relevant for 0-dim vectors
    ])
    return M


def matrix_inhomogeneous2_alternative(A, B, b, c):
    r"""
    Matrix of second alternative system for ``A x <= b, B x < c``.
    
    Returns a matrix.
    """
    m_A = A.nrows()
    m_B = B.nrows()

    M = matrix.block([
        [A.T, B.T],
        [identity_matrix(m_A), zero_matrix(m_A, m_B)],
        [zero_matrix(m_B, m_A), identity_matrix(m_B)],
        [matrix(1, -b), matrix(1, -c)], # number of rows is relevant for 0-dim vectors
        [zero_matrix(1, m_A), ones_matrix(1, m_B)]
    ])
    return M


def matrix_homogeneous_alternative(A, B, C):
    r"""
    Alternative system matrix for ``A x > 0, B x >= 0, C x = 0``.
    
    Returns a matrix.
    
    .. SEEALSO::
    
        :func:`~intervals_homogeneous_alternative`
    """
    m_A = A.nrows()
    m_B = B.nrows()
    m_C = C.nrows()

    M = matrix.block([
        [A.T, B.T, C.T],
        [identity_matrix(m_A), zero_matrix(m_A, m_B), zero_matrix(m_A, m_C)],
        [zero_matrix(m_B, m_A), identity_matrix(m_B), zero_matrix(m_B, m_C)],
        [ones_matrix(1, m_A), zero_matrix(1, m_B), zero_matrix(1, m_C)]
    ])
    return M


def intervals_inhomogeneous1_alternative(A, B, b=None, c=None):
    r"""
    Intervals of first alternative system for ``A x <= b, B x < c``
    
    Return a list of intervals.
    """
    m_A, n = A.dimensions()
    m_B = B.nrows()

    return (
        n * [interval_from_bounds(0, 0)] +
        (m_A + m_B) * [interval_from_bounds(0, 1)] +
        [interval_from_bounds(0, 1, False, True)]
    )


def intervals_inhomogeneous2_alternative(A, B, b=None, c=None):
    r"""
    Intervals of second alternative system for ``A x <= b, B x < c``.
    
    Return a list of intervals.
    """
    m_A, n = A.dimensions()
    m_B = B.nrows()

    return (
        n * [interval_from_bounds(0, 0)] +
        (m_A + m_B + 1) * [interval_from_bounds(0, 1)] +
        [interval_from_bounds(0, 1, False, True)]
    )


def intervals_homogeneous_alternative(A, B, C=None):
    r"""
    Alternative system intervals for ``A x > 0, B x >= 0, C x = 0``.
    
    Return a list of intervals.

    .. SEEALSO::
    
        :func:`~matrix_homogeneous_alternative`
    """
    m_A, n = A.dimensions()
    m_B = B.nrows()

    return (
        n * [interval_from_bounds(0, 0)] +
        (m_A + m_B) * [interval_from_bounds(0, 1)] +
        [interval_from_bounds(0, 1, False, True)]
    )


def certify_inhomogeneous(A, B, b, c):
    r"""
    Return whether the system ``A x <= b, B x < c`` has a solution and certify it.
    
    INPUT:
    
    - ``A`` -- a matrix
    
    - ``B`` -- a matrix
    
    - ``b`` -- a vector
    
    - ``c`` -- a vector
    
    OUTPUT:
    A tuple of a boolean and a list of vectors certifying the result.
    """
    M = matrix.block([[A], [B]])

    n = A.ncols()
    m_A = A.nrows()
    m_B = B.nrows()

    M1 = matrix_inhomogeneous1_alternative(A, B, b, c)
    M2 = matrix_inhomogeneous2_alternative(A, B, b, c)

    evs = elementary_vectors(M.T, generator=True)
    evs_alt1 = elementary_vectors_generator_trailing_nonzero(M1.T)
    evs_alt2 = elementary_vectors_generator_trailing_nonzero(M2.T)

    evs_end_reached = False
    evs_alt1_end_reached = False
    evs_alt2_end_reached = False

    v1_found = False
    v2_found = False
    while True:
        if not evs_end_reached:
            try:
                v = next(evs)
                if not exists_orthogonal_vector_inhomogeneous(v, b, c):
                    return False, [v]
            except StopIteration:
                evs_end_reached = True

        if not evs_alt1_end_reached and not v1_found:
            try:
                v1 = next(evs_alt1)
                if not exists_orthogonal_vector_homogeneous(
                    v1,
                    [n + m_A + m_B],
                    range(n, n + m_A + m_B)
                ):
                    if v2_found:
                        return True, [v1, v2]
                    v1_found = True
            except StopIteration:
                evs_alt1_end_reached = True

        if not evs_alt2_end_reached and not v2_found:
            try:
                v2 = next(evs_alt2)
                if not exists_orthogonal_vector_homogeneous(
                    v2,
                    [n + m_A + m_B + 1],
                    range(n, n + m_A + m_B + 1)
                ):
                    if v1_found:
                        return True, [v1, v2]
                    v2_found = True
            except StopIteration:
                evs_alt2_end_reached = True


def certify_homogeneous(A, B, C):
    r"""
    Return whether the system ``A x < 0, B x <= 0, C x = 0`` has a solution and certify it.
    
    INPUT:
    
    - ``A`` -- a matrix
    
    - ``B`` -- a matrix
    
    - ``C`` -- a matrix
    
    OUTPUT:
    A tuple of a boolean and a list of vectors certifying the result.
    """
    m_A, n = A.dimensions()
    m_B = B.nrows()
    M = matrix.block([[A.T, B.T, C.T]])

    M_alt = matrix_homogeneous_alternative(A, B, C)

    evs = elementary_vectors(M, generator=True)
    evs_alt = elementary_vectors(M_alt.T, generator=True)

    evs_end_reached = False
    evs_alt_end_reached = False
    while True:
        if not evs_end_reached:
            try:
                v = next(evs)
                if not exists_orthogonal_vector_homogeneous(
                    v,
                    range(m_A),
                    range(m_A, m_A + m_B)
                ):
                    return False, v
            except StopIteration:
                evs_end_reached = True

        if not evs_alt_end_reached:
            try:
                v = next(evs_alt)
                if not exists_orthogonal_vector_homogeneous(
                    v,
                    [n + m_A + m_B],
                    range(n, n + m_A + m_B)
                ):
                    return True, v
            except StopIteration:
                evs_alt_end_reached = True
