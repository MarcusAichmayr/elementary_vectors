r"""
Motzkin's transposition theorem

This module deals with systems of linear inequalities.
(see [Ben00]_)

.. [Ben00] A. Ben-Israel.
    Motzkin's transposition theorem, and the related theorems of Farkas, Gordan and Stiemke.
    2000.

Given are matrices ``A``, ``B`` and vectors ``b``, ``c``.
We are interested whether the system ``A x <= b, B x < c`` has a solution ``x``.
This can be checked using :func:`elementary_vectors.vectors_in_intervals.exists_vector`.
If no solution exists, this function can be used to certify the result.

To demonstrate this, consider the following example::

    sage: from elementary_vectors.transposition_theorem import *
    sage: A = matrix([[1, 2], [0, 1]])
    sage: B = matrix([[2, 3]])
    sage: b = vector([1, -1])
    sage: c = vector([2])
    sage: exists_vector(matrix.block([[A.T, B.T]]), bc_to_intervals(b, c))
    True

However, to certify existence of a solution, we need to consider the alternative system.
This system can be described by two matrices ``M1``, ``M2`` and two lists of intervals ``I1``, ``I2``::

    sage: M1 = alternative_AB_bc1_matrix(A, B, b, c)
    sage: I1 = alternative_AB_bc1_intervals(A, B, b, c)
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
    sage: M2 = alternative_AB_bc2_matrix(A, B, b, c)
    sage: I2 = alternative_AB_bc2_intervals(A, B, b, c)
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
    sage: exists_vector(M1.T, I1, certificate=True)
    (5, -2, 1, 0, 0, 2)
    sage: exists_vector(M2.T, I2, certificate=True)
    (-1, 0, 1, 0, 0, 0, 2)

The package offers a single function that certifies existence of a solution::

    sage: certify_AB_bc(A, B, b, c)
    (True, [(5, -2, 1, 0, 0, 2), (-1, 0, 1, 0, 0, 0, 2)])

We consider another example::

    sage: A = matrix([[1, 0], [1, 1]])
    sage: B = matrix([[-1, -1]])
    sage: b = vector([1, 0])
    sage: c = vector([0])
    sage: certify_AB_bc(A, B, b, c)
    (False, [(0, 1, 1)])
    
There is also a homogeneous version of Motzkin's transposition theorem.
Here, we have three matrices ``A``, ``B`` and ``C`` and deals with the system
``A x > 0, B x >= 0, C x = 0``::

    sage: A = matrix([[1, 2], [0, 1]])
    sage: B = matrix([[2, 3]])
    sage: C = matrix([[-1, 0]])
    sage: exists_vector(matrix.block([[A.T, B.T, C.T]]), intervals_ABC(A, B, C))
    True

The certify the result, we consider the alternative system
which consists of a single matrix and intervals::

    sage: M = alternative_ABC_matrix(A, B, C)
    sage: I = alternative_ABC_intervals(A, B, C)
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
    sage: exists_vector(M.T, I, certificate=True)
    (0, 1, -1, 0, -3, -1)
    sage: certify_ABC(A, B, C)
    (True, (0, 1, -1, 0, -3, -1))

We consider another example::

    sage: A = matrix([[1, 0], [0, 1]])
    sage: B = matrix([[2, -3]])
    sage: C = matrix([[-1, -1]])
    sage: certify_ABC(A, B, C)
    (False, (1, 1, 0, 1))
"""

#############################################################################
#  Copyright (C) 2023                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)        #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from .vectors_in_intervals import exists_vector, exists_orthogonal_vector, setup_interval
from .functions import elementary_vectors
from sage.matrix.constructor import matrix, zero_matrix, identity_matrix, ones_matrix
from sage.rings.infinity import Infinity


def b_to_intervals(b):
    r"""Converts the right hand side of ``A x <= b`` to intervals"""
    return [setup_interval(-Infinity, bi) for bi in b]


def c_to_intervals(c):
    r"""Converts the right hand side of ``B x < c`` to intervals"""
    return [setup_interval(-Infinity, ci, False, False) for ci in c]


def bc_to_intervals(b, c):
    r"""Converts both right hand sides of ``A x <= b, B x < c`` to intervals"""
    return b_to_intervals(b) + c_to_intervals(c)


def intervals_ABC(A, B, C):
    r"""Compute the intervals corresponding to ``A x > 0, B x >= 0, C x = 0``."""
    return (
        A.nrows() * [setup_interval(0, Infinity, False, False)] +
        B.nrows() * [setup_interval(0, Infinity)] +
        C.nrows() * [setup_interval(0, 0)]
    )


def exists_orthogonal_vector_inhomogeneous(v, b, c):
    r"""
    Return whether there exists an orthogonal vector ``z = (z_1, z_2)`` to ``v`` satisfying ``z_1 <= b, z_2 < c``.
    
    INPUT:
    
    - ``v`` -- a vector
    
    - ``b`` -- a vector
    
    - ``c`` -- a vector
    
    .. SEEALSO::
    
        :func:`elementary_vectors.vectors_in_intervals.exists_orthogonal_vector`
    
    TESTS::
    
        sage: from elementary_vectors.transposition_theorem import *
        sage: from elementary_vectors.vectors_in_intervals import *
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
            if scalarproduct <= 0 and any(v[k] != 0 for k in range(m1, m1 + m2)):
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
    
        :func:`elementary_vectors.vectors_in_intervals.exists_orthogonal_vector`
    
    TESTS::
    
        sage: from elementary_vectors.transposition_theorem import *
        sage: from elementary_vectors.vectors_in_intervals import *
        sage: A = matrix([[1, 2], [0, 1]])
        sage: B = matrix([[2, 3], [0, -1]])
        sage: C = matrix([[-1, 0], [1, 1]])
        sage: M = matrix.block([[A.T, B.T, C.T]])
        sage: I = intervals_ABC(A, B, C)
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


def alternative_AB_bc1_matrix(A, B, b, c):
    r"""
    Matrix of first alternative system for ``A x <= b, B x < c``
    
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


def alternative_AB_bc1_intervals(A, B, b, c):
    r"""
    Intervals of first alternative system for ``A x <= b, B x < c``
    
    Returns a list of intervals.
    """
    m_A, n = A.dimensions()
    m_B = B.nrows()

    I = (
        n * [setup_interval(0, 0)] +
        (m_A + m_B) * [setup_interval(0, 1)] +
        [setup_interval(0, 1, False, True)]
    )
    return I


def alternative_AB_bc2_matrix(A, B, b, c):
    r"""
    Matrix of second alternative system for ``A x <= b, B x < c``
    
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


def alternative_AB_bc2_intervals(A, B, b, c):
    r"""
    Intervals of second alternative system for ``A x <= b, B x < c``
    
    Returns a list of intervals.
    """
    m_A, n = A.dimensions()
    m_B = B.nrows()

    I = (
        n * [setup_interval(0, 0)] +
        (m_A + m_B + 1) * [setup_interval(0, 1)] +
        [setup_interval(0, 1, False, True)]
    )
    return I


def alternative_ABC_matrix(A, B, C):
    r"""
    Alternative system matrix for ``A x > 0, B x >= 0, C x = 0``
    
    Returns a matrix.
    
    .. SEEALSO::
    
        :func:`~alternative_ABC_intervals`
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


def alternative_ABC_intervals(A, B, C):
    r"""
    Alternative system intervals for ``A x > 0, B x >= 0, C x = 0``
    
    Returns a list of intervals.

    .. SEEALSO::
    
        :func:`~alternative_ABC_matrix`
    """
    m_A, n = A.dimensions()
    m_B = B.nrows()

    I = (
        n * [setup_interval(0, 0)] +
        (m_A + m_B) * [setup_interval(0, 1)] +
        [setup_interval(0, 1, False, True)]
    )
    return I


def certify_AB_bc(A, B, b, c):
    r"""
    Return whether the system ``A x <= b, B x < c`` has a solution and certifies it.
    
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

    M1 = alternative_AB_bc1_matrix(A, B, b, c)
    M2 = alternative_AB_bc2_matrix(A, B, b, c)
    
    evs = elementary_vectors(M.T, generator=True)
    evs_alt1 = elementary_vectors(M1.T, generator=True)
    evs_alt2 = elementary_vectors(M2.T, generator=True)

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


def certify_ABC(A, B, C):
    r"""
    Return whether the system ``A x < 0, B x <= 0, C x = 0`` has a solution and certifies it.
    
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
    
    M_alt = alternative_ABC_matrix(A, B, C)

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
