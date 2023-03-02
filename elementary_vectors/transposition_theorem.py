r"""
Motzkin's transposition theorem

This module deals with systems of linear inequalities.
(see [Ben00]_)

.. [Ben00] Ben-Israel, A.:
    „Motzkin's transposition theorem, and the related theorems of farkas, gordan and stiemke“.

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

    sage: M1, I1 = alternative_AB_bc1(A, B, b, c)
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
    sage: M2, I2 = alternative_AB_bc2(A, B, b, c)
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

    sage: M, I = alternative_ABC(A, B, C)
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
"""

#############################################################################
#  Copyright (C) 2023                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from .vectors_in_intervals import exists_vector, setup_interval
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


def alternative_AB_bc1(A, B, b, c):
    r"""
    First alternative system for ``A x <= b, B x < c``
    
    Returns a matrix and corresponding intervals.
    """
    m_A, n = A.dimensions()
    m_B = B.nrows()

    M = matrix.block([
        [A.T, B.T],
        [identity_matrix(m_A), zero_matrix(m_A, m_B)],
        [zero_matrix(m_B, m_A), identity_matrix(m_B)],
        [matrix(1, -b), matrix(1, -c)] # number of rows is relevant for 0-dim vectors
    ])
    I = (
        n * [setup_interval(0, 0)] +
        (m_A + m_B) * [setup_interval(0, 1)] +
        [setup_interval(0, 1, False, True)]
    )
    return M, I


def alternative_AB_bc2(A, B, b, c):
    r"""
    Second alternative system for ``A x <= b, B x < c``
    
    Returns a matrix and corresponding intervals.
    """
    m_A, n = A.dimensions()
    m_B = B.nrows()

    M = matrix.block([
        [A.T, B.T],
        [identity_matrix(m_A), zero_matrix(m_A, m_B)],
        [zero_matrix(m_B, m_A), identity_matrix(m_B)],
        [matrix(1, -b), matrix(1, -c)], # number of rows is relevant for 0-dim vectors
        [zero_matrix(1, m_A), ones_matrix(1, m_B)]
    ])
    I = (
        n * [setup_interval(0, 0)] +
        (m_A + m_B + 1) * [setup_interval(0, 1)] +
        [setup_interval(0, 1, False, True)]
    )
    return M, I


def intervals_ABC(A, B, C):
    r"""Compute the intervals corresponding to ``A x > 0, B x >= 0, C x = 0``."""
    return (
        A.nrows() * [setup_interval(0, 1, False, True)] +
        B.nrows() * [setup_interval(0, 1)] +
        C.nrows() * [setup_interval(0, 0)]
    )


def alternative_ABC(A, B, C):
    r"""
    Alternative system for ``A x > 0, B x >= 0, C x = 0``
    
    Returns a matrix and corresponding intervals.
    """
    m_A, n = A.dimensions()
    m_B = B.nrows()
    m_C = C.nrows()

    M = matrix.block([
        [A.T, B.T, C.T],
        [identity_matrix(m_A), zero_matrix(m_A, m_B), zero_matrix(m_A, m_C)],
        [zero_matrix(m_B, m_A), identity_matrix(m_B), zero_matrix(m_B, m_C)],
        [ones_matrix(1, m_A), zero_matrix(1, m_B), zero_matrix(1, m_C)]
    ])
    I = (
        n * [setup_interval(0, 0)] +
        (m_A + m_B) * [setup_interval(0, 1)] +
        [setup_interval(0, 1, False, True)]
    )
    return M, I
