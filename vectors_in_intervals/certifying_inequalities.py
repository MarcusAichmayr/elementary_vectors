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

Homogeneous systems
~~~~~~~~~~~~~~~~~~~

There is also a homogeneous version of Motzkin's transposition theorem.
Here, we have three matrices :math:`A`, :math:`B` and :math:`C` and deals with the system
:math:`A x > 0, B x \geq 0, C x = 0`::

    sage: from vectors_in_intervals import *
    sage: from vectors_in_intervals.certifying_inequalities import *
    sage: A = matrix([[1, 2], [0, 1]])
    sage: B = matrix([[2, 3]])
    sage: C = matrix([[-1, 0]])
    sage: S = HomogeneousSystem(A, B, C)
    sage: S.matrix
    [ 1  2]
    [ 0  1]
    [-----]
    [ 2  3]
    [-----]
    [-1  0]
    sage: S.intervals
    [(0, +oo), (0, +oo), [0, +oo), {0}]
    sage: exists_vector(S.matrix.T, S.intervals)
    True

The certify the result, we consider the alternative system
which consists of a single matrix and intervals::

    sage: S.matrix_alternative
    [ 1  0| 2|-1]
    [ 2  1| 3| 0]
    [-----+--+--]
    [ 1  0| 0| 0]
    [ 0  1| 0| 0]
    [-----+--+--]
    [ 0  0| 1| 0]
    [-----+--+--]
    [ 1  1| 0| 0]
    sage: S.intervals_alternative
    [{0}, {0}, [0, 1], [0, 1], [0, 1], (0, 1]]
    sage: exists_vector(S.matrix_alternative.T, S.intervals_alternative)
    False
    sage: exists_vector(S.matrix_alternative.T, S.intervals_alternative, certify=True)
    (0, 1, -1, 0, -3, -1)

There is a single command for certification::

    sage: S.certify()
    (True, (0, 1, -1, 0, -3, -1))

We consider another example::

    sage: A = matrix([[1, 0], [0, 1]])
    sage: B = matrix([[2, -3]])
    sage: C = matrix([[-1, -1]])
    sage: S = HomogeneousSystem(A, B, C)
    sage: S.certify()
    (False, (1, 1, 0, 1))

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
    sage: S = InhomogeneousSystem(A, B, b, c)
    sage: S.matrix()
    [1 2]
    [0 1]
    [---]
    [2 3]
    sage: S.intervals()
    [(-oo, 1], (-oo, -1], (-oo, 2)]
    sage: exists_vector(S.matrix().T, S.intervals())
    True

However, to certify existence of a solution, we need to consider the alternative system.
This system can be described by two matrices and two lists of intervals::

    sage: S.matrix_alternative()
    [ 1  0| 2| 0]
    [ 2  1| 3| 0]
    [-----+--+--]
    [-1  1|-2|-1]
    [-----+--+--]
    [ 1  0| 0| 0]
    [ 0  1| 0| 0]
    [-----+--+--]
    [ 0  0| 1| 0]
    [-----+--+--]
    [ 0  0| 0| 1]
    [-----+--+--]
    [ 0  0| 1| 1]
    sage: S.intervals_alternative()
    [{0}, {0}, {0}, [0, 1], [0, 1], [0, 1], [0, 1], (0, 1]]
    sage: exists_vector(S.matrix_alternative().T, S.intervals_alternative())
    False
    sage: exists_vector(S.matrix_alternative().T, S.intervals_alternative(), certify=True)
    (-4, 2, -2, -2, 0, 0, 0, -2)

The package offers a single function that certifies existence of a solution::

    sage: S.certify()
    (True, (-4, 2, -2, -2, 0, 0, 0, -2))

We consider another example::

    sage: A = matrix([[1, 0], [1, 1]])
    sage: B = matrix([[-1, -1]])
    sage: b = vector([1, 0])
    sage: c = vector([0])
    sage: S = InhomogeneousSystem(A, B, b, c)
    sage: S.certify()
    (False, [(0, 1, 1)])
"""

#############################################################################
#  Copyright (C) 2024                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

import warnings
from sage.combinat.combination import Combinations
from sage.matrix.constructor import matrix, zero_matrix, identity_matrix, ones_matrix
from sage.rings.infinity import Infinity
from sage.sets.real_set import RealSet
from sage.structure.sage_object import SageObject

from elementary_vectors import elementary_vectors
from vectors_in_intervals import exists_orthogonal_vector # TODO remove
from elementary_vectors.utility import elementary_vector_from_indices_prevent_multiples
from .utility import interval_from_bounds


def matrix_inhomogeneous(A, B, b=None, c=None):
    r"""Return a block matrix for the system ``A x <= b, B x < c``."""
    return matrix.block([[A], [B]])


def matrix_homogeneous(A, B, C):
    r"""Return a block matrix for the system ``A x > 0, B x >= 0, C x = 0``."""
    return matrix.block([[A], [B], [C]])


def b_to_intervals(b) -> list[RealSet]:
    r"""Converts the right hand side of ``A x <= b`` to intervals."""
    return [interval_from_bounds(-Infinity, bi) for bi in b]


def c_to_intervals(c) -> list[RealSet]:
    r"""Converts the right hand side of ``B x < c`` to intervals."""
    return [interval_from_bounds(-Infinity, ci, False, False) for ci in c]


def bc_to_intervals(b, c) -> list[RealSet]:
    r"""Converts both right hand sides of ``A x <= b, B x < c`` to intervals."""
    return b_to_intervals(b) + c_to_intervals(c)


def intervals_homogeneous(A, B, C) -> list[RealSet]:
    r"""Compute the intervals corresponding to ``A x > 0, B x >= 0, C x = 0``."""
    return (
        A.nrows() * [interval_from_bounds(0, Infinity, False, False)]
        + B.nrows() * [interval_from_bounds(0, Infinity)]
        + C.nrows() * [interval_from_bounds(0, 0)]
    )


def exists_orthogonal_vector_inhomogeneous(v, b, c) -> bool:
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
            scalarproduct = sum(v[k] * b[k] for k in range(m1)) + sum(
                v[k + m1] * c[k] for k in range(m2)
            )
            if scalarproduct < 0:
                return True
            if scalarproduct <= 0 and any(v[k] for k in range(m1, m1 + m2)):
                return True
        return False

    return not condition(v) and not condition(-v)


def exists_orthogonal_vector_homogeneous(v, range_strict, range_non_strict) -> bool:
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
        and (
            (
                all(v[k] >= 0 for k in range_strict)
                and all(v[k] >= 0 for k in range_non_strict)
            )
            or (
                all(v[k] <= 0 for k in range_strict)
                and all(v[k] <= 0 for k in range_non_strict)
            )
        )
    )


def elementary_vectors_generator_trailing_nonzero(M):
    r"""
    Return generator of elementary vectors that are nonzero at last component.

    INPUT:

    - ``M`` -- a matrix

    .. SEEALSO::

        :func:`elementary_vectors.functions.elementary_vectors`
    """
    try:
        M = M.matrix_from_rows(M.pivot_rows())
    except (ArithmeticError, NotImplementedError):
        warnings.warn("Could not determine rank of matrix. Expect wrong result!")

    rank, length = M.dimensions()
    minors = {}

    evs = (
        elementary_vector_from_indices_prevent_multiples(indices + [length - 1], minors, M)
        for indices in Combinations(length - 1, rank)
    )
    evs = (v for v in evs if v)

    return evs


def matrix_inhomogeneous1_alternative(A, B, b, c):
    r"""
    Matrix of first alternative system for ``A x <= b, B x < c``.

    Returns a matrix.
    """
    m_A = A.nrows()
    m_B = B.nrows()

    M = matrix.block(
        [
            [A.T, B.T],
            [identity_matrix(m_A), zero_matrix(m_A, m_B)],
            [zero_matrix(m_B, m_A), identity_matrix(m_B)],
            [
                matrix(1, -b),
                matrix(1, -c),
            ],  # number of rows is relevant for 0-dim vectors
        ]
    )
    return M


def matrix_inhomogeneous2_alternative(A, B, b, c):
    r"""
    Matrix of second alternative system for ``A x <= b, B x < c``.

    Returns a matrix.
    """
    m_A = A.nrows()
    m_B = B.nrows()

    M = matrix.block(
        [
            [A.T, B.T],
            [identity_matrix(m_A), zero_matrix(m_A, m_B)],
            [zero_matrix(m_B, m_A), identity_matrix(m_B)],
            [
                matrix(1, -b),
                matrix(1, -c),
            ],  # number of rows is relevant for 0-dim vectors
            [zero_matrix(1, m_A), ones_matrix(1, m_B)],
        ]
    )
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

    M = matrix.block(
        [
            [A.T, B.T, C.T],
            [identity_matrix(m_A), zero_matrix(m_A, m_B), zero_matrix(m_A, m_C)],
            [zero_matrix(m_B, m_A), identity_matrix(m_B), zero_matrix(m_B, m_C)],
            [ones_matrix(1, m_A), zero_matrix(1, m_B), zero_matrix(1, m_C)],
        ]
    )
    return M


def intervals_inhomogeneous1_alternative(A, B, b=None, c=None):
    r"""
    Intervals of first alternative system for ``A x <= b, B x < c``

    Return a list of intervals.
    """
    m_A, length = A.dimensions()
    m_B = B.nrows()

    return (
        length * [interval_from_bounds(0, 0)]
        + (m_A + m_B) * [interval_from_bounds(0, 1)]
        + [interval_from_bounds(0, 1, False, True)]
    )


def intervals_inhomogeneous2_alternative(A, B, b=None, c=None):
    r"""
    Intervals of second alternative system for ``A x <= b, B x < c``.

    Return a list of intervals.
    """
    m_A, length = A.dimensions()
    m_B = B.nrows()

    return (
        length * [interval_from_bounds(0, 0)]
        + (m_A + m_B + 1) * [interval_from_bounds(0, 1)]
        + [interval_from_bounds(0, 1, False, True)]
    )


def intervals_homogeneous_alternative(A, B, C=None):
    r"""
    Alternative system intervals for ``A x > 0, B x >= 0, C x = 0``.

    Return a list of intervals.

    .. SEEALSO::

        :func:`~matrix_homogeneous_alternative`
    """
    m_A, length = A.dimensions()
    m_B = B.nrows()

    return (
        length * [interval_from_bounds(0, 0)]
        + (m_A + m_B) * [interval_from_bounds(0, 1)]
        + [interval_from_bounds(0, 1, False, True)]
    )


def certify_inhomogeneous(A, B, b, c) -> tuple:
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

    length = A.ncols()
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
                    v1, [-1], range(length, length + m_A + m_B)
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
                    v2, [-1], range(length, length + m_A + m_B + 1)
                ):
                    if v1_found:
                        return True, [v1, v2]
                    v2_found = True
            except StopIteration:
                evs_alt2_end_reached = True


def certify_homogeneous(A, B, C) -> tuple:
    r"""
    Return whether the system ``A x < 0, B x <= 0, C x = 0`` has a solution and certify it.

    INPUT:

    - ``A`` -- a matrix

    - ``B`` -- a matrix

    - ``C`` -- a matrix

    OUTPUT:
    A tuple of a boolean and a list of vectors certifying the result.
    """
    m_A, length = A.dimensions()
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
                    v, range(m_A), range(m_A, m_A + m_B)
                ):
                    return False, v
            except StopIteration:
                evs_end_reached = True

        if not evs_alt_end_reached:
            try:
                v = next(evs_alt)
                if not exists_orthogonal_vector_homogeneous(
                    v, [-1], range(length, length + m_A + m_B)
                ):
                    return True, v
            except StopIteration:
                evs_alt_end_reached = True


class HomogeneousSystem(SageObject):
    r"""
    A class for setting up and certifying homogeneous linear inequality systems.

    The system is given ``A x > 0``, ``B x >= 0`` and ``C x = 0`` with matrices ``A``, ``B``, ``C``.

    EXAMPLES::

        sage: from vectors_in_intervals.certifying_inequalities import *
        sage: A = matrix([[1, 2], [0, 1]])
        sage: B = matrix([[2, 3]])
        sage: C = matrix([[-1, 0]])
        sage: S = HomogeneousSystem(A, B, C)
        sage: S.matrix
        [ 1  2]
        [ 0  1]
        [-----]
        [ 2  3]
        [-----]
        [-1  0]
        sage: S.intervals
        [(0, +oo), (0, +oo), [0, +oo), {0}]
        sage: S.matrix_alternative
        [ 1  0| 2|-1]
        [ 2  1| 3| 0]
        [-----+--+--]
        [ 1  0| 0| 0]
        [ 0  1| 0| 0]
        [-----+--+--]
        [ 0  0| 1| 0]
        [-----+--+--]
        [ 1  1| 0| 0]
        sage: S.intervals_alternative
        [{0}, {0}, [0, 1], [0, 1], [0, 1], (0, 1]]
        sage: S.certify()
        (True, (0, 1, -1, 0, -3, -1))

    We consider another example::

        sage: A = matrix([[1, 0], [0, 1]])
        sage: B = matrix([[2, -3]])
        sage: C = matrix([[-1, -1]])
        sage: S = HomogeneousSystem(A, B, C)
        sage: S.certify()
        (False, (1, 1, 0, 1))
    """
    def __init__(self, A, B, C) -> None:
        r"""
        INPUT:
        matrices with same number of columns
        """

        self.A = A
        self.B = B
        self.C = C
        self.matrix = matrix.block([[A], [B], [C]])
        m_A = A.nrows()
        m_B = B.nrows()
        m_C = C.nrows()
        length = A.ncols()
        self.intervals = (
            m_A * [interval_from_bounds(0, Infinity, False, False)]
            + m_B * [interval_from_bounds(0, Infinity)]
            + m_C * [interval_from_bounds(0, 0)]
        )
        self.matrix_alternative = matrix.block(
            [
                [A.T, B.T, C.T],
                [identity_matrix(m_A), zero_matrix(m_A, m_B), zero_matrix(m_A, m_C)],
                [zero_matrix(m_B, m_A), identity_matrix(m_B), zero_matrix(m_B, m_C)],
                [ones_matrix(1, m_A), zero_matrix(1, m_B), zero_matrix(1, m_C)],
            ]
        )
        self.intervals_alternative = (
            length * [interval_from_bounds(0, 0)]
            + (m_A + m_B) * [interval_from_bounds(0, 1)]
            + [interval_from_bounds(0, 1, False, True)]
        )

    def certify(self) -> tuple:
        r"""
        Return whether this system has a solution and certify it.

        OUTPUT:
        A tuple of a boolean and a list of vectors certifying the result.
        """

        m_A, length = self.A.dimensions()
        m_B = self.B.nrows()

        evs = elementary_vectors(self.matrix.T, generator=True)
        evs_alt = elementary_vectors(self.matrix_alternative.T, generator=True)

        evs_end_reached = False
        evs_alt_end_reached = False
        while True:
            if not evs_end_reached:
                try:
                    v = next(evs)
                    if not exists_orthogonal_vector_homogeneous(
                        v, range(m_A), range(m_A, m_A + m_B)
                    ):
                        return False, v
                except StopIteration:
                    evs_end_reached = True

            if not evs_alt_end_reached:
                try:
                    v = next(evs_alt)
                    if not exists_orthogonal_vector_homogeneous(
                        v, [-1], range(length, length + m_A + m_B)
                    ):
                        return True, v
                except StopIteration:
                    evs_alt_end_reached = True


class InhomogeneousSystem(SageObject):
    r"""
    A class for setting up and certifying an inhomogeneous linear inequality systems.

    The system is given by ``A x <= b`` and ``B x < c`` with matrices ``A``, ``B`` and vectors ``b``, ``c``.

    EXAMPLES::

        sage: from vectors_in_intervals import *
        sage: from vectors_in_intervals.certifying_inequalities import *
        sage: A = matrix([[1, 2], [0, 1]])
        sage: B = matrix([[2, 3]])
        sage: b = vector([1, -1])
        sage: c = vector([2])
        sage: S = InhomogeneousSystem(A, B, b, c)
        sage: S.matrix()
        [1 2]
        [0 1]
        [---]
        [2 3]
        sage: S.intervals()
        [(-oo, 1], (-oo, -1], (-oo, 2)]
        sage: S.matrix_alternative()
        [ 1  0| 2| 0]
        [ 2  1| 3| 0]
        [-----+--+--]
        [-1  1|-2|-1]
        [-----+--+--]
        [ 1  0| 0| 0]
        [ 0  1| 0| 0]
        [-----+--+--]
        [ 0  0| 1| 0]
        [-----+--+--]
        [ 0  0| 0| 1]
        [-----+--+--]
        [ 0  0| 1| 1]
        sage: S.intervals_alternative()
        [{0}, {0}, {0}, [0, 1], [0, 1], [0, 1], [0, 1], (0, 1]]
        sage: S.certify()
        (True, (-4, 2, -2, -2, 0, 0, 0, -2))
        sage: S.certify_using_two_systems()
        (True, [(5, -2, 1, 0, 0, 2), (-1, 0, 1, 0, 0, 0, 2)])

    We consider another example::

        sage: A = matrix([[1, 0], [1, 1]])
        sage: B = matrix([[-1, -1]])
        sage: b = vector([1, 0])
        sage: c = vector([0])
        sage: S = InhomogeneousSystem(A, B, b, c)
        sage: S.certify()
        (False, [(0, 1, 1)])
    """

    def __init__(self, A, B, b, c) -> None:
        r"""
        INPUT:

        - ``A`` -- a ``p`` times ``n`` matrix

        - ``B`` -- a ``q`` times ``n`` matrix

        - ``b`` -- a ``p`` vector

        - ``c`` -- a ``q`` vector
        """
        self.A = A
        self.B = B
        self.b = b
        self.c = c

    def matrix(self):
        if not hasattr(self, "_matrix"):
            self._matrix = matrix.block([[self.A], [self.B]])
        return self._matrix

    def intervals(self):
        if not hasattr(self, "_intervals"):
            self._intervals = (
                [interval_from_bounds(-Infinity, bi) for bi in self.b]
                + [interval_from_bounds(-Infinity, ci, False, False) for ci in self.c]
            )
        return self._intervals

    def matrix_alternative(self):
        if not hasattr(self, "_matrix_alternative"):
            m_A = self.A.nrows()
            m_B = self.B.nrows()
            length = self.A.ncols()
            self._matrix_alternative = matrix.block(
                [
                    [self.A.T, self.B.T, zero_matrix(length, 1)],
                    [matrix(1, -self.b), matrix(1, -self.c), matrix([[-1]])],
                    [identity_matrix(m_A), zero_matrix(m_A, m_B), zero_matrix(m_A, 1)],
                    [zero_matrix(m_B, m_A), identity_matrix(m_B), zero_matrix(m_B, 1)],
                    [zero_matrix(1, m_A), zero_matrix(1, m_B), matrix([[1]])],
                    [zero_matrix(1, m_A), ones_matrix(1, m_B), matrix([[1]])],
                ]
            )
        return self._matrix_alternative

    def intervals_alternative(self) -> list:
        if not hasattr(self, "_intervals_alternative"):
            m_A = self.A.nrows()
            m_B = self.B.nrows()
            length = self.A.ncols()
            self._intervals_alternative = (
                (length + 1) * [interval_from_bounds(0, 0)]
                + (m_A + m_B + 1) * [interval_from_bounds(0, 1)]
                + [interval_from_bounds(0, 1, False, True)]
            )
        return self._intervals_alternative

    def matrix_alternative_2_systems(self) -> tuple:
        if not hasattr(self, "_matrix_alternative_2_systems1"):
            m_A = self.A.nrows()
            m_B = self.B.nrows()
            length = self.A.ncols()
            self._matrix_alternative_2_systems1 = matrix.block(
                [
                    [self.A.T, self.B.T],
                    [identity_matrix(m_A), zero_matrix(m_A, m_B)],
                    [zero_matrix(m_B, m_A), identity_matrix(m_B)],
                    [
                        matrix(1, -self.b),
                        matrix(1, -self.c),
                    ],  # number of rows is relevant for 0-dim vectors
                ]
            )
        if not hasattr(self, "_matrix_alternative_2_systems2"):
            m_A = self.A.nrows()
            m_B = self.B.nrows()
            length = self.A.ncols()
            self._matrix_alternative_2_systems2 = matrix.block(
                [
                    [self.A.T, self.B.T],
                    [identity_matrix(m_A), zero_matrix(m_A, m_B)],
                    [zero_matrix(m_B, m_A), identity_matrix(m_B)],
                    [
                        matrix(1, -self.b),
                        matrix(1, -self.c),
                    ],  # number of rows is relevant for 0-dim vectors
                    [zero_matrix(1, m_A), ones_matrix(1, m_B)],
                ]
            )
        return self._matrix_alternative_2_systems1, self._matrix_alternative_2_systems2

    def intervals_alternative_2_systems(self) -> tuple:
        if not hasattr(self, "_intervals_alternative_2_systems1"):
            m_A = self.A.nrows()
            m_B = self.B.nrows()
            length = self.A.ncols()
            self._intervals_alternative_2_systems1 = (
                length * [interval_from_bounds(0, 0)]
                + (m_A + m_B) * [interval_from_bounds(0, 1)]
                + [interval_from_bounds(0, 1, False, True)]
            )
        if not hasattr(self, "_intervals_alternative_2_systems2"):
            m_A = self.A.nrows()
            m_B = self.B.nrows()
            length = self.A.ncols()
            self._intervals_alternative_2_systems2 = (
                length * [interval_from_bounds(0, 0)]
                + (m_A + m_B + 1) * [interval_from_bounds(0, 1)]
                + [interval_from_bounds(0, 1, False, True)]
            )
        return self._intervals_alternative_2_systems1, self._intervals_alternative_2_systems2

    def certify(self) -> tuple:
        r"""
        Return whether the system has a solution and certify it.

        OUTPUT:
        A tuple of a boolean and a list of vectors certifying the result.
        """
        return self.certify_using_single_system()

    def certify_using_single_system(self) -> tuple:
        r"""
        Return whether this system has a solution and certify it.

        .. NOTE::

            This implementation uses one system for the alternative.

        OUTPUT:
        A tuple of a boolean and a list of vectors certifying the result.
        """
        m_A, length = self.A.dimensions()
        m_B = self.B.nrows()

        evs = elementary_vectors(self.matrix().T, generator=True)
        evs_alt = elementary_vectors_generator_trailing_nonzero(self.matrix_alternative().T)

        evs_end_reached = False
        evs_alt_end_reached = False
        while True:
            if not evs_end_reached:
                try:
                    v = next(evs)
                    if not exists_orthogonal_vector_inhomogeneous(v, self.b, self.c):
                        return False, [v]
                except StopIteration:
                    evs_end_reached = True

            if not evs_alt_end_reached:
                try:
                    v = next(evs_alt)
                    if not exists_orthogonal_vector_homogeneous(
                        v, [-1], range(length + 1, length + m_A + m_B + 2)
                    ):
                        return True, v
                except StopIteration:
                    evs_alt_end_reached = True

    def certify_using_two_systems(self) -> tuple:
        r"""
        Return whether the system has a solution and certify it.

        .. NOTE::

            This implementation uses two systems for the alternative.

        OUTPUT:
        A tuple of a boolean and a list of vectors certifying the result.
        """
        length = self.A.ncols()
        m_A = self.A.nrows()
        m_B = self.B.nrows()

        M1, M2 = self.matrix_alternative_2_systems()

        evs = elementary_vectors(self.matrix().T, generator=True)
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
                    if not exists_orthogonal_vector_inhomogeneous(v, self.b, self.c):
                        return False, [v]
                except StopIteration:
                    evs_end_reached = True

            if not evs_alt1_end_reached and not v1_found:
                try:
                    v1 = next(evs_alt1)
                    if not exists_orthogonal_vector_homogeneous(
                        v1, [-1], range(length, length + m_A + m_B)
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
                        v2, [-1], range(length, length + m_A + m_B + 1)
                    ):
                        if v1_found:
                            return True, [v1, v2]
                        v2_found = True
                except StopIteration:
                    evs_alt2_end_reached = True
