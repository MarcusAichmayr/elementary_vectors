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
    sage: from vectors_in_intervals.inequality_systems import *
    sage: A = matrix([[1, 2], [0, 1]])
    sage: B = matrix([[2, 3]])
    sage: C = matrix([[-1, 0]])
    sage: S = AlternativesHomogeneous(A, B, C)
    sage: S.one.matrix
    [ 1  2]
    [ 0  1]
    [-----]
    [ 2  3]
    [-----]
    [-1  0]
    sage: S.one.get_intervals()
    [(0, +oo), (0, +oo), [0, +oo), {0}]
    sage: exists_vector(S.one.matrix.T, S.one.get_intervals())
    True

The certify the result, we consider the alternative system
which consists of a single matrix and intervals::

    sage: S.two.matrix
    [ 1  1  0  0]
    [-----------]
    [ 1  0  0  0]
    [ 0  1  0  0]
    [ 0  0  1  0]
    [-----------]
    [ 1  0  2 -1]
    [ 2  1  3  0]
    sage: S.two.get_intervals()
    [(0, +oo), [0, +oo), [0, +oo), [0, +oo), {0}, {0}]
    sage: exists_vector(S.two.matrix.T, S.two.get_intervals())
    False
    sage: exists_vector(S.two.matrix.T, S.two.get_intervals(), certify=True)
    (-1, -1, 0, -3, 0, 1)

There is a single command for certification::

    sage: S.certify()
    (True, (-1, -1, 0, -3, 0, 1), 2)

We consider another example::

    sage: A = matrix([[1, 0], [0, 1]])
    sage: B = matrix([[2, -3]])
    sage: C = matrix([[-1, -1]])
    sage: S = AlternativesHomogeneous(A, B, C)
    sage: S.certify()
    (False, (1, 1, 0, 1), 2)

In some cases, it is faster to randomly generate elementary vectors to certify::

    sage: S.certify(random=True) # random
    (False, (0, -5, -1, -2), 2)

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
    sage: S = AlternativesInhomogeneous(A, B, b, c)
    sage: S.one.matrix
    [1 2]
    [0 1]
    [---]
    [2 3]
    sage: S.one.get_intervals()
    [(-oo, 1], (-oo, -1], (-oo, 2)]
    sage: exists_vector(S.one.matrix.T, S.one.get_intervals())
    True

However, to certify existence of a solution, we need to consider the alternative system.
This system can be described by two matrices and two lists of intervals::

    sage: S.two.matrix
    [ 0  0  1  1]
    [-----------]
    [ 1  0  0  0]
    [ 0  1  0  0]
    [ 0  0  1  0]
    [ 0  0  0  1]
    [-----------]
    [ 1  0  2  0]
    [ 2  1  3  0]
    [-1  1 -2 -1]
    sage: S.two.get_intervals()
    [(0, +oo), [0, +oo), [0, +oo), [0, +oo), [0, +oo), {0}, {0}, {0}]
    sage: exists_vector(S.two.matrix.T, S.two.get_intervals())
    False
    sage: exists_vector(S.two.matrix.T, S.two.get_intervals(), certify=True)
    (1, 3, 0, 4, 0, 0, -1, 1)

The package offers a single function that certifies existence of a solution::

    sage: S.certify()
    (True, (2, 2, 0, 0, 0, 4, -2, 2), 12)

We consider another example::

    sage: A = matrix([[1, 0], [1, 1]])
    sage: B = matrix([[-1, -1]])
    sage: b = vector([1, 0])
    sage: c = vector([0])
    sage: S = AlternativesInhomogeneous(A, B, b, c)
    sage: S.certify()
    (False, (0, 1, 1), 1)
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
from sage.parallel.decorate import parallel
from sage.rings.infinity import Infinity
from sage.sets.real_set import RealSet
from sage.structure.sage_object import SageObject

from elementary_vectors import elementary_vectors
from elementary_vectors.utility import elementary_vector_from_indices, elementary_vector_from_indices_prevent_multiples
from vectors_in_intervals import exists_orthogonal_vector
from .utility import interval_from_bounds, CombinationsIncluding


def elementary_vectors_generator(M, fixed_elements=None, reverse=False, random=False):
    r"""
    Return generator of elementary vectors with nonzero first component.

    INPUT:

    - ``M`` -- a matrix

    - ``fixed_elements`` -- a list of indices where a minor should be used

    - ``reverse`` -- reverse the order of elements

    - ``random`` -- randomly generate elements (repetition can occur)

    .. SEEALSO::

        :func:`elementary_vectors.functions.elementary_vectors`
    """
    try:
        M = M.matrix_from_rows(M.pivot_rows())
    except (ArithmeticError, NotImplementedError):
        warnings.warn("Could not determine rank of matrix. Expect wrong result!")

    rank, length = M.dimensions()
    minors = {}

    if fixed_elements is None:
        combinations = Combinations(length, rank + 1)
    else:
        combinations = CombinationsIncluding(length, rank + 1, fixed_elements)

    if random:
        try:
            while True:
                v = elementary_vector_from_indices(
                    combinations.random_element(),
                    minors,
                    M
                )
                if v:
                    yield v
        except ValueError as exc:
            raise StopIteration("Kernel of matrix is empty. Could not generate elementary vectors.") from exc

    for indices in combinations.generator(reverse=reverse):
        v = elementary_vector_from_indices_prevent_multiples(indices, minors, M)
        if v:
            yield v


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


class LinearInequalitySystem(SageObject):
    r"""
    A class for linear inequality systems given by a matrix and intervals
    """
    __slots__ = "result", "matrix", "intervals", "elementary_vectors"
    def __init__(self, result, _matrix, intervals) -> None:
        self.result = result
        self.matrix = _matrix
        self.intervals = intervals
        self.elementary_vectors = None

    def _repr_(self) -> str:
        return str(self.matrix) + " x in " + str(self.get_intervals())

    def set_elementary_vectors(self, random=False):
        self.elementary_vectors = self.elementary_vectors_generator(random=random)

    def elementary_vectors_generator(self, random=False):
        if random:
            return elementary_vectors_generator(self.matrix.T, random=True)
        return elementary_vectors(self.matrix.T, generator=True)

    def exists_orthogonal_vector(self, v) -> bool:
        return exists_orthogonal_vector(v, self.intervals)


class HomogeneousSystem(LinearInequalitySystem):
    r"""
    A class for homogeneous linear inequality systems

    ``A x > 0``, ``B x >= 0``, ``C x = 0``
    """
    __slots__ = "strict", "nonstrict"
    def __init__(self, result, A, B, C) -> None:
        super().__init__(result, matrix.block([[A], [B], [C]]), None)
        self.strict = range(A.nrows())
        self.nonstrict = range(A.nrows(), A.nrows() + B.nrows())

    def get_intervals(self) -> list:
        return (
            len(self.strict) * [interval_from_bounds(0, Infinity, False, False)]
            + len(self.nonstrict) * [interval_from_bounds(0, Infinity)]
            + (self.matrix.nrows() - len(self.strict) - len(self.nonstrict)) * [interval_from_bounds(0, 0)]
        )

    def elementary_vectors_generator(self, random=False):
        if len(self.strict) == 1:
            return elementary_vectors_generator(
                self.matrix.T,
                self.strict,
                reverse=True, # tends to find certificates faster
                random=random
            )
        return super().elementary_vectors_generator(random=random)

    def exists_orthogonal_vector(self, v) -> bool:
        return exists_orthogonal_vector_homogeneous(v, self.strict, self.nonstrict)


class InhomogeneousSystem(LinearInequalitySystem):
    r"""
    A class for inhomogeneous linear inequality systems

    ``A x <= b``, ``B x <= c``
    """
    def __init__(self, result, A, B, b, c) -> None:
        super().__init__(result, matrix.block([[A], [B]]), None)
        self.b = b
        self.c = c

    def get_intervals(self) -> list:
        return [interval_from_bounds(-Infinity, bi) for bi in self.b] + [interval_from_bounds(-Infinity, ci, False, False) for ci in self.c]

    def exists_orthogonal_vector(self, v) -> bool:
        return exists_orthogonal_vector_inhomogeneous(v, self.b, self.c)


class Alternatives(SageObject):
    r"""
    An abstract class for certifying linear inequality systems
    """
    __slots__ = "one", "two"

    def _repr_(self) -> str:
        return f"Either\n{self.one}\nor\n{self.two}"

    def certify(self, random=False) -> tuple:
        r"""
        Return whether the first alternative has a solution.

        OUTPUT:
        A boolean, a vector certifying the result, and the number of needed iterations.
        """
        systems = [self.one, self.two]
        needed_iterations = 0
        for system in systems:
            system.set_elementary_vectors(random=random)
        while True:
            needed_iterations += 1
            for i, system in enumerate(systems):
                try:
                    v = next(system.elementary_vectors)
                    if not system.exists_orthogonal_vector(v):
                        return system.result, v, needed_iterations
                except StopIteration:
                    systems.pop(i)
                    break


class AlternativesHomogeneous(Alternatives):
    r"""
    A class for certifying homogeneous linear inequality systems
    """
    def __init__(self, A, B, C):
        self.one = HomogeneousSystem(False, A, B, C)
        m_A = A.nrows()
        m_B = B.nrows()
        m_C = C.nrows()
        self.two = HomogeneousSystem(
            True,
            matrix.block([[ones_matrix(1, m_A), zero_matrix(1, m_B), zero_matrix(1, m_C)]]),
            matrix.block([
                [identity_matrix(m_A), zero_matrix(m_A, m_B), zero_matrix(m_A, m_C)],
                [zero_matrix(m_B, m_A), identity_matrix(m_B), zero_matrix(m_B, m_C)]
            ]),
            matrix.block([[A.T, B.T, C.T]])
        )


def inhomogeneous_alternative1(A, B, b, c) -> InhomogeneousSystem:
    return InhomogeneousSystem(False, A, B, b, c)


def inhomogeneous_alternative1_homogenized(A, B, b, c) -> HomogeneousSystem:
    return HomogeneousSystem(
        False,
        matrix.block([
            [zero_matrix(1, A.ncols()), matrix([[-1]])],
            [B, -c.column()]
        ]),
        matrix.block([
            [A, -b.column()],
        ]),
        matrix(0, A.ncols() + 1)
    )


def inhomogeneous_alternative2(A, B, b, c) -> HomogeneousSystem:
    m_A = A.nrows()
    m_B = B.nrows()
    length = A.ncols()
    return HomogeneousSystem(
        True,
        matrix.block([
            [zero_matrix(1, m_A), ones_matrix(1, m_B + 1)]
        ]),
        identity_matrix(m_A + m_B + 1),
        matrix.block([
            [A.T, B.T, zero_matrix(length, 1)],
            [-b.row(), -c.row(), matrix([[-1]])]
        ])
    )


def inhomogeneous_alternative2_system1(A, B, b, c) -> HomogeneousSystem:
    m_A = A.nrows()
    m_B = B.nrows()
    return HomogeneousSystem(
        True,
        matrix.block([[-b.row(), -c.row()]]),
        identity_matrix(m_A + m_B),
        matrix.block([[A.T, B.T]])
    )


def inhomogeneous_alternative2_system2(A, B, b, c) -> HomogeneousSystem:
    m_A = A.nrows()
    m_B = B.nrows()
    return HomogeneousSystem(
        True,
        [zero_matrix(1, m_A), ones_matrix(1, m_B)],
        matrix.block([
            [-b.row(), -c.row()],
            identity_matrix(m_A + m_B),
        ]),
        matrix.block([[A.T, B.T]])
    )


class AlternativesInhomogeneous(Alternatives):
    r"""
    A class for certifying inhomogeneous linear inequality systems
    """
    def __init__(self, A, B, b, c):
        self.one = inhomogeneous_alternative1(A, B, b, c)
        self.two = inhomogeneous_alternative2(A, B, b, c)


class AlternativesInhomogeneousHomogenized(Alternatives):
    r"""
    A class for certifying inhomogeneous linear inequality systems
    """
    def __init__(self, A, B, b, c):
        self.one = inhomogeneous_alternative1_homogenized(A, B, b, c)
        self.two = inhomogeneous_alternative2(A, B, b, c)
