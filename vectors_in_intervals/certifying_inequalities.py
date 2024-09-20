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
    sage: exists_vector(S.one.matrix, S.one.get_intervals())
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
    sage: exists_vector(S.two.matrix, S.two.get_intervals())
    False
    sage: exists_vector(S.two.matrix, S.two.get_intervals(), certify=True)
    (-1, -1, 0, -3, 0, 1)

There is a single command for certification::

    sage: S.certify()
    (True, (-1, -1, 0, -3, 0, 1), 2)

We can also use elementary vectors to construct a solution::

    sage: S.one.solve()
    (0, 1)

That way, we can certify solvability using only one alternative::

    sage: S.one.certify()
    (True, (0, 1))

We consider another example::

    sage: A = matrix([[1, 0], [0, 1]])
    sage: B = matrix([[2, -3]])
    sage: C = matrix([[-1, -1]])
    sage: S = AlternativesHomogeneous(A, B, C)
    sage: S.certify()
    (False, (0, -5, -1, -2), 1)
    sage: S.one.certify()
    (False, (-2, 3, 1, 0))

In some cases, it is faster to randomly generate elementary vectors to certify::

    sage: S.certify(random=True) # random
    (False, (0, -5, -1, -2), 2)

Parallel computation is also supported::

    sage: S.certify(number_parallel=4) # random
    (False, (0, -5, -1, -2), 2)

A solution for the second system is::

    sage: S.two.solve()
    (1, 1, 0, 1)

Inhomogeneous systems
~~~~~~~~~~~~~~~~~~~~~

We deal with linear inequality systems given in the form
:math:`A x \leq b, B x < c` with matrices :math:`A`, :math:`B` and vectors :math:`b`, :math:`c`.
We are interested whether the system  has a solution :math:`x`.
This can be checked using :func:`vectors_in_intervals.existence.exists_vector`.
If no solution exists, this function can be used to certify the result.

To demonstrate this, consider the following example::

    sage: from vectors_in_intervals import *
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
    sage: exists_vector(S.one.matrix, S.one.get_intervals())
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
    sage: exists_vector(S.two.matrix, S.two.get_intervals())
    False
    sage: exists_vector(S.two.matrix, S.two.get_intervals(), certify=True)
    (1, 3, 0, 4, 0, 0, -1, 1)

The package offers a single function that certifies existence of a solution::

    sage: S.certify()
    (True, (2, 2, 0, 0, 0, 4, -2, 2), 13)

The alternatives for the inhomogeneous case can be formulated using two systems.
The resulting systems yield different certificates::

    sage: S = AlternativesInhomogeneous(A, B, b, c, two_double_system=True)
    sage: S.certify()
    (True, [(-2, 0, -1, 0, 0, 1, 0), (-2, -1, 0, 0, -5, 2)], 5)

A solution is::

    sage: S.one.solve()
    (2, -1)

We consider another example::

    sage: A = matrix([[1, 0], [1, 1]])
    sage: B = matrix([[-1, -1]])
    sage: b = vector([1, 0])
    sage: c = vector([0])
    sage: S = AlternativesInhomogeneous(A, B, b, c)
    sage: S.certify()
    (False, (0, 1, 1), 1)

We can homogenized the first alternative yielding a different system::

    sage: S = AlternativesInhomogeneous(A, B, b, c, one_homogenized=True)
    sage: S.certify()
    (False, (0, 1, 0, 1), 1)

A solution is::

    sage: S.two.solve()
    (0, 1, 1, 0)

General systems
~~~~~~~~~~~~~~~

By translating systems of the form ``M x in I`` into inhomogeneous systems,
we can certify general systems::

    sage: from vectors_in_intervals import *
    sage: M = matrix([[1, 0], [0, 1], [1, 1], [0, 1]])
    sage: lower_bounds = [2, 5, 0, -oo]
    sage: upper_bounds = [5, oo, 8, 5]
    sage: lower_bounds_closed = [True, True, False, False]
    sage: upper_bounds_closed = [False, False, False, True]
    sage: I = intervals_from_bounds(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)
    sage: exists_vector(M, I) # no certificate
    True
    sage: S = AlternativesGeneral(M, I)
    sage: S.one.matrix
    [1 0]
    [0 1]
    [1 1]
    [0 1]
    sage: S.one.get_intervals()
    [[2, 5), [5, +oo), (0, 8), (-oo, 5]]
    sage: S.two.matrix
    [ 0  0  0  1  1  1  1]
    [--------------------]
    [ 1  0  0  0  0  0  0]
    [ 0  1  0  0  0  0  0]
    [ 0  0  1  0  0  0  0]
    [ 0  0  0  1  0  0  0]
    [ 0  0  0  0  1  0  0]
    [ 0  0  0  0  0  1  0]
    [ 0  0  0  0  0  0  1]
    [--------------------]
    [-1  0  0  1 -1  1  0]
    [ 0 -1  1  0 -1  1  0]
    [ 2  5 -5 -5  0 -8 -1]
    sage: S.two.get_intervals()
    [(0, +oo), [0, +oo), [0, +oo), [0, +oo), [0, +oo), [0, +oo), [0, +oo), [0, +oo), {0}, {0}, {0}]
    sage: S.certify()
    (True, (1, 0, 0, 0, 2, 6, 0, 0, 2, 5, 1), 3)

We compute a solution using elementary vectors::

    sage: S.one.solve()
    (2, 5)
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

from copy import copy
from sage.matrix.constructor import matrix, zero_matrix, identity_matrix, ones_matrix
from sage.modules.free_module_element import vector
from sage.parallel.decorate import parallel
from sage.rings.infinity import Infinity
from sage.structure.sage_object import SageObject

from elementary_vectors.functions import ElementaryVectors
from sign_vectors import sign_vector
from vectors_in_intervals import exists_orthogonal_vector
from .construction import vector_between_sign_vectors
from .utility import interval_from_bounds, CombinationsIncluding


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

        sage: from elementary_vectors import *
        sage: from vectors_in_intervals import *
        sage: from vectors_in_intervals.certifying_inequalities import *
        sage: A = matrix([[1, 2], [0, 1], [2, 0]])
        sage: B = matrix([[2, 3], [0, 2]])
        sage: b = vector([1, -1, 2])
        sage: c = vector([2, 3])
        sage: I = [interval_from_bounds(-Infinity, bi) for bi in b] + [interval_from_bounds(-Infinity, ci, False, False) for ci in c]
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

        sage: from elementary_vectors import *
        sage: from vectors_in_intervals import *
        sage: from vectors_in_intervals.certifying_inequalities import *
        sage: A = matrix([[1, 2], [0, 1]])
        sage: B = matrix([[2, 3], [0, -1]])
        sage: C = matrix([[-1, 0], [1, 1]])
        sage: M = matrix.block([[A.T, B.T, C.T]])
        sage: I = (
        ....:     2 * [interval_from_bounds(0, Infinity, False, False)]
        ....:     + 2 * [interval_from_bounds(0, Infinity)]
        ....:     + 2 * [interval_from_bounds(0, 0)]
        ....: )
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
    __slots__ = "result", "matrix", "intervals", "evs", "elementary_vectors"

    def __init__(self, _matrix, intervals, result=None) -> None:
        self.matrix = _matrix
        self.intervals = intervals
        self.result = result
        self.evs = ElementaryVectors(self.matrix.T)

    def _repr_(self) -> str:
        return str(self.matrix) + " x in " + str(self.get_intervals())

    def get_intervals(self):
        return self.intervals

    def exists_orthogonal_vector(self, v) -> bool:
        return exists_orthogonal_vector(v, self.intervals)

    def to_homogeneous(self):
        return HomogeneousSystem(*homogeneous_from_general(self.matrix, self.intervals), result=self.result)

    def solve(self):
        homogeneous = self.to_homogeneous()
        solution = homogeneous.solve()
        return solution[:-1] / solution[-1]

    def certify_nonexistence(self, reverse=False, random=False):
        if random:
            def elementary_vectors_generator():
                while True:
                    yield self.evs.random_element()
        else:
            def elementary_vectors_generator():
                return self.evs.generator(reverse=reverse)

        for v in elementary_vectors_generator():
            if self.exists_orthogonal_vector(v):
                return v
        raise ValueError("A solution exists!")

    def certify(self):
        r"""Return a boolean and a certificate for solvability."""
        try:
            return True, self.solve()
        except ValueError:
            return False, self.certify_nonexistence()


class HomogeneousSystem(LinearInequalitySystem):
    r"""
    A class for homogeneous linear inequality systems

    ``A x > 0``, ``B x >= 0``, ``C x = 0``
    """
    __slots__ = "strict", "nonstrict"

    def __init__(self, A, B, C, result=None) -> None:
        super().__init__(matrix.block([[A], [B], [C]]), None, result=result)
        self.strict = range(A.nrows())
        self.nonstrict = range(A.nrows(), A.nrows() + B.nrows())

        if len(self.strict) == 1:
            self.evs.combinations = CombinationsIncluding(self.evs.length, self.evs.rank + 1, self.strict)

    def get_intervals(self) -> list:
        self.intervals = [
            interval_from_bounds(0, Infinity, False, False)
            if i in self.strict else
            (
                interval_from_bounds(0, Infinity)
                if i in self.nonstrict else
                interval_from_bounds(0, 0)
            )
            for i in range(self.matrix.nrows())
        ]
        return self.intervals

    def exists_orthogonal_vector(self, v) -> bool:
        return exists_orthogonal_vector_homogeneous(v, self.strict, self.nonstrict)

    def to_homogeneous(self):
        return self

    def solve(self):
        r"""
        Compute a solution if existent.

        This approach sums up positive elementary vectors in the row space.
        """
        try:
            return self.matrix.solve_right(
                vector_between_sign_vectors(
                    self.matrix.T,
                    sign_vector(
                        len(self.strict) * [1] + (self.matrix.nrows() - len(self.strict)) * [0]
                    ),
                    sign_vector(
                        (len(self.strict) + len(self.nonstrict)) * [1]
                        + (self.matrix.nrows() - len(self.strict) - len(self.nonstrict)) * [0]
                    )
                )
            )
        except ValueError as exc:
            raise ValueError("No solution exists!") from exc


class InhomogeneousSystem(LinearInequalitySystem):
    r"""
    A class for inhomogeneous linear inequality systems

    ``A x <= b``, ``B x <= c``
    """
    def __init__(self, A, B, b, c, result=None) -> None:
        super().__init__(matrix.block([[A], [B]]), None, result=result)
        self.b = b
        self.c = c

    def get_intervals(self) -> list:
        self.intervals = [interval_from_bounds(-Infinity, bi) for bi in self.b] + [interval_from_bounds(-Infinity, ci, False, False) for ci in self.c]
        return self.intervals

    def exists_orthogonal_vector(self, v) -> bool:
        return exists_orthogonal_vector_inhomogeneous(v, self.b, self.c)

    def to_homogeneous(self):
        return HomogeneousSystem(
            *homogeneous_from_inhomogeneous(
                self.matrix.submatrix(0, 0, len(self.b)),
                self.matrix.submatrix(len(self.b), 0, len(self.c)),
                self.b,
                self.c
            ),
            result=self.result
        )


class Alternatives(SageObject):
    r"""
    An abstract class for certifying linear inequality systems.
    """
    __slots__ = "one", "two"

    def _repr_(self) -> str:
        return f"Either\n{self.one}\nor\n{self.two}"

    def certify(self, reverse=True, random=False, number_parallel=1) -> tuple:
        r"""
        Certify whether the first alternative has a solution.

        - ``reverse`` -- reverses the order of elementary vectors

        - ``random`` -- randomizes computed elementary vectors (repetition possible)

        - ``number_parallel`` -- number of parallel computations

        OUTPUT:
        A boolean, a vector certifying the result, and the number of needed iterations.
        """
        systems = [self.one, self.two]
        needed_iterations = 0
        for system in systems:
            system.elementary_vectors = system.evs.generator(reverse=reverse)

        if number_parallel <= 1:
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

        results = []

        @parallel("reference")
        def check_random_ev(system):
            if results:
                return
            try:
                v = next(system.elementary_vectors)
            except StopIteration:
                return
            if not system.exists_orthogonal_vector(v):
                results.append((system.result, v, -1))

        while not results:
            for _ in check_random_ev([*systems] * number_parallel):
                pass

        return results[0]


class AlternativesHomogeneous(Alternatives):
    r"""
    A class for certifying homogeneous linear inequality systems.

    ``A x > 0``, ``B x >= 0``, ``C x = 0``
    """
    def __init__(self, A, B, C):
        self.one = HomogeneousSystem(A, B, C, result=False)
        m_A = A.nrows()
        m_B = B.nrows()
        m_C = C.nrows()
        self.two = HomogeneousSystem(
            matrix.block([[ones_matrix(1, m_A), zero_matrix(1, m_B), zero_matrix(1, m_C)]]),
            matrix.block([
                [identity_matrix(m_A), zero_matrix(m_A, m_B), zero_matrix(m_A, m_C)],
                [zero_matrix(m_B, m_A), identity_matrix(m_B), zero_matrix(m_B, m_C)]
            ]),
            matrix.block([[A.T, B.T, C.T]]),
            result=True
        )


def inhomogeneous_alternative1(A, B, b, c) -> InhomogeneousSystem:
    r"""
    A standard inhomogeneous linear inequality system.
    
    ``A x <= b``, ``B x < c``
    """
    return InhomogeneousSystem(A, B, b, c, result=False)


def inhomogeneous_alternative1_homogenized(A, B, b, c) -> HomogeneousSystem:
    r"""
    Homogenization of a standard inhomogeneous linear inequality system.

    The system has this form::

        [A -b]
        [B -c] x in [0, oo)^p x (0, oo)^(q + 1)
        [0 -1]
    """
    return HomogeneousSystem(
        matrix.block([
            [zero_matrix(1, A.ncols()), matrix([[-1]])],
            [B, -c.column()]
        ]),
        matrix.block([
            [A, -b.column()],
        ]),
        matrix(0, A.ncols() + 1),
        result=False
    )


def inhomogeneous_alternative2(A, B, b, c) -> HomogeneousSystem:
    r"""
    Alternative of a standard inhomogeneous linear inequality system.
    """
    m_A = A.nrows()
    m_B = B.nrows()
    length = A.ncols()
    return HomogeneousSystem(
        matrix.block([
            [zero_matrix(1, m_A), ones_matrix(1, m_B + 1)]
        ]),
        identity_matrix(m_A + m_B + 1),
        matrix.block([
            [A.T, B.T, zero_matrix(length, 1)],
            [-b.row(), -c.row(), matrix([[-1]])]
        ]),
        result=True
    )


def inhomogeneous_alternative2_system1(A, B, b, c) -> HomogeneousSystem:
    r"""
    Alternative of a standard inhomogeneous linear inequality system given by two systems.
    """
    m_A = A.nrows()
    m_B = B.nrows()
    return HomogeneousSystem(
        matrix.block([[-b.row(), -c.row()]]),
        identity_matrix(m_A + m_B),
        matrix.block([[A.T, B.T]]),
        result=True
    )


def inhomogeneous_alternative2_system2(A, B, b, c) -> HomogeneousSystem:
    r"""
    Alternative of a standard inhomogeneous linear inequality system given by two systems.
    """
    m_A = A.nrows()
    m_B = B.nrows()
    return HomogeneousSystem(
        matrix.block([[zero_matrix(1, m_A), ones_matrix(1, m_B)]]),
        matrix.block([
            [matrix.block([[-b.row(), -c.row()]])],
            [identity_matrix(m_A + m_B)],
        ]),
        matrix.block([[A.T, B.T]]),
        result=True
    )


def inhomogeneous_from_general(M, I) -> tuple:
    r"""
    Translate a general system into an inhomogeneous system.

    INPUT:

    - ``M`` -- a matrix with m rows

    - ``I`` -- a list of m intervals

    - To use a homogenized representation of the first alternative, pass ``one_homogenized=True``.

    - To use two systems for the second alternative instead of a homogenized system, pass ``two_double_system=True``.

    OUTPUT:
    Matrices ``A``, ``B`` and vectors ``b``, ``c`` describing the equivalent system ``A x <= b``, ``B x < c``.
    """
    A_list = []
    B_list = []
    b_list = []
    c_list = []

    for line, interval in zip(M, I):
        if interval.inf() != -Infinity:
            if interval.inf() in interval:
                A_list.append(-line)
                b_list.append(-interval.inf())
            else:
                B_list.append(-line)
                c_list.append(-interval.inf())
        if interval.sup() != Infinity:
            if interval.sup() in interval:
                A_list.append(line)
                b_list.append(interval.sup())
            else:
                B_list.append(line)
                c_list.append(interval.sup())

    return (
        matrix(len(A_list), M.ncols(), A_list),
        matrix(len(B_list), M.ncols(), B_list),
        vector(b_list),
        vector(c_list)
    )


def homogeneous_from_general(M, I) -> tuple:
    r"""
    Convert a general system to a homogeneous system.

    INPUT:

    - ``M`` -- a matrix with m rows

    - ``I`` -- a list of m intervals

    OUTPUT:
    Matrices ``A``, ``B``, ``C`` describing the homogeneous system.

    EXAMPLE::

        sage: from vectors_in_intervals import *
        sage: M = matrix([[1, 0], [0, 1], [1, 1]])
        sage: lower_bounds = [2, 5, 0]
        sage: upper_bounds = [5, oo, 0]
        sage: lower_bounds_closed = [True, True, True]
        sage: upper_bounds_closed = [False, False, True]
        sage: I = intervals_from_bounds(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)
        sage: homogeneous_from_general(M, I)
        (
        [ 1  0 -5]  [-1  0  2]         
        [ 0  0 -1], [ 0 -1  5], [1 1 0]
        )
    """
    A_list = []
    B_list = []
    C_list = []

    length = M.ncols()

    for line, interval in zip(M, I):
        if interval.inf() == interval.sup():
            C_list.append(list(line) + [-interval.inf()])
            continue
        if interval.inf() != -Infinity:
            if interval.inf() in interval:
                B_list.append(list(-line) + [interval.inf()])
            else:
                A_list.append(list(-line) + [interval.inf()])
        if interval.sup() != Infinity:
            if interval.sup() in interval:
                B_list.append(list(line) + [-interval.sup()])
            else:
                A_list.append(list(line) + [-interval.sup()])

    A_list.append([0] * length + [-1])

    return (
        matrix(len(A_list), length + 1, A_list),
        matrix(len(B_list), length + 1, B_list),
        matrix(len(C_list), length + 1, C_list)
    )


def homogeneous_from_inhomogeneous(A, B, b, c) -> tuple:
    r"""
    Convert an inhomogeneous system to a homogeneous system.

    OUTPUT:
    Matrices ``A``, ``B``, ``C`` describing the homogeneous system.
    """
    return (
        matrix.block([[B, matrix(len(c), 1, -c)], [zero_matrix(1, A.ncols()), matrix([[-1]])]]),
        matrix.block([[A, matrix(len(b), 1, -b)]]),
        matrix(0, A.ncols() + 1)
    )


class AlternativesInhomogeneous(Alternatives):
    r"""
    A class for certifying inhomogeneous linear inequality systems.

    ``A x <= b``, ``B x < c``

    INPUT:

    - To use a homogenized representation of the first alternative, pass ``one_homogenized=True``.

    - To use two systems for the second alternative instead of a homogenized system, pass ``two_double_system=True``.

    TESTS::

        sage: from vectors_in_intervals import *
        sage: A = matrix([[1]])
        sage: B = matrix(0, 1)
        sage: b = vector([1])
        sage: c = vector([])
        sage: S = AlternativesInhomogeneous(A, B, b, c)
        sage: S.certify()
        (True, (-1, 0, 0, -1, -1), 2)
        sage: S.certify(random=True)[0]
        True
    """
    def __init__(self, A, B, b, c, one_homogenized=False, two_double_system=False):
        if one_homogenized:
            self.one = inhomogeneous_alternative1_homogenized(A, B, b, c)
        else:
            self.one = inhomogeneous_alternative1(A, B, b, c)
        
        if two_double_system:
            self.two = inhomogeneous_alternative2_system1(A, B, b, c)
            self.three = inhomogeneous_alternative2_system2(A, B, b, c)
        else:
            self.two = inhomogeneous_alternative2(A, B, b, c)

    def certify(self, reverse=True, random=False, number_parallel=1) -> tuple:
        r"""
        Certify whether the first alternative has a solution.

        OUTPUT:
        A boolean, vectors certifying the result, and the number of needed iterations.
        """
        if not hasattr(self, "three"):
            return super().certify(reverse=reverse, random=random, number_parallel=number_parallel)

        if number_parallel > 1:
            raise NotImplementedError("Parallel certification for 2-system alternative not implemented.")

        systems = {0: self.one, 1: self.two, 2: self.three}
        needed_iterations = 0
        for system in systems.values():
            system.elementary_vectors = system.evs.generator(reverse=reverse)

        certificates = {}
        while True:
            needed_iterations += 1
            for i, system in copy(systems).items():
                try:
                    v = next(system.elementary_vectors)
                    if not system.exists_orthogonal_vector(v):
                        if i == 0:
                            return system.result, v, needed_iterations
                        certificates[i] = v
                        if len(certificates) == 2:
                            return system.result, list(certificates.values()), needed_iterations
                        systems.pop(i)
                except StopIteration:
                    systems.pop(i)
                    break


class AlternativesGeneral(Alternatives):
    r"""Alternatives for a general system ``M x in I``."""
    def __init__(self, M, I):
        self.one = LinearInequalitySystem(M, I, result=False)
        self.two = inhomogeneous_alternative2(*inhomogeneous_from_general(M, I))
