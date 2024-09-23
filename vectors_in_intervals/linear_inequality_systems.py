r"""
Linear inequality systems
=========================

EXAMPLES::

    sage: from vectors_in_intervals import *
    sage: A = matrix([[1, 2], [0, 1]])
    sage: B = matrix([[2, 3]])
    sage: C = matrix([[-1, 0]])
    sage: S = HomogeneousSystem(A, B, C)
    sage: S.solve()
    (0, 1)
    sage: S.certify_nonexistence()
    Traceback (most recent call last):
    ...
    ValueError: A solution exists!
    sage: S.certify()
    (True, (0, 1))
    sage: S.certify(reverse=True)
    (True, (0, 1))
    sage: S.certify_parallel()
    (True, (0, 1))
    sage: S.certify_parallel(reverse=True) # random
    (True, (0, 1))
    sage: S.certify_parallel(random=True) # random
    (True, (0, 1))

We consider another system::

    sage: M = matrix([[1, 0], [0, 1], [1, 1], [0, 1]])
    sage: lower_bounds = [2, 5, 0, -oo]
    sage: upper_bounds = [5, oo, 8, 5]
    sage: lower_bounds_closed = [True, True, False, False]
    sage: upper_bounds_closed = [False, False, False, True]
    sage: I = intervals_from_bounds(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)
    sage: S = LinearInequalitySystem(M, I)
    sage: S.certify()
    (True, (2, 5))
    sage: S.certify(reverse=True)
    (True, (5/2, 5))
    sage: # S.certify_parallel() # TODO SignalError: Segmentation Fault
    (True, (2, 5))

We consider yet another system::

    sage: A = matrix([[1, 0], [1, 1]])
    sage: B = matrix([[-1, -1]])
    sage: b = vector([1, 0])
    sage: c = vector([0])
    sage: S = InhomogeneousSystem(A, B, b, c)
    sage: S.certify()
    (False, (0, 1, 1))
    sage: S.certify_parallel()
    (False, (0, 1, 1))
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

import concurrent.futures

from collections.abc import Generator
from sage.matrix.constructor import matrix, zero_matrix
from sage.modules.free_module_element import vector
from sage.rings.infinity import Infinity
from sage.structure.sage_object import SageObject

from elementary_vectors.functions import ElementaryVectors
from sign_vectors import sign_vector
from vectors_in_intervals import exists_orthogonal_vector
from .construction import vector_between_sign_vectors
from .utility import interval_from_bounds, CombinationsIncluding


class LinearInequalitySystem(SageObject):
    r"""
    A class for linear inequality systems given by a matrix and intervals
    """
    __slots__ = "result", "matrix", "intervals", "evs", "elementary_vectors", "solvable"

    def __init__(self, _matrix, intervals, result=None) -> None:
        self.matrix = _matrix
        self.intervals = intervals
        self.result = result
        self.evs = ElementaryVectors(self.matrix.T)
        self.solvable = None

    def _repr_(self) -> str:
        return str(self.matrix) + " x in " + str(self.get_intervals())

    def get_intervals(self):
        r"""Return the corresponding intervals."""
        return self.intervals

    def exists_orthogonal_vector(self, v) -> bool:
        r"""Check if an orthogonal vector exists in the intervals."""
        return exists_orthogonal_vector(v, self.intervals)

    def elementary_vectors_generator(self, reverse: bool = False, random: bool = False) -> Generator:
        r"""Return a generator of elementary vectors."""
        if random:
            while True:
                yield self.evs.random_element()
        else:
            yield from self.evs.generator(reverse=reverse)

    def to_homogeneous(self):
        r"""Return the equivalent homogeneous system."""
        return HomogeneousSystem(*homogeneous_from_general(self.matrix, self.intervals), result=self.result)

    def solve(self, reverse: bool = False):
        r"""
        Compute a solution for this linear inequality system.

        If no solution exists, a ``ValueError`` is raised.
        """
        homogeneous = self.to_homogeneous()
        solution = homogeneous.solve(reverse=reverse)
        return solution[:-1] / solution[-1]

    def certify_nonexistence(self, reverse: bool = False, random: bool = False):
        r"""
        Certifies nonexistence if no solution exists.

        Otherwise, a ``ValueError`` is raised.

        .. NOTE::

            If a solution exists and ``random`` is set to true, this method will never finish.
        """
        for v in self.elementary_vectors_generator(reverse=reverse, random=random):
            if self.solvable:
                break
            if not self.exists_orthogonal_vector(v):
                self.solvable = False
                return v
        self.solvable = True
        raise ValueError("A solution exists!")

    def certify(self, reverse: bool = False):
        r"""Return a boolean and a certificate for solvability."""
        try:
            return False, self.certify_nonexistence(reverse=reverse)
        except ValueError:
            return True, self.solve(reverse=reverse)

    def certify_parallel(self, reverse: bool = False, random: bool = False):
        r"""
        Return a boolean and a certificate for solvability.

        Attempts to find a solution and certify nonexistence in parallel.
        """
        with concurrent.futures.ThreadPoolExecutor() as executor:
            done, not_done = concurrent.futures.wait(
                [
                    executor.submit(lambda: (False, self.certify_nonexistence(reverse=not reverse, random=random))),
                    executor.submit(lambda: (True, self.solve(reverse=reverse)))
                ],
                return_when=concurrent.futures.FIRST_COMPLETED
            )

            try:
                result = done.pop().result()
                for task in not_done:
                    task.cancel()

                return result
            except ValueError:
                for task in not_done:
                    return task.result()


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
        return not (
            any(v[k] != 0 for k in self.strict)
            and (
                (
                    all(v[k] >= 0 for k in self.strict)
                    and all(v[k] >= 0 for k in self.nonstrict)
                )
                or (
                    all(v[k] <= 0 for k in self.strict)
                    and all(v[k] <= 0 for k in self.nonstrict)
                )
            )
        )

    def to_homogeneous(self):
        return self

    def solve(self, reverse: bool = False):
        r"""
        Compute a solution if existent.

        This approach sums up positive elementary vectors in the row space.
        """
        try:
            result = self.matrix.solve_right(
                vector_between_sign_vectors(
                    self.evs.generator(kernel=False, reverse=reverse),
                    sign_vector(
                        len(self.strict) * [1] + (self.matrix.nrows() - len(self.strict)) * [0]
                    ),
                    sign_vector(
                        (len(self.strict) + len(self.nonstrict)) * [1]
                        + (self.matrix.nrows() - len(self.strict) - len(self.nonstrict)) * [0]
                    )
                )
            )
            self.solvable = True
            return result
        except ValueError as exc:
            self.solvable = False
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
        len_b = len(self.b)
        len_c = len(self.c)

        def condition(v):
            if all(vk >= 0 for vk in v):
                scalarproduct = sum(v[k] * self.b[k] for k in range(len_b)) + sum(
                    v[k + len_b] * self.c[k] for k in range(len_c)
                )
                if scalarproduct < 0:
                    return True
                if scalarproduct <= 0 and any(v[k] for k in range(len_b, len_b + len_c)):
                    return True
            return False

        return not condition(v) and not condition(-v)

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
