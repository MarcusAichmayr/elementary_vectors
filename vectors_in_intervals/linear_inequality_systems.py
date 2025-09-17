r"""
Linear inequality systems
=========================

EXAMPLES::

    sage: from vectors_in_intervals import *
    sage: A = matrix([[1, 2], [0, 1]])
    sage: B = matrix([[2, 3]])
    sage: C = matrix([[-1, 0]])
    sage: S = HomogeneousSystem(A, B, C)
    sage: S.intervals
    [(0, +oo), (0, +oo), [0, +oo), {0}]
    sage: S.solve()
    (0, 1)
    sage: S.certify_existence()
    (2, 1, 3, 0)
    sage: S.certify_nonexistence()
    Traceback (most recent call last):
    ...
    ValueError: A solution exists!
    sage: S.certify()
    (True, (2, 1, 3, 0))
    sage: S.certify(reverse=True)
    (True, (2, 1, 3, 0))
    sage: S.certify_parallel()
    (True, (2, 1, 3, 0))
    sage: S.certify_parallel(reverse=True)
    (True, (2, 1, 3, 0))
    sage: S.certify_parallel(random=True)
    (True, (2, 1, 3, 0))

We consider another system::

    sage: M = matrix([[1, 0], [0, 1], [1, 1], [0, 1]])
    sage: lower_bounds = [2, 5, 0, -oo]
    sage: upper_bounds = [5, oo, 8, 5]
    sage: lower_bounds_closed = [True, True, False, False]
    sage: upper_bounds_closed = [False, False, False, True]
    sage: I = Intervals.from_bounds(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)
    sage: S = LinearInequalitySystem(M, I)
    sage: S.solve()
    (5/2, 5)
    sage: S.certify_existence()
    (5, 15, 1, 2, 1, 0, 0)
    sage: S.certify()
    (True, (5, 15, 1, 2, 1, 0, 0))
    sage: S.certify(reverse=True)
    (True, (3, 7, 1, 1, 0, 0, 0))
    sage: # S.certify_parallel() # TODO SignalError: Segmentation Fault
    (True, (3, 7, 1, 1, 0, 0, 0))

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
    sage: S.certify_parallel(random=True)
    (False, (0, 1, 1))

In the case of homogeneous systems, we can use cocircuits to certify::

    sage: A = matrix([[1, 2], [0, 1]])
    sage: B = matrix([[2, 3]])
    sage: C = matrix([[-1, 0]])
    sage: S = HomogeneousSystemCocircuits(A, B, C) # TODO: not implemented
    sage: S.solve() # TODO: not implemented
    Traceback (most recent call last):
    ...
    ValueError: Can't solve using cocircuits!
    sage: S.certify_existence() # TODO: not implemented
    (+++0)
    sage: # S.certify() # TODO

Now, we consider the example::

    sage: A = matrix([[1, 0], [0, 1]])
    sage: B = matrix([[2, -3]])
    sage: C = matrix([[-1, -1]])
    sage: S = HomogeneousSystemCocircuits(A, B, C) # TODO: not implemented
    sage: S.certify_nonexistence() # TODO: not implemented
    (++0+)
"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from __future__ import annotations

import concurrent.futures

from collections.abc import Generator
from sage.matrix.constructor import Matrix, zero_matrix
from sage.modules.free_module_element import vector, zero_vector
from sage.rings.infinity import Infinity
from sage.structure.sage_object import SageObject

from sign_vectors import SignVector, sign_vector, zero_sign_vector
from sign_vectors.oriented_matroids import OrientedMatroid
from vectors_in_intervals import exists_orthogonal_vector, Intervals, Interval
from elementary_vectors.functions import ElementaryVectors
from .utility import CombinationsIncluding, solve_without_division


class LinearInequalitySystem(SageObject):
    r"""
    A class for linear inequality systems given by a matrix and intervals
    """
    __slots__ = "result", "matrix", "_intervals", "_evs", "elementary_vectors", "_solvable"

    def __init__(self, matrix: Matrix, intervals: Intervals = None, result: bool = None) -> None:
        self.matrix = matrix
        self._intervals = intervals
        self.result = result
        self._evs = ElementaryVectors(self.matrix.T)
        self._solvable = None

    def _repr_(self) -> str:
        return str(self.matrix) + " x in " + str(self.intervals)

    @property
    def intervals(self) -> Intervals:
        r"""Return the corresponding intervals."""
        if self._intervals is None:
            self._intervals = self._compute_intervals()
        return self._intervals

    def _compute_intervals(self) -> Intervals:
        raise NotImplementedError("This method should be implemented in subclasses.")

    def exists_orthogonal_vector(self, v) -> bool:
        r"""Check if an orthogonal vector exists in the intervals."""
        return exists_orthogonal_vector(v, self._intervals)

    def candidate_generator(self, dual: bool = True, reverse: bool = False, random: bool = False) -> Generator:
        r"""Return a generator of elementary vectors."""
        if random:
            while True:
                yield self._evs.random_element(dual=dual)
        else:
            yield from self._evs.generator(dual=dual, reverse=reverse)

    def to_homogeneous(self) -> HomogeneousSystem:
        r"""Return the equivalent homogeneous system."""
        return homogeneous_from_general(self)

    def certify_nonexistence(self, reverse: bool = False, random: bool = False):
        r"""
        Certify nonexistence of solutions.

        Otherwise, a ``ValueError`` is raised.

        .. NOTE::

            If a solution exists and ``random`` is set to true, this method will never finish.
        """
        for v in self.candidate_generator(reverse=reverse, random=random):
            if self._solvable:
                break
            if not self.exists_orthogonal_vector(v):
                self._solvable = False
                return v
        self._solvable = True
        raise ValueError("A solution exists!")

    def certify_existence(self, reverse: bool = False, random: bool = False):
        r"""
        Certify existence of a solution if one exists.

        Otherwise, a ``ValueError`` is raised.

        .. NOTE::

            If a solution exists and ``random`` is set to true, this method will never finish.
        """
        return self.to_homogeneous().certify_existence(reverse=reverse, random=random)

    def solve(self, reverse: bool = False, random: bool = False):
        r"""
        Compute a solution for this linear inequality system.

        If no solution exists, a ``ValueError`` is raised.
        """
        solution = self.to_homogeneous().solve(reverse=reverse, random=random)
        return solution[:-1] / solution[-1]

    def certify(self, reverse: bool = False) -> tuple:
        r"""Return a boolean and a certificate for solvability."""
        try:
            return False, self.certify_nonexistence(reverse=reverse)
        except ValueError:
            return True, self.certify_existence(reverse=reverse)

    def certify_parallel(self, reverse: bool = False, random: bool = False) -> tuple:
        r"""
        Return a boolean and a certificate for solvability.

        Attempts to find a solution and certify nonexistence in parallel.
        """
        with concurrent.futures.ThreadPoolExecutor() as executor:
            done, not_done = concurrent.futures.wait(
                [
                    executor.submit(lambda: (False, self.certify_nonexistence(reverse=reverse, random=random))),
                    executor.submit(lambda: (True, self.certify_existence(reverse=reverse, random=random)))
                ],
                return_when=concurrent.futures.FIRST_COMPLETED
            )

            for task in done:
                try:
                    result = task.result()
                    for task in not_done:
                        task.cancel()

                    return result
                except ValueError:
                    pass
            return not_done.pop().result()


class InhomogeneousSystem(LinearInequalitySystem):
    r"""
    A class for inhomogeneous linear inequality systems

    ``A x <= b``, ``B x <= c``
    """
    def __init__(self, A: Matrix, B: Matrix, b: vector, c: vector, result: bool = None) -> None:
        super().__init__(Matrix.block([[A], [B]]), None, result=result)
        self.A = A
        self.B = B
        self.b = b
        self.c = c

    def _compute_intervals(self) -> Intervals:
        return [Interval(-Infinity, bi) for bi in self.b] + [Interval(-Infinity, ci, False, False) for ci in self.c]

    def exists_orthogonal_vector(self, v) -> bool:
        len_b = len(self.b)
        len_c = len(self.c)

        def condition(v) -> bool:
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

    def to_homogeneous(self) -> HomogeneousSystem:
        return homogeneous_from_inhomogeneous(self)


class HomogeneousSystem(LinearInequalitySystem):
    r"""
    A class for homogeneous linear inequality systems

    ``A x > 0``, ``B x >= 0``, ``C x = 0``

    TESTS::

        sage: from vectors_in_intervals import *
        sage: A = matrix([[0, 1], [0, 1], [0, 1]])
        sage: B = zero_matrix(0, 2)
        sage: C = matrix([[1, 1], [0, 0]])
        sage: S = HomogeneousSystem(A, B, C)
        sage: S.certify_existence()
        (1, 1, 1, 0, 0)
    """
    __slots__ = "_positive", "_nonnegative", "_zero"

    def __init__(self, A: Matrix, B: Matrix, C: Matrix, result: bool = None) -> None:
        super().__init__(Matrix.block([[A], [B], [C]]), None, result=result)
        self._positive = range(A.nrows())
        self._nonnegative = range(A.nrows() + B.nrows())
        self._zero = range(A.nrows() + B.nrows(), self.matrix.nrows())

        # self._evs._set_combinations_row_space(Combinations(range(A.nrows() + B.nrows()), self._evs.length - self._evs.rank + 1))

        if len(self._positive) == 1:
            self._evs._set_combinations_kernel(CombinationsIncluding(self._evs.length, self._evs.rank + 1, self._positive))

    def _compute_intervals(self) -> Intervals:
        return [
            Interval(0, Infinity, False, False)
            if i in self._positive else
            (
                Interval(0, Infinity)
                if i in self._nonnegative else
                Interval(0, 0)
            )
            for i in range(self.matrix.nrows())
        ]

    def exists_orthogonal_vector(self, v) -> bool:
        return not (
            any(v[k] for k in self._positive)
            and (
                all(v[k] >= 0 for k in self._nonnegative)
                or all(v[k] <= 0 for k in self._nonnegative)
            )
        )

    def to_homogeneous(self) -> HomogeneousSystem:
        return self

    def certify_existence(self, reverse: bool = False, random: bool = False):
        result = zero_vector(self.matrix.base_ring(), self.matrix.nrows())

        if self._positive.stop == 0:
            return result
        for v in self.candidate_generator(dual=False, reverse=reverse, random=random):
            if self._solvable is False:
                raise ValueError("No solution exists!")
            for w in [v, -v]:
                if all(w[i] >= 0 for i in self._nonnegative) and all(w[i] == 0 for i in self._zero):
                    result += w
                    if all(result[i] > 0 for i in self._positive):
                        self._solvable = True
                        return result
                    break

        self._solvable = False
        raise ValueError("No solution exists!")

    def solve(self, reverse: bool = False, random: bool = False):
        r"""
        Compute a solution if existent.

        This approach sums up positive elementary vectors in the row space.
        It doesn't use division.

        .. NOTE::

            If no solution exists, and ``random`` is true, this method will never finish.
        """
        return solve_without_division(self.matrix, self.certify_existence(reverse=reverse, random=random))


# class HomogeneousSystemCocircuits(HomogeneousSystem):
#     r"""
#     A class for homogeneous linear inequality systems

#     ``A x > 0``, ``B x >= 0``, ``C x = 0``

#     Certifying makes use of cocircuits instead of elementary vectors
#     """
#     def __init__(self, A: Matrix, B: Matrix, C: Matrix, result: bool = None) -> None:
#         super().__init__(A, B, C, result=result)

#         # self._evs = Cocircuits(self.matrix.T)
#         self.om = OrientedMatroid(self.matrix.T)

#         if len(self.positive) == 1:
#             self._evs.set_combinations(CombinationsIncluding(self._evs.length, self._evs.rank + 1, self.positive))

#     def exists_orthogonal_vector(self, v) -> bool:
#         return not (
#             any(v[k] for k in self.positive)
#             and all(v[k] >= 0 for k in self.nonnegative)
#         )

#     def certify_existence(self, reverse: bool = False, random: bool = False) -> SignVector:
#         r"""
#         Compute a solution if existent.

#         This approach sums up positive elementary vectors in the row space.

#         .. NOTE::

#             If no solution exists, and ``random`` is true, this method will never finish.
#         """
#         lower = sign_vector(
#             len(self.positive) * [1] + (self.matrix.nrows() - len(self.positive)) * [0]
#         )
#         upper = sign_vector(
#             len(self.nonnegative) * [1] + (self.matrix.nrows() - len(self.nonnegative)) * [0]
#         )
#         result = zero_sign_vector(self.matrix.nrows())

#         if result >= lower:
#             return result
#         for v in self.candidate_generator(dual=False, reverse=reverse, random=random):
#             if self._solvable is False:
#                 raise ValueError("No solution exists!")
#             if v <= upper:
#                 result &= v
#                 if result >= lower:
#                     self._solvable = True
#                     return result

#         self._solvable = False
#         raise ValueError("No solution exists!")

#     def solve(self, reverse: bool = False, random: bool = False):
#         raise ValueError("Can't solve using cocircuits!")


def inhomogeneous_from_general(system: LinearInequalitySystem) -> InhomogeneousSystem:
    r"""Translate a general system into an inhomogeneous system."""
    A_list = []
    B_list = []
    b_list = []
    c_list = []

    for line, interval in zip(system.matrix, system._intervals):
        if interval.infimum() != -Infinity:
            if interval.infimum() in interval:
                A_list.append(-line)
                b_list.append(-interval.infimum())
            else:
                B_list.append(-line)
                c_list.append(-interval.infimum())
        if interval.supremum() != Infinity:
            if interval.supremum() in interval:
                A_list.append(line)
                b_list.append(interval.supremum())
            else:
                B_list.append(line)
                c_list.append(interval.supremum())

    return InhomogeneousSystem(
        Matrix(len(A_list), system.matrix.ncols(), A_list),
        Matrix(len(B_list), system.matrix.ncols(), B_list),
        vector(b_list),
        vector(c_list),
        result=system.result
    )


def homogeneous_from_general(system: LinearInequalitySystem) -> HomogeneousSystem:
    r"""
    Convert a general system to a homogeneous system.

    EXAMPLE::

        sage: from vectors_in_intervals import *
        sage: M = matrix([[1, 0], [0, 1], [1, 1]])
        sage: lower_bounds = [2, 5, 0]
        sage: upper_bounds = [5, oo, 0]
        sage: lower_bounds_closed = [True, True, True]
        sage: upper_bounds_closed = [False, False, True]
        sage: I = Intervals.from_bounds(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)
        sage: S = LinearInequalitySystem(M, I)
        sage: homogeneous_from_general(S)
        [ 1  0 -5]
        [ 0  0 -1]
        [--------]
        [-1  0  2]
        [ 0 -1  5]
        [--------]
        [ 1  1  0] x in [(0, +oo), (0, +oo), [0, +oo), [0, +oo), {0}]
    """
    A_list = []
    B_list = []
    C_list = []

    length = system.matrix.ncols()

    for line, interval in zip(system.matrix, system._intervals):
        if interval.infimum() == interval.supremum():
            C_list.append(list(line) + [-interval.infimum()])
            continue
        if interval.infimum() != -Infinity:
            if interval.infimum() in interval:
                B_list.append(list(-line) + [interval.infimum()])
            else:
                A_list.append(list(-line) + [interval.infimum()])
        if interval.supremum() != Infinity:
            if interval.supremum() in interval:
                B_list.append(list(line) + [-interval.supremum()])
            else:
                A_list.append(list(line) + [-interval.supremum()])

    A_list.append([0] * length + [-1])

    return HomogeneousSystem(
        Matrix(len(A_list), length + 1, A_list),
        Matrix(len(B_list), length + 1, B_list),
        Matrix(len(C_list), length + 1, C_list),
        result=system.result
    )


def homogeneous_from_inhomogeneous(system: InhomogeneousSystem) -> HomogeneousSystem:
    r"""Convert an inhomogeneous system to a homogeneous system."""
    return HomogeneousSystem(
        Matrix.block([[system.B, Matrix(len(system.c), 1, -system.c)], [zero_matrix(1, system.A.ncols()), Matrix([[-1]])]]),
        Matrix.block([[system.A, Matrix(len(system.b), 1, -system.b)]]),
        Matrix(0, system.A.ncols() + 1),
        result=system.result
    )
