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
    sage: S.certify()
    (True, (5, 15, 1, 2, 1, 0, 0))
    sage: S.certify(reverse=True)
    (True, (3, 7, 1, 1, 0, 0, 0))
    sage: # S.certify_parallel() # TODO SignalError: Segmentation Fault
    (True, (3, 7, 1, 1, 0, 0, 0))

::

    sage: S.to_inhomogeneous()
    [ 1  0]
    [-1 -1]
    [ 1  1]
    [-----]
    [-1  0]
    [ 0 -1]
    [ 0  1] x in [(-oo, 5), (-oo, 0), (-oo, 8), (-oo, -2], (-oo, -5], (-oo, 5]]
    sage: S.to_homogeneous()
    [ 1  0 -5]
    [-1 -1  0]
    [ 1  1 -8]
    [ 0  0 -1]
    [--------]
    [-1  0  2]
    [ 0 -1  5]
    [ 0  1 -5]
    [--------] x in [(0, +oo), (0, +oo), (0, +oo), (0, +oo), [0, +oo), [0, +oo), [0, +oo)]

We consider yet another system::

    sage: A = matrix([[-1, -1]])
    sage: B = matrix([[1, 0], [1, 1]])
    sage: a = vector([0])
    sage: b = vector([1, 0])
    sage: S = InhomogeneousSystem(A, B, a, b)
    sage: S.certify()
    (False, (1, 0, 1))
    sage: S.certify_parallel()
    (False, (1, 0, 1))
    sage: S.certify_parallel(random=True)
    (False, (1, 0, 1))

::

    sage: S.to_homogeneous()
    [-1 -1  0]
    [ 0  0 -1]
    [--------]
    [ 1  0 -1]
    [ 1  1  0]
    [--------] x in [(0, +oo), (0, +oo), [0, +oo), [0, +oo)]

In the case of homogeneous systems, we can use cocircuits to certify::

    sage: A = matrix([[1, 2], [0, 1]])
    sage: B = matrix([[2, 3]])
    sage: C = matrix([[-1, 0]])
    sage: S = HomogeneousSystemCocircuits(A, B, C) # TODO: not implemented
    sage: S.solve() # TODO: not implemented
    Traceback (most recent call last):
    ...
    ValueError: Can't solve using cocircuits!
    sage: # S.certify() # TODO

Now, we consider the example::

    sage: A = matrix([[1, 0], [0, 1]])
    sage: B = matrix([[2, -3]])
    sage: C = matrix([[-1, -1]])
    sage: S = HomogeneousSystemCocircuits(A, B, C) # TODO: not implemented
    sage: S.certify() # TODO: not implemented
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

# from sign_vectors import SignVector, sign_vector, zero_sign_vector
# from sign_vectors.oriented_matroids import OrientedMatroid
from . import Intervals, Interval
from elementary_vectors.functions import ElementaryVectors
from .utility import CombinationsIncluding, solve_without_division


class LinearInequalitySystem(SageObject):
    r"""
    A class for linear inequality systems given by a matrix and intervals
    """
    def __init__(self, matrix: Matrix, intervals: Intervals = None, result: bool = None) -> None:
        if intervals is not None and matrix.nrows() != len(intervals):
            raise ValueError("Matrix row count and number of intervals must agree!")
        self._matrix = matrix
        self._intervals = intervals
        self.result = result
        self._evs = ElementaryVectors(self.matrix.T)
        self._solvable = None

    def _repr_(self) -> str:
        return str(self.matrix) + " x in " + str(self.intervals)

    @property
    def matrix(self) -> Matrix:
        r"""Return the corresponding matrix."""
        return self._matrix

    @property
    def intervals(self) -> Intervals:
        r"""Return the corresponding intervals."""
        if self._intervals is None:
            self._intervals = self._compute_intervals()
        return self._intervals

    def _compute_intervals(self) -> Intervals:
        raise NotImplementedError("This method should be implemented in subclasses.")

    def to_homogeneous(self) -> HomogeneousSystem:
        r"""Return the equivalent homogeneous system."""
        matrix1_list = []
        matrix2_list = []
        matrix3_list = []

        length = self.matrix.ncols()

        for line, interval in zip(self.matrix, self.intervals):
            if interval.infimum() == interval.supremum():
                matrix3_list.append(list(line) + [-interval.infimum()])
                continue
            if interval.infimum() != -Infinity:
                if interval.infimum() in interval:
                    matrix2_list.append(list(-line) + [interval.infimum()])
                else:
                    matrix1_list.append(list(-line) + [interval.infimum()])
            if interval.supremum() != Infinity:
                if interval.supremum() in interval:
                    matrix2_list.append(list(line) + [-interval.supremum()])
                else:
                    matrix1_list.append(list(line) + [-interval.supremum()])

        matrix1_list.append([0] * length + [-1])

        return HomogeneousSystem(
            Matrix(len(matrix1_list), length + 1, matrix1_list),
            Matrix(len(matrix2_list), length + 1, matrix2_list),
            Matrix(len(matrix3_list), length + 1, matrix3_list),
            result=self.result
        )

    def to_inhomogeneous(self) -> InhomogeneousSystem:
        r"""Return the equivalent inhomogeneous system."""
        matrix_strict_list = []
        matrix_nonstrict_list = []
        vector_strict_list = []
        vector_nonstrict_list = []

        for line, interval in zip(self.matrix, self.intervals):
            if interval.infimum() != -Infinity:
                if interval.infimum() in interval:
                    matrix_nonstrict_list.append(-line)
                    vector_nonstrict_list.append(-interval.infimum())
                else:
                    matrix_strict_list.append(-line)
                    vector_strict_list.append(-interval.infimum())
            if interval.supremum() != Infinity:
                if interval.supremum() in interval:
                    matrix_nonstrict_list.append(line)
                    vector_nonstrict_list.append(interval.supremum())
                else:
                    matrix_strict_list.append(line)
                    vector_strict_list.append(interval.supremum())

        return InhomogeneousSystem(
            Matrix(len(matrix_strict_list), self.matrix.ncols(), matrix_strict_list),
            Matrix(len(matrix_nonstrict_list), self.matrix.ncols(), matrix_nonstrict_list),
            vector(vector_strict_list),
            vector(vector_nonstrict_list),
            result=self.result
        )

    def dual(self) -> LinearInequalitySystem:
        r"""Return the dual linear inequality system."""
        return self.to_inhomogeneous().dual()

    def with_intervals(self, intervals: Intervals) -> LinearInequalitySystem:
        r"""
        Return a copy of this system with different intervals.

        TESTS::

            sage: from vectors_in_intervals import *
            sage: M = matrix([[1, 0], [0, 1], [1, 1], [0, 1]])
            sage: lower_bounds = [2, 5, 0, -oo]
            sage: upper_bounds = [5, oo, 8, 5]
            sage: lower_bounds_closed = [True, True, False, False]
            sage: upper_bounds_closed = [False, False, False, True]
            sage: I = Intervals.from_bounds(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)
            sage: S = LinearInequalitySystem(M, I)
            sage: S.certify()
            (True, (5, 15, 1, 2, 1, 0, 0))
            sage: S.with_intervals(I).certify()
            (True, (5, 15, 1, 2, 1, 0, 0))
            sage: S.with_intervals(Intervals.from_bounds([2, 6, 0, -oo], [5, oo, 8, 5])).certify()
            (False, (0, -1, 0, 1))
            sage: S.with_intervals(Intervals.from_bounds([2, 5, 0, -oo], [5, 5, 8, 5])).certify()
            (True, (1, 0, 3, 7, 1, 0, 0))
        """
        system = LinearInequalitySystem(self.matrix, intervals, result=self.result)
        if self.__class__ is LinearInequalitySystem:
            system._evs = self._evs
        return system

    def _evs_generator(self, dual: bool = True, reverse: bool = False, random: bool = False) -> Generator:
        r"""Return a generator of elementary vectors."""
        if random:
            while True:
                yield self._evs.random_element(dual=dual)
        else:
            yield from self._evs.generator(dual=dual, reverse=reverse)

    def _exists_orthogonal_vector(self, v: vector) -> bool:
        r"""Check if an orthogonal vector exists in the intervals."""
        lower_product = 0
        upper_product = 0
        lower_product_attainable = True
        upper_product_attainable = True

        for entry, interval in zip(v, self.intervals):
            if interval.is_empty():
                return False
            if not entry:
                continue
            bound1 = interval.infimum() if entry > 0 else interval.supremum()
            bound2 = interval.supremum() if entry > 0 else interval.infimum()
            lower_product += entry * bound1
            upper_product += entry * bound2
            lower_product_attainable &= bound1 in interval
            upper_product_attainable &= bound2 in interval

        if lower_product > 0:
            return False
        if upper_product < 0:
            return False
        if lower_product == 0 and not lower_product_attainable:
            return False
        if upper_product == 0 and not upper_product_attainable:
            return False
        return True

    def _certify_nonexistence(self, reverse: bool = False, random: bool = False):
        r"""
        Certify nonexistence of solutions.

        Otherwise, a ``ValueError`` is raised.

        .. NOTE::

            If a solution exists and ``random`` is set to true, this method will never finish.
        """
        for v in self._evs_generator(reverse=reverse, random=random):
            if self._solvable:
                break
            if not self._exists_orthogonal_vector(v):
                self._solvable = False
                return v
        self._solvable = True
        raise ValueError("A solution exists!")

    def _certify_existence(self, reverse: bool = False, random: bool = False):
        r"""
        Certify existence of a solution if one exists.

        Otherwise, a ``ValueError`` is raised.

        .. NOTE::

            If no solution exists and ``random`` is set to true, this method will never finish.
        """
        return self.to_homogeneous()._certify_existence(reverse=reverse, random=random)

    def certify(self, reverse: bool = False) -> tuple:
        r"""Return a boolean and a certificate for solvability."""
        try:
            return False, self._certify_nonexistence(reverse=reverse)
        except ValueError:
            return True, self._certify_existence(reverse=reverse)

    def certify_parallel(self, reverse: bool = False, random: bool = False) -> tuple:
        r"""
        Return a boolean and a certificate for solvability.

        Attempts to find a solution and certify nonexistence in parallel.
        """
        with concurrent.futures.ThreadPoolExecutor() as executor:
            done, not_done = concurrent.futures.wait(
                [
                    executor.submit(lambda: (False, self._certify_nonexistence(reverse=reverse, random=random))),
                    executor.submit(lambda: (True, self._certify_existence(reverse=reverse, random=random)))
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

    def has_solution(self, reverse: bool = False, random: bool = False) -> bool:
        r"""
        Check whether a solution exists.

        If ``random`` is true, this method will never finish if a solution exists.
        """
        try:
            self._certify_nonexistence(reverse=reverse, random=random)
            return False
        except ValueError:
            return True

    def solve(self, reverse: bool = False, random: bool = False):
        r"""
        Compute a solution for this linear inequality system.

        If no solution exists, a ``ValueError`` is raised.
        """
        solution = self.to_homogeneous().solve(reverse=reverse, random=random)
        return solution[:-1] / solution[-1]


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
        sage: S.certify()
        (True, (1, 1, 1, 0, 0))
    """
    def __init__(self, matrix_strict: Matrix, matrix_nonstrict: Matrix, matrix_zero: Matrix, result: bool = None) -> None:
        super().__init__(Matrix.block([[matrix_strict], [matrix_nonstrict], [matrix_zero]]), None, result=result)
        self._length_strict = matrix_strict.nrows()
        self._length_nonstrict = matrix_nonstrict.nrows()
        self._length_zero = matrix_zero.nrows()

        # self._evs._set_combinations_row_space(Combinations(range(A.nrows() + B.nrows()), self._evs.length - self._evs.rank + 1))

        if self._length_strict == 1:
            self._evs._set_combinations_kernel(CombinationsIncluding(self._evs.length, self._evs.rank + 1, range(self._length_strict)))

    def _compute_intervals(self) -> Intervals:
        return [
            Interval.open(0, Infinity)
            if i < self._length_strict else
            (Interval.closed(0, Infinity) if i < self._length_strict + self._length_nonstrict else Interval.closed(0, 0))
            for i in range(self.matrix.nrows())
        ]

    def to_homogeneous(self) -> HomogeneousSystem:
        return self

    def dual(self) -> HomogeneousSystem:
        return HomogeneousSystem(
            Matrix.block([[Matrix.ones(1, self._length_strict), Matrix.zero(1, self._length_nonstrict), Matrix.zero(1, self._length_zero)]]),
            Matrix.block([
                [Matrix.identity(self._length_strict), Matrix.zero(self._length_strict, self._length_nonstrict), Matrix.zero(self._length_strict, self._length_zero)],
                [Matrix.zero(self._length_nonstrict, self._length_strict), Matrix.identity(self._length_nonstrict), Matrix.zero(self._length_nonstrict, self._length_zero)]
            ]),
            self.matrix.T,
            result=not self.result
        )

    def _exists_orthogonal_vector(self, v: vector) -> bool:
        if all(v[k] == 0 for k in range(self._length_strict)):
            return True
        if all(v[k] >= 0 for k in range(self._length_strict + self._length_nonstrict)):
            return False
        if all(v[k] <= 0 for k in range(self._length_strict + self._length_nonstrict)):
            return False
        return True

    def _certify_existence(self, reverse: bool = False, random: bool = False):
        certificate = zero_vector(self.matrix.base_ring(), self.matrix.nrows())

        if self._length_strict == 0:
            return certificate
        for v in self._evs_generator(dual=False, reverse=reverse, random=random):
            if self._solvable is False:
                raise ValueError("System is marked as unsolvable!")
            for w in [v, -v]:
                if any(w[i] < 0 for i in range(self._length_strict + self._length_nonstrict)):
                    continue
                if any(w[i] != 0 for i in range(self._length_strict + self._length_nonstrict, self.matrix.nrows())):
                    continue
                certificate += w
                if all(certificate[i] > 0 for i in range(self._length_strict)):
                    self._solvable = True
                    return certificate
                break

        self._solvable = False
        raise ValueError("Couldn't construct a solution. No solution exists!")

    def solve(self, reverse: bool = False, random: bool = False):
        r"""
        Compute a solution if existent.

        This approach sums up positive elementary vectors in the row space.
        It doesn't use division.

        .. NOTE::

            If no solution exists, and ``random`` is true, this method will never finish.
        """
        return solve_without_division(self.matrix, self._certify_existence(reverse=reverse, random=random))


class InhomogeneousSystem(LinearInequalitySystem):
    r"""
    A class for inhomogeneous linear inequality systems

    ``A x < a``, ``B x <= b``
    """
    def __init__(
            self,
            matrix_strict: Matrix,
            matrix_nonstrict: Matrix,
            vector_strict: vector,
            vector_nonstrict: vector,
            result: bool = None
    ) -> None:
        super().__init__(Matrix.block([[matrix_strict], [matrix_nonstrict]]), None, result=result)
        self._matrix_strict = matrix_strict
        self._matrix_nonstrict = matrix_nonstrict
        self._vector_strict = vector_strict
        self._vector_nonstrict = vector_nonstrict

    def _compute_intervals(self) -> Intervals:
        return [Interval.open(-Infinity, ai) for ai in self._vector_strict] + [Interval.closed(-Infinity, bi) for bi in self._vector_nonstrict]

    def to_homogeneous(self) -> HomogeneousSystem:
        return HomogeneousSystem(
            Matrix.block([[self._matrix_strict, Matrix(len(self._vector_strict), 1, -self._vector_strict)], [zero_matrix(1, self._matrix_nonstrict.ncols()), Matrix([[-1]])]]),
            Matrix.block([[self._matrix_nonstrict, Matrix(len(self._vector_nonstrict), 1, -self._vector_nonstrict)]]),
            Matrix(0, self._matrix_nonstrict.ncols() + 1),
            result=self.result
        )

    def to_inhomogeneous(self):
        return self

    def dual(self) -> HomogeneousSystem:
        length_strict = self._matrix_strict.nrows()
        length_nonstrict = self._matrix_nonstrict.nrows()
        length = self._matrix_nonstrict.ncols()
        return HomogeneousSystem(
            Matrix.block([
                [Matrix.zero(1, length_nonstrict), Matrix.ones(1, length_strict + 1)]
            ]),
            Matrix.identity(length_nonstrict + length_strict + 1),
            Matrix.block([
                [self._matrix_nonstrict.T, self._matrix_strict.T, Matrix.zero(length, 1)],
                [-self._vector_nonstrict.row(), -self._vector_strict.row(), Matrix([[-1]])]
            ]),
            result=not self.result
        )

    def _exists_orthogonal_vector(self, v: vector) -> bool:
        length_strict = len(self._vector_strict)
        length_nonstrict = len(self._vector_nonstrict)

        if v == 0:
            return True

        positive = None
        if all(vk >= 0 for vk in v):
            positive = True
        elif all(vk <= 0 for vk in v):
            positive = False
        if positive is None:
            return True

        v_strict = vector(v[k] for k in range(length_strict))
        v_nonstrict = vector(v[k + length_strict] for k in range(length_nonstrict))

        scalarproduct = v_strict * self._vector_strict + v_nonstrict * self._vector_nonstrict

        if positive and scalarproduct < 0:
            return False
        if not positive and scalarproduct > 0:
            return False
        if scalarproduct == 0 and v_strict != 0:
            return False
        return True


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

#     def _certify_existence(self, reverse: bool = False, random: bool = False) -> SignVector:
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
#         certificate = zero_sign_vector(self.matrix.nrows())

#         if certificate >= lower:
#             return certificate
#         for v in self._evs_generator(dual=False, reverse=reverse, random=random):
#             if self._solvable is False:
#                 raise ValueError("No solution exists!")
#             if v <= upper:
#                 certificate &= v
#                 if certificate >= lower:
#                     self._solvable = True
#                     return certificate

#         self._solvable = False
#         raise ValueError("No solution exists!")

#     def solve(self, reverse: bool = False, random: bool = False):
#         raise ValueError("Can't solve using cocircuits!")
