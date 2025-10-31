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
    (0, +oo) x (0, +oo) x [0, +oo) x {0}
    sage: S.find_solution()
    (0, 1)
    sage: S.certify()
    (True, (0, 1))
    sage: S.certify(random=True)
    (True, (0, 1))

We consider another system::

    sage: M = matrix([[1, 0], [0, 1], [1, 1], [0, 1]])
    sage: lower_bounds = [2, 5, 0, -oo]
    sage: upper_bounds = [5, oo, 8, 5]
    sage: lower_bounds_closed = [True, True, False, False]
    sage: upper_bounds_closed = [False, False, False, True]
    sage: I = Intervals.from_bounds(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)
    sage: S = LinearInequalitySystem(M, I)
    sage: S.find_solution()
    (5/2, 5)
    sage: S.certify()
    (True, (5/2, 5))

::

    sage: S.to_inhomogeneous()
    [ 1  0]
    [-1 -1]
    [ 1  1]
    [-----]
    [-1  0]
    [ 0 -1]
    [ 0  1] x in (-oo, 5) x (-oo, 0) x (-oo, 8) x (-oo, -2] x (-oo, -5] x (-oo, 5]
    sage: S.to_homogeneous()
    [ 1  0 -5]
    [-1 -1  0]
    [ 1  1 -8]
    [ 0  0 -1]
    [--------]
    [-1  0  2]
    [ 0 -1  5]
    [ 0  1 -5]
    [--------] x in (0, +oo) x (0, +oo) x (0, +oo) x (0, +oo) x [0, +oo) x [0, +oo) x [0, +oo)

We consider yet another system::

    sage: A = matrix([[-1, -1]])
    sage: B = matrix([[1, 0], [1, 1]])
    sage: a = vector([0])
    sage: b = vector([1, 0])
    sage: S = InhomogeneousSystem(A, B, a, b)
    sage: S.certify()
    (False, (1, 0, 1))
    sage: S.certify(random=True)
    (False, (1, 0, 1))
    sage: S.certify_nonexistence()
    (1, 0, 1)

::

    sage: S.to_homogeneous()
    [-1 -1  0]
    [ 0  0 -1]
    [--------]
    [ 1  0 -1]
    [ 1  1  0]
    [--------] x in (0, +oo) x (0, +oo) x [0, +oo) x [0, +oo)

TESTS::

    sage: A = matrix([[1, 1]])
    sage: B = matrix([[0, 1]])
    sage: C = matrix([[1, -1]])
    sage: S = HomogeneousSystem(A, B, C)
    sage: S
    [ 1  1]
    [-----]
    [ 0  1]
    [-----]
    [ 1 -1] x in (0, +oo) x [0, +oo) x {0}
    sage: S.to_inhomogeneous()
    [-1 -1]
    [-----]
    [ 0 -1]
    [-1  1]
    [ 1 -1] x in (-oo, 0) x (-oo, 0] x (-oo, 0] x (-oo, 0]
    sage: S.to_inhomogeneous().to_homogeneous()
    [-1 -1  0]
    [ 0  0 -1]
    [--------]
    [ 0 -1  0]
    [-1  1  0]
    [ 1 -1  0]
    [--------] x in (0, +oo) x (0, +oo) x [0, +oo) x [0, +oo) x [0, +oo)
    sage: S.to_inhomogeneous().to_homogeneous().to_inhomogeneous()
    [ 1  1  0]
    [ 0  0  1]
    [--------]
    [ 0  1  0]
    [ 1 -1  0]
    [-1  1  0] x in (-oo, 0) x (-oo, 0) x (-oo, 0] x (-oo, 0] x (-oo, 0]
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

from typing import Iterator

from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Manager

from sage.matrix.constructor import Matrix, zero_matrix
from sage.modules.free_module_element import vector, zero_vector
from sage.rings.infinity import Infinity
from sage.structure.sage_object import SageObject

from . import Intervals, Interval
from elementary_vectors.functions import ElementaryVectors
from .utility import CombinationsIncluding, solve_without_division


class LinearInequalitySystem(SageObject):
    r"""A class for linear inequality systems given by a matrix and intervals."""
    def __init__(self, matrix: Matrix, intervals: Intervals = None) -> None:
        if intervals is not None and matrix.nrows() != len(intervals):
            raise ValueError("Matrix row count and number of intervals must agree!")
        self._matrix = matrix
        self._intervals = intervals
        self._evs = ElementaryVectors(self.matrix.T)

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
            Matrix(len(matrix3_list), length + 1, matrix3_list)
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
            vector(vector_nonstrict_list)
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
            (True, (5/2, 5))
            sage: S.with_intervals(I).certify()
            (True, (5/2, 5))
            sage: S.with_intervals(Intervals.from_bounds([2, 6, 0, -oo], [5, oo, 8, 5])).certify()
            (False, (0, 1, 0, -1))
            sage: S.with_intervals(Intervals.from_bounds([2, 5, 0, -oo], [5, 5, 8, 5])).certify()
            (True, (2, 5))
        """
        system = LinearInequalitySystem(self.matrix, intervals)
        if self.__class__ is LinearInequalitySystem:
            system._evs = self._evs
        return system

    def _evs_generator(self, kernel: bool = True, random: bool = False, reverse: bool = False) -> Iterator[vector]:
        r"""Return a generator of elementary vectors."""
        if random:
            while True:
                yield self._evs.random_element(kernel=kernel)
        else:
            yield from self._evs.generator(kernel=kernel, reverse=reverse)

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

    def certify_nonexistence(self, random: bool = False, iteration_limit: int = 10000) -> vector:
        r"""
        Certify nonexistence of a solution if no solution exists.

        INPUT:

        - ``random`` -- if true, tries random elementary vectors
        - ``iteration_limit`` -- maximum number of iterations (by default 10000). If -1, unlimited.

        OUTPUT:
        A vector certifying that no solution exists.

        .. NOTE::

            - If the iteration limit is reached, a ``MaxIterationsReachedError`` is raised.
            - If a solution exists, a ``ValueError`` is raised.
            - If a solution exists, the iteration limit is ``-1`` _and_ ``random`` is true, this leads to an endless loop.
        """
        return self._certify_nonexistence(random, reverse=True, iteration_limit=iteration_limit, stop_event=None)

    def _certify_nonexistence(self, random: bool, reverse: bool, iteration_limit: int, stop_event=None) -> vector:
        r"""
        Certify nonexistence of a solution.

        Otherwise, a ``ValueError`` is raised.

        .. NOTE::

            Raises an exception if the maximum number of iterations is reached.
        """
        for i, v in enumerate(self._evs_generator(random=random, reverse=reverse)):
            if stop_event is not None and stop_event.is_set():
                raise ProcessStoppedError("Process was stopped because another process found a solution.")
            if iteration_limit != -1 and i >= iteration_limit:
                raise MaxIterationsReachedError("Reached maximum number of iterations! Does a solution exist?")
            if v is None:
                continue
            if not self._exists_orthogonal_vector(v):
                return v
        raise ValueError("A solution exists!")

    def find_solution(self, random: bool = False, iteration_limit: int = 10000) -> vector:
        r"""
        Compute a solution if existent.

        INPUT:

        - ``random`` -- if true, the returned sum consists of random elementary vectors
        - ``iteration_limit`` -- maximum number of iterations (by default 10000). If -1, unlimited.

        OUTPUT:
        A vector (as a sum of elementary vectors) solving the system.

        .. NOTE::

            - If the iteration limit is reached, a ``MaxIterationsReachedError`` is raised.
            - If no solution exists, a ``ValueError`` is raised.
            - If no solution exists, the iteration limit is ``-1`` _and_ ``random`` is true, this leads to an endless loop.

        .. SEEALSO::

            :meth:`certify`
        """
        return self._find_solution(random, reverse=False, iteration_limit=iteration_limit, stop_event=None)

    def _find_solution(self, random: bool, reverse: bool, iteration_limit: int, stop_event=None) -> vector:
        solution = self.to_homogeneous()._find_solution(random=random, reverse=reverse, iteration_limit=iteration_limit, stop_event=stop_event)
        return solution[:-1] / solution[-1]

    def certify(self, random: bool = False, iteration_limit: int = -1) -> tuple[bool, vector]:
        r"""
        Return a boolean and a certificate for solvability.

        Both existence and nonexistence are checked in parallel.

        INPUT:

        - ``random`` -- if true, elementary vectors are generated randomly
        - ``iteration_limit`` -- maximum number of iterations for each process (by default unlimited)

        OUTPUT:
        A tuple ``(exists, certificate)`` where ``exists`` is a boolean indicating whether a solution exists,
        and ``certificate`` is either a solution (if ``exists`` is true) or a vector certifying nonexistence (if ``exists`` is false).
        """
        return self._certify_parallel(random=random, iteration_limit=iteration_limit)

    def _certify_parallel(self, random: bool, iteration_limit: int) -> tuple[bool, vector]:
        r"""Return a boolean and a certificate for solvability in parallel."""
        stop_event = Manager().Event()
        with ProcessPoolExecutor(max_workers=2) as executor:
            futures = {
                executor.submit(self._certify_nonexistence, random=random, reverse=True, iteration_limit=iteration_limit, stop_event=stop_event): False,
                executor.submit(self._find_solution, random=random, reverse=False, iteration_limit=iteration_limit, stop_event=stop_event): True,
            }
            for future in as_completed(futures):
                try:
                    result = future.result()
                    stop_event.set()
                    return (futures[future], result)
                except (ValueError, MaxIterationsReachedError, ProcessStoppedError):
                    pass

        raise MaxIterationsReachedError("Both processes exceeded the maximum number of iterations.")


class HomogeneousSystem(LinearInequalitySystem):
    r"""
    A class for homogeneous linear inequality systems.

    ``A x > 0``, ``B x >= 0``, ``C x = 0``

    TESTS::

        sage: from vectors_in_intervals import *
        sage: A = matrix([[0, 1], [0, 1], [0, 1]])
        sage: B = zero_matrix(0, 2)
        sage: C = matrix([[1, 1], [0, 0]])
        sage: S = HomogeneousSystem(A, B, C)
        sage: S.certify()
        (True, (-1, 1))
    """
    def __init__(self, matrix_strict: Matrix, matrix_nonstrict: Matrix, matrix_zero: Matrix) -> None:
        super().__init__(Matrix.block([[matrix_strict], [matrix_nonstrict], [matrix_zero]]), None)
        self._length_strict = matrix_strict.nrows()
        self._length_nonstrict = matrix_nonstrict.nrows()
        self._length_zero = matrix_zero.nrows()

        if self._length_strict == 1:
            self._evs._set_combinations_kernel(CombinationsIncluding(self._evs.length, self._evs.rank + 1, [0]))

    def _compute_intervals(self) -> Intervals:
        return Intervals([
            Interval.open(0, Infinity)
            if i < self._length_strict else
            (Interval.closed(0, Infinity) if i < self._length_strict + self._length_nonstrict else Interval.closed(0, 0))
            for i in range(self.matrix.nrows())
        ])

    def to_homogeneous(self) -> HomogeneousSystem:
        return self

    def dual(self) -> HomogeneousSystem:
        return HomogeneousSystem(
            Matrix.block([[Matrix.ones(1, self._length_strict), Matrix.zero(1, self._length_nonstrict), Matrix.zero(1, self._length_zero)]]),
            Matrix.block([
                [Matrix.identity(self._length_strict), Matrix.zero(self._length_strict, self._length_nonstrict), Matrix.zero(self._length_strict, self._length_zero)],
                [Matrix.zero(self._length_nonstrict, self._length_strict), Matrix.identity(self._length_nonstrict), Matrix.zero(self._length_nonstrict, self._length_zero)]
            ]),
            self.matrix.T
        )

    def _exists_orthogonal_vector(self, v: vector) -> bool:
        if all(v[k] == 0 for k in range(self._length_strict)):
            return True
        if all(v[k] >= 0 for k in range(self._length_strict + self._length_nonstrict)):
            return False
        if all(v[k] <= 0 for k in range(self._length_strict + self._length_nonstrict)):
            return False
        return True

    def _find_solution(self, random: bool, reverse: bool, iteration_limit: int, stop_event=None) -> vector:
        return solve_without_division(self.matrix, self._certify_existence(random=random, reverse=reverse, iteration_limit=iteration_limit, stop_event=stop_event))

    def _certify_existence(self, random: bool, reverse: bool, iteration_limit: int, stop_event=None) -> vector:
        certificate = zero_vector(self.matrix.base_ring(), self.matrix.nrows())

        if self._length_strict == 0:
            return certificate
        for i, v in enumerate(self._evs_generator(kernel=False, random=random, reverse=reverse)):
            if stop_event is not None and stop_event.is_set():
                raise ProcessStoppedError("Process was stopped because another process found a certificate for nonexistence.")
            if iteration_limit != -1 and i >= iteration_limit:
                raise MaxIterationsReachedError("Reached maximum number of iterations! Is system unsolvable?")
            if v is None:
                continue
            for w in [v, -v]:
                if any(w[i] < 0 for i in range(self._length_strict + self._length_nonstrict)):
                    continue
                if any(w[i] != 0 for i in range(self._length_strict + self._length_nonstrict, self.matrix.nrows())):
                    continue
                certificate += w
                if all(certificate[i] > 0 for i in range(self._length_strict)):
                    return certificate
                break

        raise ValueError("Couldn't construct a solution. No solution exists!")


class InhomogeneousSystem(LinearInequalitySystem):
    r"""
    A class for inhomogeneous linear inequality systems.

    ``A x < a``, ``B x <= b``
    """
    def __init__(self, matrix_strict: Matrix, matrix_nonstrict: Matrix, vector_strict: vector, vector_nonstrict: vector) -> None:
        super().__init__(Matrix.block([[matrix_strict], [matrix_nonstrict]]), None)
        self._matrix_strict = matrix_strict
        self._matrix_nonstrict = matrix_nonstrict
        self._vector_strict = vector_strict
        self._vector_nonstrict = vector_nonstrict

    def _compute_intervals(self) -> Intervals:
        return Intervals([Interval.open(-Infinity, ai) for ai in self._vector_strict] + [Interval.closed(-Infinity, bi) for bi in self._vector_nonstrict])

    def to_homogeneous(self) -> HomogeneousSystem:
        return HomogeneousSystem(
            Matrix.block([[self._matrix_strict, Matrix(len(self._vector_strict), 1, -self._vector_strict)], [zero_matrix(1, self._matrix_nonstrict.ncols()), Matrix([[-1]])]]),
            Matrix.block([[self._matrix_nonstrict, Matrix(len(self._vector_nonstrict), 1, -self._vector_nonstrict)]]),
            Matrix(0, self._matrix_nonstrict.ncols() + 1)
        )

    def to_inhomogeneous(self) -> InhomogeneousSystem:
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
            ])
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


class MaxIterationsReachedError(Exception):
    """Raised when the maximum number of iterations is reached."""


class ProcessStoppedError(Exception):
    """Raised when the process was stopped by another process."""
