r"""
Circuits and cocircuits
=======================

Elementary vectors are support minimal vectors of a subspace.
They are also called *circuits*.
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

from typing import Optional, List, Iterator

from sage.combinat.combination import Combinations
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import zero_vector, vector
from sage.structure.sage_object import SageObject

from .utility import is_constant


def circuits(matrix: Matrix, prevent_multiples: bool = True) -> List[vector]:
    r"""
    Compute the circuits of a matrix.

    INPUT:

    - ``matrix`` -- a matrix
    - ``prevent_multiples`` -- a boolean (default: ``True``)

    OUTPUT:
    Return the circuits of this matrix.
    These are the nonzero support-minimal elements in the kernel.

    .. SEEALSO::

        - :func:`cocircuits`

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
        sage: M
        [1 2 0 0]
        [0 1 2 3]
        sage: circuits(M)
        [(4, -2, 1, 0), (6, -3, 0, 1), (0, 0, -3, 2)]
        sage: circuits(M, prevent_multiples=False)
        [(4, -2, 1, 0), (6, -3, 0, 1), (0, 0, -3, 2), (0, 0, -6, 4)]

    Variables are also supported::

        sage: var('a, b')
        (a, b)
        sage: M = matrix([[1, 2, a, 0], [0, 1, 2, b]])
        sage: M
        [1 2 a 0]
        [0 1 2 b]
        sage: circuits(M)
        [(-a + 4, -2, 1, 0), (2*b, -b, 0, 1), (a*b, 0, -b, 2), (0, a*b, -2*b, -a + 4)]

    Matrices over the polynomial ring work, too::

        sage: R = PolynomialRing(ZZ, "x")
        sage: x = R.gen()
        sage: M = matrix([[1, 2, x, 0], [0, 1, 2, x]])
        sage: M
        [1 2 x 0]
        [0 1 2 x]
        sage: circuits(M)
        [(-x + 4, -2, 1, 0), (2*x, -x, 0, 1), (x^2, 0, -x, 2), (0, x^2, -2*x, -x + 4)]
        sage: R = PolynomialRing(ZZ, "x, y")
        sage: x, y = R.gens()
        sage: M = matrix([[x, y, 0, 0], [0, 1, 2, 3]])
        sage: M
        [x y 0 0]
        [0 1 2 3]
        sage: circuits(M)
        [(2*y, -2*x, x, 0), (3*y, -3*x, 0, x), (0, 0, -3*x, 2*x)]
    """
    return CircuitEnumerator(matrix).circuits(prevent_multiples=prevent_multiples)


def cocircuits(matrix: Matrix, prevent_multiples: bool = True) -> List[vector]:
    r"""
    Compute the cocircuits of a matrix.

    INPUT:
    - ``matrix`` -- a matrix
    - ``prevent_multiples`` -- a boolean (default: ``True``)

    OUTPUT:
    Return the cocircuits of this matrix.
    These are the nonzero support-minimal elements in the row space.

    .. SEEALSO::

        - :func:`circuits`

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
        sage: M
        [1 2 0 0]
        [0 1 2 3]
        sage: cocircuits(M)
        [(0, -1, -2, -3), (1, 0, -4, -6), (2, 4, 0, 0)]
        sage: cocircuits(M, prevent_multiples=False)
        [(0, -1, -2, -3), (1, 0, -4, -6), (2, 4, 0, 0), (3, 6, 0, 0)]
    """
    return CircuitEnumerator(matrix).cocircuits(prevent_multiples=prevent_multiples)


def circuit_generator(matrix: Matrix, prevent_multiples: bool = True, reverse: bool = False) -> Iterator[vector]:
    r"""
    Generator of circuits of a matrix.

    INPUT:
    - ``matrix`` -- a matrix
    - ``prevent_multiples`` -- a boolean (default: ``True``)
    - ``reverse`` -- a boolean (default: ``False``)

    OUTPUT:
    A generator of the circuits of this matrix.
    These are the nonzero support-minimal elements in the kernel.

    EXAMPLES::

        sage: from elementary_vectors.circuits import circuit_generator
        sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
        sage: M
        [1 2 0 0]
        [0 1 2 3]
        sage: list(circuit_generator(M))
        [(4, -2, 1, 0), (6, -3, 0, 1), (0, 0, -3, 2)]
        sage: list(circuit_generator(M, reverse=True))
        [(0, 0, -6, 4), (6, -3, 0, 1), (4, -2, 1, 0)]
    """
    return CircuitEnumerator(matrix).circuit_generator(prevent_multiples=prevent_multiples, reverse=reverse)


def cocircuit_generator(matrix: Matrix, prevent_multiples: bool = True, reverse: bool = False) -> Iterator[vector]:
    r"""
    Generator of cocircuits of a matrix.

    INPUT:
    - ``matrix`` -- a matrix
    - ``prevent_multiples`` -- a boolean (default: ``True``)
    - ``reverse`` -- a boolean (default: ``False``)

    OUTPUT:
    A generator of the cocircuits of this matrix.
    These are the nonzero support-minimal elements in the row space.

    EXAMPLES::

        sage: from elementary_vectors.circuits import cocircuit_generator
        sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
        sage: M
        [1 2 0 0]
        [0 1 2 3]
        sage: list(cocircuit_generator(M))
        [(0, -1, -2, -3), (1, 0, -4, -6), (2, 4, 0, 0)]
        sage: list(cocircuit_generator(M, reverse=True))
        [(3, 6, 0, 0), (1, 0, -4, -6), (0, -1, -2, -3)]
    """
    return CircuitEnumerator(matrix).cocircuit_generator(prevent_multiples=prevent_multiples, reverse=reverse)


def circuit_kernel_matrix(matrix: Matrix) -> Matrix:
    """
    Right kernel matrix based on circuits.

    OUTPUT:
    A right kernel matrix.
    Each row is a circuit of the input matrix.
    Therefore, the output matrix is division-free.
    Also works for symbolic matrices.

    .. NOTE::

        Raises a ``ValueError`` if the matrix has no constant nonzero maximal minor.

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[1, 0, 1, -1, 0], [0, 1, 1, 1, -1]])
        sage: M
        [ 1  0  1 -1  0]
        [ 0  1  1  1 -1]
        sage: circuit_kernel_matrix(M)
        [1 0 0 1 1]
        [0 1 0 0 1]
        [0 0 1 1 2]
        sage: M = matrix([[1, 1, -1, -1, 0], [2, 1, -1, 0, 1], [1, 1, 1, 1, 1]])
        sage: M
        [ 1  1 -1 -1  0]
        [ 2  1 -1  0  1]
        [ 1  1  1  1  1]
        sage: circuit_kernel_matrix(M)
        [-1  0  0 -1  2]
        [ 0 -1  1 -2  2]
        sage: var('a')
        a
        sage: M = matrix([[1, 0, 1, -1, 0], [0, 1, a, 1, -1]])
        sage: M
        [ 1  0  1 -1  0]
        [ 0  1  a  1 -1]
        sage: circuit_kernel_matrix(M)
        [    1     0     0     1     1]
        [    0     1     0     0     1]
        [    0     0     1     1 a + 1]

    TESTS::

        sage: circuit_kernel_matrix(identity_matrix(3, 3))
        []
        sage: _.dimensions()
        (0, 3)
        sage: circuit_kernel_matrix(matrix(0, 3))
        [1 0 0]
        [0 1 0]
        [0 0 1]

    ::

        sage: M = matrix([[0, 1], [0, 1]])
        sage: circuit_kernel_matrix(M)
        [1 0]
    """
    ce = CircuitEnumerator(matrix)

    rank = ce.rank
    length = ce.length

    if rank == length:
        return Matrix(matrix.base_ring(), 0, length)

    for indices_minor in Combinations(range(length - 1, -1, -1), rank):
        minor = ce.minor(indices_minor)
        if minor != 0 and is_constant(minor):
            return Matrix(ce.circuit(indices) for indices in ce._index_sets_from_minor(indices_minor, kernel=True))
    raise ValueError("Matrix has no constant nonzero maximal minor.")


def degenerate_circuits(matrix: Matrix) -> list[vector]:
    r"""
    Compute degenerate circuits of a matrix.

    OUTPUT:
    Return a list of degenerate circuits of the matrix.
    These are the nonzero support-minimal elements in the kernel
    with support smaller than rank + 1.

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[1, 0, 1, 0], [0, 0, 1, 1]])
        sage: M
        [1 0 1 0]
        [0 0 1 1]
        sage: degenerate_circuits(M)
        [(0, -1, 0, 0)]
        sage: M = matrix([[1, 1, 1, 0], [0, 1, 1, 1]])
        sage: M
        [1 1 1 0]
        [0 1 1 1]
        sage: degenerate_circuits(M)
        [(0, -1, 1, 0)]
    """
    return CircuitEnumerator(matrix).degenerate_circuits()


def degenerate_cocircuits(matrix: Matrix) -> list[vector]:
    r"""
    Compute degenerate cocircuits of a matrix.

    OUTPUT:
    Return a list of degenerate cocircuits of the matrix.
    These are the nonzero support-minimal elements in the row space
    with support smaller than rank - 1.

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[1, 0, 1, 0], [0, 0, 1, 1]])
        sage: M
        [1 0 1 0]
        [0 0 1 1]
        sage: degenerate_cocircuits(M)
        [(0, 0, -1, -1), (1, 0, 0, -1), (1, 0, 1, 0)]
        sage: M = matrix([[1, 1, 1, 0], [0, 1, 1, 1]])
        sage: M
        [1 1 1 0]
        [0 1 1 1]
        sage: degenerate_cocircuits(M)
        [(1, 0, 0, -1)]
    """
    return CircuitEnumerator(matrix).degenerate_cocircuits()


class CircuitEnumerator(SageObject):
    r"""
    A class used to compute circuits and cocircuits.

    If you want to compute circuits *and* cocircuits of the same matrix,
    use this class instead of the functions :func:`circuits` and :func:`cocircuits`
    since they share computed maximal minors.
    Also use this class if you want to compute individual circuits or cocircuits.

    .. NOTE::

        - Whenever a maximal minor is computed, it is stored in a dictionary for efficient reuse.
        - If the provided matrix is not full rank, it is replaced by a full-rank one.

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[1, 2, 4, 1, -1], [0, 1, 2, 3, 4]])
        sage: ce = CircuitEnumerator(M)
        sage: ce
        Circuit enumerator of 2x5 matrix
        sage: ce.circuits()
        [(0, -2, 1, 0, 0),
         (5, -3, 0, 1, 0),
         (9, -4, 0, 0, 1),
         (10, 0, -3, 2, 0),
         (18, 0, -4, 0, 2),
         (7, 0, 0, -4, 3),
         (0, 7, 0, -9, 5),
         (0, 0, 7, -18, 10)]
        sage: ce.circuits(prevent_multiples=False)
        [(0, -2, 1, 0, 0),
         (5, -3, 0, 1, 0),
         (9, -4, 0, 0, 1),
         (10, 0, -3, 2, 0),
         (18, 0, -4, 0, 2),
         (7, 0, 0, -4, 3),
         (0, 10, -5, 0, 0),
         (0, 18, -9, 0, 0),
         (0, 7, 0, -9, 5),
         (0, 0, 7, -18, 10)]
        sage: ce.cocircuits()
        [(0, -1, -2, -3, -4), (1, 0, 0, -5, -9), (3, 5, 10, 0, -7), (4, 9, 18, 7, 0)]
        sage: ce.cocircuits(prevent_multiples=False)
        [(0, -1, -2, -3, -4),
         (1, 0, 0, -5, -9),
         (2, 0, 0, -10, -18),
         (3, 5, 10, 0, -7),
         (4, 9, 18, 7, 0)]

    We compute individual elements::

        sage: ce.minor([0, 2])
        2
        sage: ce.circuit([0, 2, 3])
        (10, 0, -3, 2, 0)
        sage: ce.cocircuit([0])
        (0, -1, -2, -3, -4)
        sage: ce.random_circuit() # random
        (0, 0, 7, -18, 10)
        sage: ce.random_cocircuit() # random
        (3, 5, 10, 0, -7)

    Now, we consider an example that involves many zero minors::

        sage: M = matrix([[1, 2, 4, 0], [0, 1, 2, 0]])
        sage: M.minors(2)
        [1, 2, 0, 0, 0, 0]
        sage: ce = CircuitEnumerator(M)
        sage: ce.circuits()
        [(0, -2, 1, 0), (0, 0, 0, 1)]
    """
    def __init__(self, matrix: Matrix) -> None:
        try:
            self.matrix = matrix.matrix_from_rows(matrix.pivot_rows())
        except NotImplementedError as exc:
            if all(
                matrix.matrix_from_columns(indices).det() == 0
                for indices in Combinations(matrix.ncols(), matrix.nrows())
            ):
                raise ValueError("Provide a matrix with maximal rank.") from exc
            self.matrix = matrix
        self.rank, self.length = self.matrix.dimensions()
        self.ring = matrix.base_ring()
        self._minors = {}

        self._set_combinations_kernel()
        self._set_combinations_row_space()
        self._reset_set_for_preventing_multiples()

    def _repr_(self) -> str:
        return f"Circuit enumerator of {self.rank}x{self.length} matrix"

    def _set_combinations_kernel(self, combinations: Optional[Combinations] = None) -> None:
        r"""Set or reset combinations for elements in the kernel."""
        if combinations is None:
            self._combinations_kernel = Combinations(self.length, self.rank + 1)
        else:
            self._combinations_kernel = combinations

    def _set_combinations_row_space(self, combinations: Optional[Combinations] = None) -> None:
        r"""Set or reset combinations for elements in the row space."""
        if combinations is None:
            self._combinations_row_space = Combinations(self.length, self.rank - 1)
        else:
            self._combinations_row_space = combinations

    def _reset_set_for_preventing_multiples(self) -> None:
        self._zero_minors = set()
        self._marked_minors = set()

    def minor(self, indices: List[int], mark_if_zero: bool = False):
        r"""
        Compute a minor given by (sorted) indices.

        .. NOTE::

            The minor is cached for efficient reuse.

        TESTS::

            sage: from elementary_vectors import *
            sage: M = matrix([[1, 2, 4, 1, -1], [0, 1, 2, 3, 4]])
            sage: ce = CircuitEnumerator(M)
            sage: ce._minors
            {}
            sage: ce.minor([0, 1])
            1
            sage: ce._minors
            {(0, 1): 1}
            sage: ce.minor([2, 4])
            18
            sage: ce._minors
            {(0, 1): 1, (2, 4): 18}
            sage: ce.minor([0, 1, 2])
            Traceback (most recent call last):
            ...
            ValueError: Indices (0, 1, 2) should have size 2 and not 3.
        """
        indices = tuple(indices)
        minor = self._minors.get(indices)
        if minor is None:
            try:
                minor = self._compute_minor(indices)
                self._minors[indices] = minor
            except ValueError as e:
                raise ValueError(f"Indices {indices} should have size {self.rank} and not {len(indices)}.") from e
        if mark_if_zero and minor == 0:
            self._zero_minors.add(indices)
        return minor

    def _compute_minor(self, indices: tuple[int]):
        return self.matrix.matrix_from_columns(indices).det()

    def compute_minors(self) -> None:
        r"""Compute all maximal minors of the matrix."""
        for indices in Combinations(self.length, self.rank):
            self.minor(indices)

    def circuit(self, indices: List[int], prevent_multiple: bool = False) -> vector:
        r"""
        Compute the circuit corresponding to a list of indices.

        INPUT:

        - ``indices`` -- a list of ``rank + 1`` integers
        - ``prevent_multiple`` -- a boolean

        Return the circuit corresponding to the given indices.

        If ``prevent_multiple`` is true, a ``ValueError`` is raised if a multiple
        of this element has been computed before.

        .. NOTE::

            Raises a ``ValueError`` if the indices correspond to the zero vector.

        EXAMPLES::

            sage: from elementary_vectors import *
            sage: M = matrix([[1, 2, 4, 0], [0, 1, 2, 0]])
            sage: ce = CircuitEnumerator(M)

        In this example, circuits require 3 indices::

            sage: ce.circuit([0, 1, 2])
            (0, -2, 1, 0)
            sage: ce.circuit([1, 2, 3])
            Traceback (most recent call last):
            ...
            ValueError: The indices [1, 2, 3] correspond to the zero vector.
        """
        return self._element(indices, kernel=True, mark_zeros=prevent_multiple)

    def cocircuit(self, indices: List[int], prevent_multiple: bool = False) -> vector:
        r"""
        Compute the cocircuit corresponding to a list of indices.

        INPUT:

        - ``indices`` -- a list of ``rank - 1`` integers
        - ``prevent_multiple`` -- a boolean

        Return the cocircuit corresponding to the given indices.

        If ``prevent_multiple`` is true, a ``ValueError`` is raised if a multiple
        of this element has been computed before.

        .. NOTE::

            Raises a ``ValueError`` if the indices correspond to the zero vector.

        EXAMPLES::

            sage: from elementary_vectors import *
            sage: M = matrix([[1, 2, 4, 0], [0, 1, 2, 0]])
            sage: ce = CircuitEnumerator(M)

        In this example, cocircuits require 1 index::

            sage: ce.cocircuit([0])
            (0, -1, -2, 0)
            sage: ce.cocircuit([3])
            Traceback (most recent call last):
            ...
            ValueError: The indices [3] correspond to the zero vector.
        """
        return self._element(indices, kernel=False, mark_zeros=prevent_multiple)

    def _element(self, indices: List[int], kernel: bool, mark_zeros: bool = False) -> vector:
        r"""
        Compute the elementary vector corresponding to a list of indices.

        INPUT:

        - ``indices`` -- a list of ``rank - 1`` (for elements in the row space)
                         or ``rank + 1`` (for elements in the kernel) integers
        - ``kernel`` -- a boolean
        - ``prevent_multiple`` -- a boolean

        If ``kernel`` is true, return an elementary vector in the kernel
        and otherwise in the row space.

        If ``prevent_multiple`` is true, a ``ValueError`` is raised if a multiple
        of this element has been computed before.

        .. NOTE::

            Raises a ``ValueError`` if the indices correspond to the zero vector.

        EXAMPLES::

            sage: from elementary_vectors import *
            sage: M = matrix([[1, 2, 4, 0], [0, 1, 2, 0]])
            sage: ce = CircuitEnumerator(M)
            sage: ce._element_kernel([1, 2, 3], mark_zeros=False)
            (0, 0, 0, 0)
            sage: ce._element_row_space([3], mark_zeros=False)
            (0, 0, 0, 0)
        """
        if kernel:
            element = self._element_kernel(indices, mark_zeros=mark_zeros)
        else:
            element = self._element_row_space(indices, mark_zeros=mark_zeros)
        if self._is_element_zero(element):
            self._clear_zero_minors()
            raise ValueError(f"The indices {indices} correspond to the zero vector.")

        if mark_zeros and self._mark_zero_minors():
            raise MultipleException(f"Indices {indices} produce a nonzero multiple of a previously computed elementary vector.")
        return element

    def _element_kernel(self, indices: List[int], mark_zeros: bool) -> vector:
        element = self._zero_element()
        for pos in range(self.rank + 1):
            indices_minor = indices.copy()
            i = indices_minor.pop(pos)
            minor = self.minor(indices_minor, mark_if_zero=mark_zeros)
            if minor != 0:
                # check whether pos is even or odd
                self._set_element_entry(element, i, -minor if (pos & 1) else minor)
        return element

    def _element_row_space(self, indices: List[int], mark_zeros: bool) -> vector:
        element = self._zero_element()
        pos = 0
        for i in range(self.length):
            if i in indices:
                pos += 1
                continue
            minor = self.minor(sorted(indices + [i]), mark_if_zero=mark_zeros)
            if minor != 0:
                # check whether pos is even or odd
                self._set_element_entry(element, i, -minor if (pos & 1) else minor)
        return element

    def random_circuit(self) -> Optional[vector]:
        r"""
        Return a random circuit

        .. NOTE::

            If no circuit exists or the zero vector has been generated, ``None`` is returned.
        """
        try:
            return self.circuit(self._combinations_kernel.random_element())
        except ValueError: # no circuits exist or generated zero vector
            return

    def random_cocircuit(self) -> Optional[vector]:
        r"""
        Return a random cocircuit

        .. NOTE::

            If no cocircuit exists or the zero vector has been generated, ``None`` is returned.
        """
        try:
            return self.cocircuit(self._combinations_row_space.random_element())
        except ValueError: # no cocircuits exist or generated zero vector
            return

    def _zero_element(self) -> vector:
        return zero_vector(self.ring, self.length)

    @staticmethod
    def _set_element_entry(element: vector, index: int, value) -> None:
        element.set(index, value)

    @staticmethod
    def _is_element_zero(element: vector) -> bool:
        return element == 0

    def _mark_zero_minors(self) -> bool:
        r"""Return whether a marked minor is encountered again."""
        detected_marked_minor = not self._marked_minors.isdisjoint(self._zero_minors)
        self._marked_minors.update(self._zero_minors)
        self._clear_zero_minors()
        return detected_marked_minor

    def _clear_zero_minors(self) -> None:
        self._zero_minors.clear()

    def _index_sets_from_minor(self, indices_minor: List[int], kernel: bool) -> Iterator[List[int]]:
        r"""Generator of index sets corresponding to elementary vectors involving given minor."""
        if kernel:
            for i in range(self.length):
                if i in indices_minor:
                    continue
                yield sorted(indices_minor + [i])
        else:
            yield from Combinations(indices_minor, self.rank - 1)

    def circuit_generator(self, prevent_multiples: bool = True, reverse: bool = False) -> Iterator[vector]:
        r"""Return a generator of circuits"""
        return self._generator(kernel=True, prevent_multiples=prevent_multiples, reverse=reverse)

    def cocircuit_generator(self, prevent_multiples: bool = True, reverse: bool = False) -> Iterator[vector]:
        r"""Return a generator of cocircuits"""
        return self._generator(kernel=False, prevent_multiples=prevent_multiples, reverse=reverse)

    def _generator(self, kernel: bool, prevent_multiples: bool, reverse: bool) -> Iterator[vector]:
        r"""Return a generator of elementary vectors"""
        if prevent_multiples:
            self._reset_set_for_preventing_multiples()
        if kernel:
            combinations = self._combinations_kernel
        elif self.rank == 0:
            return
        else:
            combinations = self._combinations_row_space
        if reverse:
            combinations = reversed(combinations)
        for indices in combinations:
            try:
                yield self._element(indices, kernel=kernel, mark_zeros=prevent_multiples)
            except ValueError:
                pass

    def circuits(self, prevent_multiples: bool = True) -> List[vector]:
        r"""
        Return a list of circuits

        EXAMPLES::

            sage: from elementary_vectors import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: M
            [1 2 0 0]
            [0 1 2 3]
            sage: ce = CircuitEnumerator(M)
            sage: ce.circuits()
            [(4, -2, 1, 0), (6, -3, 0, 1), (0, 0, -3, 2)]

        TESTS::

            sage: ce = CircuitEnumerator(matrix(0, 4))
            sage: ce.circuits()
            [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
            sage: ce = CircuitEnumerator(matrix([[1, 0, 0], [0, 0, 0]]))
            sage: ce.circuits()
            [(0, -1, 0), (0, 0, -1)]
        """
        return list(self.circuit_generator(prevent_multiples=prevent_multiples))

    def cocircuits(self, prevent_multiples: bool = True) -> List[vector]:
        r"""
        Return a list of cocircuits

        EXAMPLES::

            sage: from elementary_vectors import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: M
            [1 2 0 0]
            [0 1 2 3]
            sage: ce = CircuitEnumerator(M)
            sage: ce.cocircuits()
            [(0, -1, -2, -3), (1, 0, -4, -6), (2, 4, 0, 0)]
        """
        return list(self.cocircuit_generator(prevent_multiples=prevent_multiples))

    def degenerate_circuits(self) -> List[vector]:
        r"""
        Return a list of degenerate circuits

        EXAMPLES::

            sage: from elementary_vectors import *
            sage: M = matrix([[1, 0, 1, 0], [0, 0, 1, 1]])
            sage: ce = CircuitEnumerator(M)
            sage: ce.degenerate_circuits()
            [(0, -1, 0, 0)]

        ::

            sage: M = matrix([[1, 1, 1, 0], [0, 1, 1, 1]])
            sage: ce = CircuitEnumerator(M)
            sage: ce.degenerate_circuits()
            [(0, -1, 1, 0)]

        We consider an example with 4 zero minors.
        There are six multiples that involve 2 of them each::

            sage: M = matrix([[1, -1, 0, 0, 1, 1], [0, 0, 1, 0, 1, 2], [0, 0, 0, 1, 1, 3]])
            sage: M
            [ 1 -1  0  0  1  1]
            [ 0  0  1  0  1  2]
            [ 0  0  0  1  1  3]
            sage: ce = CircuitEnumerator(M)
            sage: ce.degenerate_circuits()
            [(-1, -1, 0, 0, 0, 0)]
        """
        return list(self._degenerate_elements(kernel=True))

    def degenerate_cocircuits(self) -> List[vector]:
        r"""
        Return a list of degenerate cocircuits

        EXAMPLES::

            sage: from elementary_vectors import *
            sage: M = matrix([[1, 0, 1, 0], [0, 0, 1, 1]])
            sage: M
            [1 0 1 0]
            [0 0 1 1]
            sage: ce = CircuitEnumerator(M)
            sage: ce.degenerate_cocircuits()
            [(0, 0, -1, -1), (1, 0, 0, -1), (1, 0, 1, 0)]
        """
        return list(self._degenerate_elements(kernel=False))

    def _degenerate_elements(self, kernel: bool) -> Iterator[vector]:
        r"""
        Generator of elementary vectors with smaller-than-usual support.
        """
        self._reset_set_for_preventing_multiples()
        for indices_minor in Combinations(self.length, self.rank):
            if self.minor(indices_minor) != 0:
                continue
            if tuple(indices_minor) in self._marked_minors:
                continue
            for indices in self._index_sets_from_minor(indices_minor, kernel=kernel):
                try:
                    if kernel:
                        yield self.circuit(indices, prevent_multiple=True)
                    else:
                        yield self.cocircuit(indices, prevent_multiple=True)
                    break
                except MultipleException:
                    break
                except ValueError: # zero vector
                    continue


class MultipleException(ValueError):
    r"""Raised when a multiple of a previously computed elementary vector is detected."""
