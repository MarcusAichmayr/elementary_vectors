r"""
Computing elementary vectors

Elementary vectors are support minimal vectors of a subspace.
They are also called _circuits_.
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

from typing import Optional, Union, List, Iterator

from sage.combinat.combination import Combinations
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import zero_vector, vector
from sage.structure.sage_object import SageObject

from .utility import is_symbolic


def circuits(matrix, prevent_multiples: bool = True, generator: bool = False) -> Union[List[vector], Iterator[vector]]:
    r"""
    Compute the circuits of a matrix.

    OUTPUT:
    Return a list of circuits of the matrix.
    These are the nonzero support-minimal elements in the kernel.

    .. SEEALSO::

        - - :func:`cocircuits`
        - - :func:`elementary_vectors`

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
        sage: M
        [1 2 0 0]
        [0 1 2 3]
        sage: circuits(M)
        [(4, -2, 1, 0), (6, -3, 0, 1), (0, 0, -3, 2)]
    """
    return elementary_vectors(matrix, kernel=True, prevent_multiples=prevent_multiples, generator=generator)


def cocircuits(matrix: Matrix, prevent_multiples: bool = True, generator: bool = False) -> Union[List[vector], Iterator[vector]]:
    r"""
    Compute the cocircuits of a matrix.

    OUTPUT:
    Return a list of cocircuits of the matrix.
    These are the nonzero support-minimal elements in the row space.

    .. SEEALSO::

        - - :func:`circuits`
        - - :func:`elementary_vectors`

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
        sage: M
        [1 2 0 0]
        [0 1 2 3]
        sage: cocircuits(M)
        [(0, -1, -2, -3), (1, 0, -4, -6), (2, 4, 0, 0)]
    """
    return elementary_vectors(matrix, kernel=False, prevent_multiples=prevent_multiples, generator=generator)


def elementary_vectors(matrix, kernel: bool = True, prevent_multiples: bool = True, generator: bool = False) -> Union[List[vector], Iterator[vector]]:
    r"""
    Compute elementary vectors of a subspace determined by a matrix.

    INPUT:

    - ``matrix`` -- a matrix
    - ``kernel`` -- a boolean (default: ``True``)
    - ``prevent_multiples`` -- a boolean (default: ``True``)
    - ``generator`` -- a boolean (default: ``False``)

    OUTPUT:
    Return a list of elementary vectors in the kernel of a given matrix.
    To compute the elementary vectors in the row space, pass ``False`` for ``kernel``.

    - If ``generator`` is ``True``, the output will be a generator object instead of a list.

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
        sage: M
        [1 2 0 0]
        [0 1 2 3]
        sage: elementary_vectors(M)
        [(4, -2, 1, 0), (6, -3, 0, 1), (0, 0, -3, 2)]
        sage: elementary_vectors(M, prevent_multiples=False)
        [(4, -2, 1, 0), (6, -3, 0, 1), (0, 0, -3, 2), (0, 0, -6, 4)]

    By default, the elementary vectors in the kernel are computed.
    To consider the row space, pass ``kernel=False``::

        sage: elementary_vectors(M, kernel=False)
        [(0, -1, -2, -3), (1, 0, -4, -6), (2, 4, 0, 0)]

    We can also compute elementary vectors over a finite field::

        sage: B = matrix(GF(7), [[1, 2, 3, 4, 0], [0, 5, 2, 3, 3]])
        sage: B
        [1 2 3 4 0]
        [0 5 2 3 3]
        sage: elementary_vectors(B)
        [(3, 5, 5, 0, 0),
         (0, 4, 0, 5, 0),
         (6, 4, 0, 0, 5),
         (1, 0, 4, 2, 0),
         (2, 0, 4, 0, 2),
         (5, 0, 0, 4, 3),
         (0, 2, 1, 0, 3),
         (0, 0, 5, 5, 1)]

    Variables are also supported::

        sage: var('a, b')
        (a, b)
        sage: M = matrix([[1, 2, a, 0], [0, 1, 2, b]])
        sage: M
        [1 2 a 0]
        [0 1 2 b]
        sage: elementary_vectors(M)
        [(-a + 4, -2, 1, 0), (2*b, -b, 0, 1), (a*b, 0, -b, 2), (0, a*b, -2*b, -a + 4)]

    Matrices over the polynomial ring work, too::

        sage: R = PolynomialRing(ZZ, "x")
        sage: x = R.gen()
        sage: M = matrix([[1, 2, x, 0], [0, 1, 2, x]])
        sage: M
        [1 2 x 0]
        [0 1 2 x]
        sage: elementary_vectors(M)
        [(-x + 4, -2, 1, 0), (2*x, -x, 0, 1), (x^2, 0, -x, 2), (0, x^2, -2*x, -x + 4)]
        sage: R = PolynomialRing(ZZ, "x, y")
        sage: x, y = R.gens()
        sage: M = matrix([[x, y, 0, 0], [0, 1, 2, 3]])
        sage: M
        [x y 0 0]
        [0 1 2 3]
        sage: elementary_vectors(M)
        [(2*y, -2*x, x, 0), (3*y, -3*x, 0, x), (0, 0, -3*x, 2*x)]

    TESTS::

        sage: elementary_vectors(random_matrix(QQ, 0, 4))
        [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
        sage: elementary_vectors(matrix([[1, 0, 0], [0, 0, 0]]))
        [(0, -1, 0), (0, 0, -1)]
    """
    if generator:
        return ElementaryVectors(matrix).generator(kernel=kernel, prevent_multiples=prevent_multiples)
    return ElementaryVectors(matrix).elements(kernel=kernel, prevent_multiples=prevent_multiples)


def kernel_matrix_using_elementary_vectors(matrix: Matrix) -> Matrix:
    """
    Division-free right kernel matrix based on elementary vectors.

    OUTPUT:
    A right kernel matrix.
    It also works for symbolic matrices.

    .. NOTE::

        Raises a ``ValueError`` if the matrix has no constant nonzero maximal minor.

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[1, 0, 1, -1, 0], [0, 1, 1, 1, -1]])
        sage: M
        [ 1  0  1 -1  0]
        [ 0  1  1  1 -1]
        sage: kernel_matrix_using_elementary_vectors(M)
        [1 0 0 1 1]
        [0 1 0 0 1]
        [0 0 1 1 2]
        sage: M = matrix([[1, 1, -1, -1, 0], [2, 1, -1, 0, 1], [1, 1, 1, 1, 1]])
        sage: M
        [ 1  1 -1 -1  0]
        [ 2  1 -1  0  1]
        [ 1  1  1  1  1]
        sage: kernel_matrix_using_elementary_vectors(M)
        [-1  0  0 -1  2]
        [ 0 -1  1 -2  2]
        sage: var('a')
        a
        sage: M = matrix([[1, 0, 1, -1, 0], [0, 1, a, 1, -1]])
        sage: M
        [ 1  0  1 -1  0]
        [ 0  1  a  1 -1]
        sage: kernel_matrix_using_elementary_vectors(M)
        [    1     0     0     1     1]
        [    0     1     0     0     1]
        [    0     0     1     1 a + 1]

    TESTS::

        sage: kernel_matrix_using_elementary_vectors(identity_matrix(3, 3))
        []
        sage: _.dimensions()
        (0, 3)
        sage: kernel_matrix_using_elementary_vectors(matrix(0, 3))
        [1 0 0]
        [0 1 0]
        [0 0 1]

    ::

        sage: M = matrix([[0, 1], [0, 1]])
        sage: kernel_matrix_using_elementary_vectors(M)
        [1 0]
    """
    evs = ElementaryVectors(matrix)

    rank = evs.rank
    length = evs.length

    if rank == length:
        return Matrix(matrix.base_ring(), 0, length)

    for indices_minor in Combinations(range(length - 1, -1, -1), rank):
        minor = evs.minor(indices_minor)
        if minor != 0 and not is_symbolic(minor):
            return Matrix(evs.element(indices) for indices in evs._index_sets_from_minor(indices_minor, kernel=True))
    raise ValueError("Matrix has no constant nonzero maximal minor.")


class ElementaryVectors(SageObject):
    r"""
    A class used to compute elementary vectors.

    Whenever a maximal minor is computed, it is stored in a dictionary for efficient reuse.
    Supports elementary vectors in the kernel and row space.

    EXAMPLES::

        sage: from elementary_vectors.functions import ElementaryVectors
        sage: M = matrix([[1, 2, 4, 1, -1], [0, 1, 2, 3, 4]])
        sage: evs = ElementaryVectors(M)
        sage: evs.elements()
        [(0, -2, 1, 0, 0),
         (5, -3, 0, 1, 0),
         (9, -4, 0, 0, 1),
         (10, 0, -3, 2, 0),
         (18, 0, -4, 0, 2),
         (7, 0, 0, -4, 3),
         (0, 7, 0, -9, 5),
         (0, 0, 7, -18, 10)]
        sage: evs.elements(prevent_multiples=False)
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
        sage: evs.elements(kernel=False)
        [(0, -1, -2, -3, -4), (1, 0, 0, -5, -9), (3, 5, 10, 0, -7), (4, 9, 18, 7, 0)]
        sage: evs.elements(kernel=False, prevent_multiples=False)
        [(0, -1, -2, -3, -4),
         (1, 0, 0, -5, -9),
         (2, 0, 0, -10, -18),
         (3, 5, 10, 0, -7),
         (4, 9, 18, 7, 0)]

    ::

        sage: evs.minor([0, 2])
        2
        sage: evs.element([0, 2, 3])
        (10, 0, -3, 2, 0)
        sage: evs.element([0])
        (0, -1, -2, -3, -4)
        sage: evs.element([0, 2, 3], kernel=True)
        (10, 0, -3, 2, 0)
        sage: evs.element([0], kernel=False)
        (0, -1, -2, -3, -4)
        sage: evs.random_element() # random
        (0, 0, 7, -18, 10)
        sage: evs.random_element(kernel=False) # random
        (3, 5, 10, 0, -7)

    We consider an example that involves many zero minors::

        sage: M = matrix([[1, 2, 4, 0], [0, 1, 2, 0]])
        sage: M.minors(2)
        [1, 2, 0, 0, 0, 0]
        sage: evs = ElementaryVectors(M)
        sage: evs.elements()
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

        The minor is cached for efficient reuse.

        TESTS::

            sage: from elementary_vectors.functions import ElementaryVectors
            sage: M = matrix([[1, 2, 4, 1, -1], [0, 1, 2, 3, 4]])
            sage: evs = ElementaryVectors(M)
            sage: evs._minors
            {}
            sage: evs.minor([0, 1])
            1
            sage: evs._minors
            {(0, 1): 1}
            sage: evs.minor([2, 4])
            18
            sage: evs._minors
            {(0, 1): 1, (2, 4): 18}
            sage: evs.minor([0, 1, 2])
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

    def element(self, indices: List[int], kernel: Optional[bool] = None, prevent_multiple: bool = False) -> vector:
        r"""
        Compute the elementary vector corresponding to a list of indices.

        INPUT:

        - ``indices`` -- a list of ``rank - 1`` (for elements in the row space)
                         or ``rank + 1`` (for elements in the kernel) integers
        - ``kernel`` -- a boolean
        - ``prevent_multiple`` -- a boolean

        If ``kernel`` is true, return an elementary vector in the kernel
        and otherwise in the row space.
        If not specified, this is determined from the number of indices.

        If ``prevent_multiple`` is true, a ``ValueError`` is raised if a multiple
        of this element has been computed before.

        .. NOTE::

            Raises a ``ValueError`` if the indices correspond to the zero vector.

        EXAMPLES::

            sage: from elementary_vectors import *
            sage: M = matrix([[1, 2, 4, 0], [0, 1, 2, 0]])
            sage: evs = ElementaryVectors(M)

        Elementary vectors in the kernel require 3 indices::

            sage: evs.element([0, 1, 2])
            (0, -2, 1, 0)
            sage: evs.element([1, 2, 3])
            Traceback (most recent call last):
            ...
            ValueError: The indices [1, 2, 3] correspond to the zero vector!

        For the row space, we need 1 element::

            sage: evs.element([0])
            (0, -1, -2, 0)

        TESTS::

            sage: evs.element([1, 2])
            Traceback (most recent call last):
            ...
            ValueError: The number of indices should be 1 or 3 but got 2.

        ::

            evs._element_kernel([1, 2, 3])
            (0, 0, 0, 0)
            evs._element_row_space([3])
            (0, 0, 0, 0)
        """
        if kernel is None:
            if len(indices) == self.rank + 1:
                kernel = True
            elif len(indices) == self.rank - 1:
                kernel = False
            else:
                raise ValueError(f"The number of indices should be {self.rank - 1} or {self.rank + 1} but got {len(indices)}.")

        if kernel:
            element = self._element_kernel(indices, mark_zeros=prevent_multiple)
        else:
            element = self._element_row_space(indices, mark_zeros=prevent_multiple)
        if element == 0:
            self._clear_zero_minors()
            raise ValueError(f"The indices {indices} correspond to the zero vector!")

        if prevent_multiple and self._mark_zero_minors():
            raise MultipleException(f"Indices {indices} produce a nonzero multiple of a previously computed elementary vector!")
        return element

    def _element_kernel(self, indices: List[int], mark_zeros: bool) -> vector:
        element = self._zero_element()
        for pos in range(self.rank + 1):
            indices_minor = indices.copy()
            i = indices_minor.pop(pos)
            minor = self.minor(indices_minor, mark_if_zero=mark_zeros)
            if minor != 0:
                # check whether pos is even or odd
                element[i] = -minor if (pos & 1) else minor
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
                element[i] = -minor if (pos & 1) else minor
        return element

    def random_element(self, kernel: bool = True) -> Optional[vector]:
        r"""
        Return a random elementary vector

        .. NOTE::

            If no elementary vector exists or the zero vector has been generated, ``None`` is returned.
        """
        try:
            if kernel:
                return self.element(self._combinations_kernel.random_element(), kernel=kernel)
            return self.element(self._combinations_row_space.random_element(), kernel=kernel)
        except ValueError: # no elementary vectors exist or generated zero vector
            return

    def _zero_element(self) -> tuple:
        return zero_vector(self.ring, self.length)

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

    def generator(self, kernel: bool = True, prevent_multiples: bool = True, reverse: bool = False) -> Iterator[vector]:
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
                yield self.element(indices, kernel=kernel, prevent_multiple=prevent_multiples)
            except ValueError:
                pass

    def elements(self, kernel: bool = True, prevent_multiples: bool = True) -> List[vector]:
        r"""Return a list of elementary vectors"""
        return list(self.generator(kernel=kernel, prevent_multiples=prevent_multiples))

    def degenerate_elements(self, kernel: bool = True) -> Iterator[vector]:
        r"""
        Generator of elementary vectors with smaller-than-usual support.

        EXAMPLES::

            sage: from elementary_vectors import *
            sage: M = matrix([[1, 0, 1, 0], [0, 0, 1, 1]])
            sage: evs = ElementaryVectors(M)
            sage: list(evs.degenerate_elements())
            [(0, -1, 0, 0)]
            sage: list(evs.degenerate_elements(kernel=False))
            [(0, 0, -1, -1), (1, 0, 0, -1), (1, 0, 1, 0)]

        ::

            sage: M = matrix([[1, 1, 1, 0], [0, 1, 1, 1]])
            sage: evs = ElementaryVectors(M)
            sage: list(evs.degenerate_elements())
            [(0, -1, 1, 0)]
            sage: list(evs.degenerate_elements(kernel=False))
            [(1, 0, 0, -1)]

        We consider an example with 4 zero minors.
        There are six multiples that involve 2 of them each::

            sage: M = matrix([[1, -1, 0, 0, 1, 1], [0, 0, 1, 0, 1, 2], [0, 0, 0, 1, 1, 3]])
            sage: M
            [ 1 -1  0  0  1  1]
            [ 0  0  1  0  1  2]
            [ 0  0  0  1  1  3]
            sage: evs = ElementaryVectors(M)
            sage: list(evs.degenerate_elements())
            [(-1, -1, 0, 0, 0, 0)]
        """
        self._reset_set_for_preventing_multiples()
        for indices_minor in Combinations(self.length, self.rank):
            if self.minor(indices_minor) != 0:
                continue
            if tuple(indices_minor) in self._marked_minors:
                continue
            for indices in self._index_sets_from_minor(indices_minor, kernel=kernel):
                try:
                    yield self.element(indices, kernel=kernel, prevent_multiple=True)
                    break
                except MultipleException:
                    break
                except ValueError: # zero vector
                    continue


class MultipleException(ValueError):
    r"""Raised when a multiple of a previously computed elementary vector is detected."""
