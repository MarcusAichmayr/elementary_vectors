r"""Computing elementary vectors"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from collections.abc import Generator
from sage.combinat.combination import Combinations
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import zero_vector
from sage.structure.sage_object import SageObject

from .utility import is_symbolic


def elementary_vectors(M, dual: bool = True, prevent_multiples: bool = True, generator: bool = False):
    r"""
    Compute elementary vectors of a subspace determined by a matrix or a list of maximal minors.

    INPUT:

    - ``M`` -- a matrix

    - ``dual`` -- a boolean (default: ``True``)

    - ``generator`` -- a boolean (default: ``False``)

    OUTPUT:
    Return a list of elementary vectors in the kernel of a matrix given by the matrix ``M``.
    To compute the elementary vectors in the row space, pass false for ``dual``.

    - If ``generator`` is true, the output will be a generator object instead of a list.

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
    To consider the row space, pass ``dual=False``::

        sage: elementary_vectors(M, dual=False)
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

    TESTS::

        sage: elementary_vectors(random_matrix(QQ, 0, 4))
        [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
        sage: elementary_vectors(matrix([[1, 0, 0], [0, 0, 0]]))
        [(0, -1, 0), (0, 0, -1)]
    """
    if generator:
        return ElementaryVectors(M).generator(dual=dual, prevent_multiples=prevent_multiples)
    return ElementaryVectors(M).elements(dual=dual, prevent_multiples=prevent_multiples)


def kernel_matrix_using_elementary_vectors(M):
    """
    Division-free computation of the right kernel based on elementary vectors.

    INPUT:

    - ``M`` -- a matrix

    OUTPUT:
    A right kernel matrix of ``M``.
    Also works for symbolic matrices.

    .. NOTE::

        Raises a ``ValueError`` if ``M`` has no constant nonzero maximal minor.

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
    """
    rank, length = M.dimensions()

    if rank == length:
        return matrix(M.base_ring(), 0, length)

    evs = ElementaryVectors(M)
    for indices_minor in Combinations(range(length - 1, -1, -1), rank):
        minor = evs.minor(indices_minor)
        if minor and not is_symbolic(minor):
            return matrix(evs.element(indices) for indices in evs.index_sets_from_minor(indices_minor))
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
        sage: evs.elements(dual=False)
        [(0, -1, -2, -3, -4), (1, 0, 0, -5, -9), (3, 5, 10, 0, -7), (4, 9, 18, 7, 0)]
        sage: evs.elements(dual=False, prevent_multiples=False)
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
        sage: evs.element([0, 2, 3], dual=True)
        (10, 0, -3, 2, 0)
        sage: evs.element([0], dual=False)
        (0, -1, -2, -3, -4)
        sage: evs.random_element() # random
        (0, 0, 7, -18, 10)
        sage: evs.random_element(dual=False) # random
        (3, 5, 10, 0, -7)

    We consider an example that involves many zero minors::

        sage: M = matrix([[1, 2, 4, 0], [0, 1, 2, 0]])
        sage: M.minors(2)
        [1, 2, 0, 0, 0, 0]
        sage: evs = ElementaryVectors(M)
        sage: evs.elements()
        [(0, -2, 1, 0), (0, 0, 0, 1)]
    """
    def __init__(self, M) -> None:
        self.matrix = M.matrix_from_rows(M.pivot_rows())
        self.rank, self.length = self.matrix.dimensions()
        self.ring = M.base_ring()
        self.minors = {}
        self._zero_minors = set()

        self.set_combinations_kernel()
        self.set_combinations_row_space()
        self._reset_set_for_preventing_multiples()

    def minor(self, indices: list[int], mark_if_zero: bool = False):
        r"""
        Compute a minor given by indices

        The minor is cached for efficient reuse.

        TESTS::

            sage: from elementary_vectors.functions import ElementaryVectors
            sage: M = matrix([[1, 2, 4, 1, -1], [0, 1, 2, 3, 4]])
            sage: evs = ElementaryVectors(M)
            sage: evs.minors
            {}
            sage: evs.minor([0, 1])
            1
            sage: evs.minors
            {(0, 1): 1}
            sage: evs.minor([2, 4])
            18
            sage: evs.minors
            {(0, 1): 1, (2, 4): 18}
            sage: evs.minor([0, 1, 2])
            Traceback (most recent call last):
            ...
            ValueError: Indices (0, 1, 2) should have size 2 and not 3.
        """
        indices = tuple(indices)
        minor = self.minors.get(indices)
        if minor is None:
            try:
                minor = self.matrix.matrix_from_columns(indices).det()
                self.minors[indices] = minor
            except ValueError as e:
                raise ValueError(f"Indices {indices} should have size {self.rank} and not {len(indices)}.") from e
        if mark_if_zero and minor == 0:
            self._zero_minors.add(indices)
        return minor

    def set_combinations_kernel(self, combinations=None) -> None:
        r"""Set or reset combinations for elements in the kernel."""
        if combinations is None:
            self._combinations_kernel = Combinations(self.length, self.rank + 1)
        else:
            self._combinations_kernel = combinations

    def set_combinations_row_space(self, combinations=None) -> None:
        r"""Set or reset combinations for elements in the row space."""
        if combinations is None:
            self._combinations_row_space = Combinations(self.length, self.rank - 1)
        else:
            self._combinations_row_space = combinations

    def element(self, indices: list[int], dual: bool = None, prevent_multiple: bool = False):
        r"""
        Compute the elementary vector corresponding to a list of indices.

        INPUT:

        - ``indices`` -- a list of ``rank - 1`` (for elements in the row space)
                         or ``rank + 1`` (for elements in the kernel) integers

        - ``dual`` -- a boolean

        - ``prevent_multiple`` -- a boolean

        If ``dual`` is true, return an elementary vector in the kernel
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
            ValueError: Indices [1, 2, 3] correspond to zero vector!

        For the row space, we need 1 element::

            sage: evs.element([0])
            (0, -1, -2, 0)

        TESTS::

            sage: evs.element([1, 2])
            Traceback (most recent call last):
            ...
            ValueError: Number of indices should be 1 or 3 but got 2.

        ::

            evs._element_kernel([1, 2, 3])
            (0, 0, 0, 0)
            evs._element_row_space([3])
            (0, 0, 0, 0)
        """
        if dual is None:
            if len(indices) == self.rank + 1:
                dual = True
            elif len(indices) == self.rank - 1:
                dual = False
            else:
                raise ValueError(f"Number of indices should be {self.rank - 1} or {self.rank + 1} but got {len(indices)}.")

        if dual:
            element = self._element_kernel(indices, mark_zeros=prevent_multiple)
        else:
            element = self._element_row_space(indices, mark_zeros=prevent_multiple)
        if element == 0:
            self._zero_minors.clear()
            raise ValueError(f"Indices {indices} correspond to zero vector!")

        if prevent_multiple:
            multiple = False
            while self._zero_minors:
                minor = self._zero_minors.pop()
                if minor in self._marked_minors:
                    multiple = True
                    continue
                self._marked_minors.add(minor)

            if multiple:
                raise ValueError(f"Indices {indices} produce a multiple of a previously computed element!")
        return element

    def _element_kernel(self, indices: list[int], mark_zeros: bool = False):
        element = self._zero_element()
        for pos in range(self.rank + 1):
            indices_minor = indices.copy()
            i = indices_minor.pop(pos)
            minor = self.minor(indices_minor, mark_if_zero=mark_zeros)
            if minor == 0:
                continue
            element[i] = (-1) ** pos * minor
        return element

    def _element_row_space(self, indices: list[int], mark_zeros: bool = False):
        element = self._zero_element()
        pos = 0
        for i in range(self.length):
            if i in indices:
                pos += 1
                continue
            minor = self.minor(sorted(indices + [i]), mark_if_zero=mark_zeros)
            if minor == 0:
                continue
            element[i] = (-1) ** pos * minor
        return element

    def random_element(self, dual: bool = True):
        r"""
        Return a random elementary vector

        .. NOTE::

            If no elementary vector exists or the zero vector has been generated, ``None`` is returned.
        """
        try:
            return self.element(
                (self._combinations_kernel if dual else self._combinations_row_space).random_element(),
                dual=dual)
        except ValueError: # no elementary vectors exist or generated zero vector
            return

    def _zero_element(self):
        return zero_vector(self.ring, self.length)

    def _reset_set_for_preventing_multiples(self) -> None:
        self._marked_minors = set()

    def index_sets_from_minor(self, indices_minor: list[int], dual: bool = True) -> Generator[list]:
        r"""Generator of index sets corresponding to elementary vectors involving given minor."""
        if dual:
            for i in range(self.length):
                if i in indices_minor:
                    continue
                yield sorted(indices_minor + [i])
        else:
            yield from Combinations(indices_minor, self.rank - 1)

    def elements_with_smaller_support(self, dual: bool = True) -> Generator:
        r"""
        Generator of elementary vectors with smaller than usual support.

        EXAMPLES::

            sage: from elementary_vectors import *
            sage: M = matrix([[1, 0, 1, 0], [0, 0, 1, 1]])
            sage: evs = ElementaryVectors(M)
            sage: list(evs.elements_with_smaller_support())
            [(0, -1, 0, 0)]
            sage: list(evs.elements_with_smaller_support(dual=False))
            [(0, 0, -1, -1), (1, 0, 0, -1), (1, 0, 1, 0)]

        ::

            sage: M = matrix([[1, 1, 1, 0], [0, 1, 1, 1]])
            sage: evs = ElementaryVectors(M)
            sage: list(evs.elements_with_smaller_support())
            [(0, -1, 1, 0)]
            sage: list(evs.elements_with_smaller_support(dual=False))
            [(1, 0, 0, -1)]
        """
        self._reset_set_for_preventing_multiples()
        for indices_minor in Combinations(self.length, self.rank):
            if self.minor(indices_minor) != 0:
                continue
            for indices in self.index_sets_from_minor(indices_minor, dual=dual):
                try:
                    yield self.element(indices, dual=dual, prevent_multiple=True)
                except ValueError:
                    pass

    def generator(
        self,
        dual: bool = True,
        prevent_multiples: bool = True,
        reverse: bool = False
    ) -> Generator:
        r"""Return a generator of elementary vectors"""
        if prevent_multiples:
            self._reset_set_for_preventing_multiples()
        if dual:
            combinations = self._combinations_kernel
        else:
            if self.rank == 0:
                return
            combinations = self._combinations_row_space
        if reverse:
            combinations = reversed(combinations)
        for indices in combinations:
            try:
                yield self.element(indices, dual=dual, prevent_multiple=prevent_multiples)
            except ValueError:
                pass

    def elements(self, dual: bool = True, prevent_multiples: bool = True) -> list:
        r"""Return a list of elementary vectors"""
        return list(self.generator(dual=dual, prevent_multiples=prevent_multiples))
