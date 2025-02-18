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

import warnings
from collections.abc import Generator
from sage.combinat.combination import Combinations
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import zero_vector
from sage.structure.sage_object import SageObject

from .utility import is_symbolic


def elementary_vectors(M, dual: bool = True, generator: bool = False):
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
        return ElementaryVectors(M).generator(dual=dual)
    return ElementaryVectors(M).elements(dual=dual)


def kernel_matrix_using_elementary_vectors(M):
    """
    Division-free computation of the right kernel based on elementary vectors.

    INPUT:

    - ``M`` -- a matrix

    OUTPUT:
    A right kernel matrix of ``M``.
    Also works for symbolic matrices.

    .. NOTE::

        A ``ValueError`` is returned if ``M`` has no constant nonzero maximal minor.

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
    constant_minor_found = False
    for indices in Combinations(reversed(range(length)), rank):
        minor = evs.minor(indices)
        if minor and not is_symbolic(minor):
            constant_minor_found = True
            break
    if not constant_minor_found:
        raise ValueError("Could not find a constant nonzero maximal minor.")

    return matrix(evs.element(set(indices).union([k])) for k in range(length) if k not in indices)


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

        sage: evs = ElementaryVectors(M)
        sage: evs.minor([0, 2])
        2
        sage: evs.element([0, 2, 3])
        (10, 0, -3, 2, 0)
        sage: evs.element([0])
        (0, -1, -2, -3, -4)
        sage: evs.element_kernel([0, 2, 3])
        (10, 0, -3, 2, 0)
        sage: evs.element_row_space([0])
        (0, -1, -2, -3, -4)
        sage: evs.random_element() # random
        (0, 0, 7, -18, 10)
        sage: evs.random_element(dual=False) # random
        (3, 5, 10, 0, -7)
    """
    def __init__(self, M) -> None:
        try:
            self.matrix = M.matrix_from_rows(M.pivot_rows())  # fails for polynomial matrices
        except (ArithmeticError, NotImplementedError):
            self.matrix = M
            warnings.warn("Could not determine rank of matrix. Expect wrong result!")

        self.length = self.matrix.ncols()
        self.rank = self.matrix.nrows()
        self.ring = self.matrix.base_ring()
        self.minors = {}
        self.set_combinations()
        self.set_combinations_dual()

        self._reset_set_for_preventing_multiples()

    def minor(self, indices):
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
        """
        indices = tuple(indices)
        try:
            return self.minors[indices]
        except KeyError:
            self.minors[indices] = self.matrix.matrix_from_columns(indices).det()
            return self.minors[indices]

    def set_combinations(self, combinations=None) -> None:
        r"""Set or reset combinations."""
        if combinations is None:
            self._combinations = Combinations(self.length, self.rank + 1)
        else:
            self._combinations = combinations

    def set_combinations_dual(self, combinations=None) -> None:
        r"""Set or reset combinations."""
        if combinations is None:
            self._combinations_dual = Combinations(self.length, self.rank - 1)
        else:
            self._combinations_dual = combinations

    def element(self, indices: list, dual: bool = None, prevent_multiple: bool = False):
        r"""
        Compute the elementary vector corresponding to a list of indices.

        .. NOTE::

            Raises a ``ValueError`` if the indices correspond to the zero vector.
        """
        if dual is None:
            if len(indices) == self.rank + 1:
                dual = True
            elif len(indices) == self.rank - 1:
                dual = False
            else:
                raise ValueError("Number of indices does not fit!")
        if dual:
            return self.element_kernel(indices, prevent_multiple=prevent_multiple)
        return self.element_row_space(indices, prevent_multiple=prevent_multiple)

    def element_kernel(self, indices: list, prevent_multiple: bool = False):
        r"""
        Compute the elementary vector corresponding to a list of indices.

        .. NOTE::

            Raises a ``ValueError`` if the indices correspond to the zero vector.
        """
        if prevent_multiple:
            return self._element_kernel_prevent_multiple(indices)
        element = self._zero_element()
        nonzero_detected = False
        for pos, i in enumerate(indices):
            minor = self.minor(tuple(j for j in indices if j != i))
            if minor == 0:
                continue
            nonzero_detected = True
            element[i] = (-1) ** pos * minor
        if nonzero_detected:
            return element
        raise ValueError("Indices correspond to zero vector!")

    def element_row_space(self, indices: list, prevent_multiple: bool = False):
        """
        Compute the elementary vector in the row space corresponding to the given indices.

        INPUT::

        - ``indices`` -- a list of ``rank - 1``elements
        """
        if prevent_multiple:
            return self._element_row_space_prevent_multiple(indices)
        element = self._zero_element()
        nonzero_detected = False
        pos = 0
        for i in range(self.length):
            if i in indices:
                pos += 1
                continue
            minor = self.minor(tuple(set(indices + [i])))
            if minor == 0:
                continue
            nonzero_detected = True
            element[i] = (-1) ** pos * minor
        if nonzero_detected:
            return element
        raise ValueError("Indices correspond to zero vector!")

    def _element_kernel_prevent_multiple(self, indices: list):
        r"""
        Compute the elementary vector corresponding to a list of indices.

        .. NOTE::

            If this results in a multiple of a previous element, a ``ValueError`` is raised.
        """
        element = self._zero_element()
        nonzero_detected = False
        zero_minors = []
        multiple_detected = False
        for pos, i in enumerate(indices):
            indices_minor = tuple(j for j in indices if j != i)
            if indices_minor in self.marked_minors:
                multiple_detected = True
                continue
            minor = self.minor(indices_minor)
            if minor == 0:
                zero_minors.append(indices_minor)
                continue
            nonzero_detected = True
            element[i] = (-1) ** pos * minor
        if nonzero_detected:
            for marked_minor in zero_minors:
                self.marked_minors.add(marked_minor)
            if multiple_detected:
                raise ValueError("Multiple detected!")
            return element
        raise ValueError("Indices correspond to zero vector!")

    def _element_row_space_prevent_multiple(self, indices: list):
        r"""
        Compute the elementary vector corresponding to a list of indices.

        .. NOTE::

            If this results in a multiple of a previous element, a ``ValueError`` is raised.
        """
        element = self._zero_element()
        pos = 0
        nonzero_detected = False
        zero_minors = []
        multiple_detected = False
        for k in range(self.length):
            if k in indices:
                pos += 1
                continue
            indices_minor = tuple(set(indices + [k]))
            if indices_minor in self.marked_minors:
                multiple_detected = True
                continue
            minor = self.minor(indices_minor)
            if minor == 0:
                zero_minors.append(indices_minor)
                continue
            nonzero_detected = True
            element[k] = (-1) ** pos * minor
        if nonzero_detected:
            for marked_minor in zero_minors:
                self.marked_minors.add(marked_minor)
            if multiple_detected:
                raise ValueError("Multiple detected!")
            return element
        raise ValueError("Indices correspond to zero vector!")

    def random_element(self, dual: bool = True):
        r"""
        Return a random elementary vector

        .. NOTE::

            If no elementary vector exists or the zero vector has been generated, ``None`` is returned.
        """
        try:
            if dual:
                return self.element_kernel(self._combinations.random_element())
            return self.element_row_space(self._combinations_dual.random_element())
        except ValueError: # no elementary vectors exist or generated zero vector
            return

    def _zero_element(self):
        return zero_vector(self.ring, self.length)

    def _reset_set_for_preventing_multiples(self) -> None:
        self.marked_minors = set()

    def index_sets_from_minor(self, indices: list, dual: bool = True) -> Generator[list]:
        r"""Generator of index sets corresponding to elementary vectors involving given minor."""
        if dual:
            for i in range(self.length):
                if i in indices:
                    continue
                yield list(set(indices + [i]))
        else:
            for i in indices:
                yield [j for j in indices if j != i]

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
            if self.minor(indices_minor):
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
            combinations = self._combinations
        else:
            combinations = self._combinations_dual
            if self.rank == 0:
                return
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
