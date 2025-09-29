"""Utility functions"""

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
from typing import Iterator

from sage.combinat.combination import Combinations
from sage.functions.generalized import sign
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector, zero_vector
from sage.structure.sage_object import SageObject

from sign_vectors import sign_vector, SignVector
from elementary_vectors.functions import ElementaryVectors


def vector_from_sign_vector(data, sv: SignVector) -> vector:
    r"""
    Find a vector in the row space of a matrix that has given signs.

    INPUT:

    - ``data`` -- either a real matrix with ``n`` columns or a list of
                elementary vectors of length ``n``
    - ``sv`` -- a sign vector of length ``n``

    OUTPUT:
    Return a conformal sum of elementary vectors that lies in the given subspace.

    If ``data`` is a matrix, the elementary vectors in the kernel of this matrix are used for the result.
    If ``data`` is a list of elementary vectors, those are used.

    .. NOTE::

        A ``ValueError`` is raised if no solution exists.

    EXAMPLES::

        sage: from vectors_in_intervals.utility import vector_from_sign_vector
        sage: from elementary_vectors import *
        sage: from sign_vectors import *
        sage: M = matrix([[1, 0, 2, 0], [0, 1, 1, 0], [0, 0, 0, 1]])
        sage: vector_from_sign_vector(M, zero_sign_vector(4))
        (0, 0, 0, 0)
        sage: vector_from_sign_vector(M, sign_vector("+-+0"))
        (2, -2, 2, 0)
        sage: vector_from_sign_vector(M, sign_vector("+0+0"))
        (1, 0, 2, 0)
        sage: vector_from_sign_vector(M, sign_vector("+-0+"))
        (1, -2, 0, 1)
        sage: evs = elementary_vectors(M, dual=False)
        sage: vector_from_sign_vector(evs, sign_vector("+-0+"))
        (1, -2, 0, 1)
        sage: vector_from_sign_vector(M, sign_vector("+0-0"))
        Traceback (most recent call last):
        ...
        ValueError: Cannot find vector corresponding to given sign vector.
        sage: vector_from_sign_vector([], zero_sign_vector(4))
        (0, 0, 0, 0)
    """
    if isinstance(data, list):
        evs = data
        try:
            result = data[0].parent().zero_vector()
        except IndexError:
            result = zero_vector(sv.length())
    elif isinstance(data, Generator):
        evs = data
        result = zero_vector(sv.length())
    else:
        evs_object = ElementaryVectors(data)
        # evs_object.set_combinations_dual(Combinations(upper.support(), evs_object.length - evs_object.rank + 1))
        evs = evs_object.generator(dual=False)
        result = zero_vector(data.base_ring(), sv.length())

    if sign_vector(result) == sv:
        return result
    for v in evs:
        for w in [v, -v]:
            if sign_vector(w) <= sv:
                result += w
                if sign_vector(result) == sv:
                    return result
                break

    raise ValueError("Cannot find vector corresponding to given sign vector.")


def solve_without_division(matrix: Matrix, rhs: vector):
    r"""
    Solve a linear system of equations without division.

    Compute ``x`` for ``A x = c b`` where ``c`` is any positive constant.
    Uses an elementary vector.

    EXAMPLES::

        sage: from vectors_in_intervals.utility import solve_without_division
        sage: A = matrix([[1, 2], [0, 1], [1, -1]])
        sage: b = vector([1, 0, 1])
        sage: solve_without_division(A, b)
        (1, 0)
        sage: A = matrix([[1, 4], [0, 2], [1, -2]])
        sage: b = vector([6, 2, 0])
        sage: solve_without_division(A, b)
        (4, 2)
        sage: A.solve_right(b)
        (2, 1)
        sage: A = matrix([[1, 1, 1], [0, 1, 2]])
        sage: b = vector([2, 3])
        sage: solve_without_division(A, b)
        (0, 1, 1)
    """
    ev = next(ElementaryVectors(
        Matrix.block([[matrix, Matrix.column(rhs)]])
    ).generator(reverse=True))
    return -sign(ev[-1]) * ev[:-1]


class CombinationsIncluding(SageObject):
    r"""
    Combinatorial object of all combinations that include given elements

    EXAMPLES:

    We generate all subsets of ``range(4)`` with ``2`` elements that include the element ``2``::

        sage: from vectors_in_intervals.utility import CombinationsIncluding
        sage: C = CombinationsIncluding(4, 2, [2])
        sage: list(C)
        [[0, 2], [1, 2], [2, 3]]
        sage: list(reversed(C))
        [[2, 3], [1, 2], [0, 2]]
    """
    def __init__(self, mset, k, elements=None):
        self.combinations = Combinations(set(range(mset)) - set(elements), k - len(elements))
        if elements is None:
            self.elements = []
        else:
            self.elements = elements

    def __iter__(self) -> Iterator[list]:
        for combination in self.combinations:
            yield sorted(list(self.elements) + combination)
    
    def __reversed__(self) -> Iterator[list]:
        for combination in reversed(self.combinations):
            yield sorted(list(self.elements) + combination)

    def random_element(self) -> list:
        return sorted(list(self.elements) + self.combinations.random_element())
