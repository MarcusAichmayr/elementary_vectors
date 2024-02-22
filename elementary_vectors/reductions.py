r"""Reducing and simplifying lists of vectors."""

#############################################################################
#  Copyright (C) 2024                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)        #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

import operator
from sage.arith.misc import gcd
from sage.modules.free_module_element import vector
from sage.symbolic.ring import SR


def simplify_using_equalities(value, equalities):
    r"""
    Simplifies an expression using a list of equalities.
    
    .. NOTE::

        Only equalities are considered. Other expressions (e.g. inequalities) are ignored.

    EXAMPLES::

        sage: from elementary_vectors.reductions import simplify_using_equalities
        sage: var('a')
        a
        sage: expr = a + 1
        sage: simplify_using_equalities(expr, [a == 0])
        1
        sage: simplify_using_equalities(expr, [])
        a + 1

    Now, we consider a polynomial expression::

        sage: R = PolynomialRing(ZZ, 'x')
        sage: x = R.gen()
        sage: assume(SR(x) == 0)
        sage: simplify_using_equalities(x + 1, assumptions())
        1
    """
    expressions = []
    for equation in equalities:
        try:
            if equation.operator() is operator.eq:
                expressions.append(equation)
        except AttributeError:
            pass

    if equalities:
        # casting because substitute would work differently for polynomials
        expr = SR(value).substitute(expressions)
        try:
            return value.base_ring()(expr)
        except TypeError:
            return expr
    return value


def reduce_factor(iterable):
    r"""
    Cancel a common factor of each entry.
    
    INPUT:
    
    - ``iterable`` -- a vector or a list

    EXAMPLES::

        sage: from elementary_vectors.reductions import reduce_factor
        sage: var('a')
        a
        sage: v = vector([5*a, 10*a])
        sage: reduce_factor(v)
        (1, 2)
        sage: w = vector([4, 6])
        sage: reduce_factor(w)
        (2, 3)

    When we cancel a common factor, we expect to have an integer vector again::

        sage: type(reduce_factor(w))
        <class 'sage.modules.vector_integer_dense.Vector_integer_dense'>

    We can also cancel common factors from lists::

        sage: reduce_factor([4, 6])
        [2, 3]
        sage: reduce_factor([-4, -6])
        [-2, -3]

    The function also cancels denominators::

        sage: v = vector([1/10, 0, 1/3, 1/4])
        sage: v
        (1/10, 0, 1/3, 1/4)
        sage: reduce_factor(v)
        (6, 0, 20, 15)

    TESTS::

        sage: from elementary_vectors.reductions import reduce_factor
        sage: v = vector([2, 4])
        sage: w = reduce_factor(v)
        sage: w
        (1, 2)
        sage: w.base_ring()
        Integer Ring
        sage: v = [2, 4]
        sage: w = reduce_factor(v)
        sage: w
        [1, 2]
        sage: w[0].base_ring()
        Integer Ring
    """
    divisor = gcd(iterable)
    try:
        if isinstance(iterable, list):
            return type(iterable)(entry // divisor for entry in iterable)
        return (iterable / divisor).change_ring(iterable.base_ring())
    except ZeroDivisionError:
        return iterable


def reduce_vector_using_equalities(iterable, equalities):
    r"""Use a list of equalities ``equalities`` to simplify expressions in a vector ``iterable``."""
    if equalities:
        return vector(iterable.base_ring(), (simplify_using_equalities(entry, equalities=equalities) for entry in iterable))
    return iterable


def reduce_vector(element, equalities=None, cancel_factor=True):
    r"""
    Reduce this vector by canceling common factors and applying equalities to variables.

    INPUT:

    - ``element`` -- a vector

    - ``equalities`` -- a list of equalities (default: ``None``)

    - ``cancel_factor`` -- a boolean (default: ``True``). If true, cancels a common factor.

    EXAMPLES::

        sage: from elementary_vectors.reductions import reduce_vector
        sage: var('a')
        a
        sage: assume(a == 0)
        sage: v = vector([5*a, 10*a])
        sage: reduce_vector(v, equalities=assumptions())
        (0, 0)
        sage: reduce_vector(v)
        (1, 2)
        sage: reduce_vector(v, equalities=[a == 0], cancel_factor=False)
        (0, 0)
        sage: reduce_vector(v, cancel_factor=False)
        (5*a, 10*a)
    """
    if equalities:
        element = reduce_vector_using_equalities(element, equalities=equalities)
    if cancel_factor:
        element = reduce_factor(element)
    return element


def reduce_vectors_support(vectors, generator=False):
    r"""
    Return a ist of vectors where each element has distinct support.

    INPUT:

    - ``vectors`` -- an iterable of vectors
    
    - ``generator`` -- an optional boolean

    OUTPUT:

    Only the first element with a specific support is kept.
    Return a generator if ``generator`` is true.

    EXAMPLES::

        sage: from elementary_vectors.reductions import reduce_vectors_support
        sage: l = [vector([1, 3, 2]), vector([0, 0, 1]), vector([2, 2, 0]), vector([0, 0, -5])]
        sage: reduce_vectors_support(l)
        [(1, 3, 2), (0, 0, 1), (2, 2, 0)]
        sage: reduce_vectors_support([zero_vector(5)])
        [(0, 0, 0, 0, 0)]
    """
    def reduce_vectors_support_generator(vectors):
        r"""Return a generator of vectors where elements with same support are removed."""
        checked_supports = set()
        for element in vectors:
            support = frozenset(element.support())
            if support not in checked_supports:
                checked_supports.add(support)
                yield element

    if generator:
        return reduce_vectors_support_generator(vectors)
    return list(reduce_vectors_support_generator(vectors))


def remove_zero_vectors(vectors):
    r"""
    Remove all zero vectors from this list.

    EXAMPLES::

        sage: from elementary_vectors.reductions import remove_zero_vectors
        sage: var('a')
        a
        sage: remove_zero_vectors([vector([5, 2, 3]), zero_vector(3), vector([a*0, 0, 0])])
        [(5, 2, 3)]
    """
    return [v for v in vectors if v]


def reduce_vectors(vectors, equalities=None, cancel_factors=False, reduce_support=True, remove_zeros=True):
    r"""
    Reduces this list of vectors.

    INPUT:

    - ``vectors`` -- a list of vectors

    - ``equalities`` -- a list of equalities (default: ``None``)

    - ``cancel_factors`` -- a boolean (default: ``False``). If true, cancels common factors of each vector.

    - ``reduce_support`` -- a boolean (default: ``True``). Keeps only the first vector for each different support.

    - ``remove_zeros`` -- a boolean (default: ``True``). If true, removes all zero vectors.

    EXAMPLES::

        sage: from elementary_vectors.reductions import reduce_vectors
        sage: var('a')
        a
        sage: assume(a == 0)
        sage: l = [vector([5*a, 10*a, 0]), vector([5*a, 2*a, a]), vector([4, 6, 0])]
        sage: reduce_vectors(l, equalities=assumptions(), cancel_factors=True)
        [(2, 3, 0)]
        sage: reduce_vectors(l, cancel_factors=True)
        [(1, 2, 0), (5, 2, 1)]
        sage: reduce_vectors(l, equalities=[a == 0], cancel_factors=False)
        [(4, 6, 0)]
        sage: reduce_vectors(l, equalities=assumptions(), reduce_support=True, cancel_factors=True)
        [(2, 3, 0)]
        sage: reduce_vectors(l, equalities=assumptions(), reduce_support=False, remove_zeros=False, cancel_factors=True)
        [(0, 0, 0), (0, 0, 0), (2, 3, 0)]
    """
    vectors = [reduce_vector(v, equalities=equalities, cancel_factor=cancel_factors) for v in vectors]
    if remove_zeros:
        vectors = remove_zero_vectors(vectors)
    if reduce_support:
        vectors = reduce_vectors_support(vectors)
    return vectors
