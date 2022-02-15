r"""
Reducing lists of vectors.

This module offers several functions to reduce a list of vector.
"""

#############################################################################
#  Copyright (C) 2022                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

import operator
from sage.modules.free_module_element import vector
from sage.arith.misc import gcd
from sage.symbolic.ring import SR


def simplify_using_equalities(a, eq):
    r"""
    Simplifies the expression ``a`` using a list of equalities ``eq``. Only considers equalities in ``eq``. Other expressions (e.g. inequalities) are ignored.

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
        sage: simplify_using_equalities(x+1, assumptions())
        1
    """
    # There might be inequalities in ``eq`` or other expressions like ``a is real``.
    # We get rid of them:
    l = []
    for s in eq:
        try:
            if s.operator() is operator.eq:
                l.append(s)
        except AttributeError:
            pass

#    eq = [s for s in eq if s.operator() is operator.eq]
    if eq:
        # We cast ``a`` to SR. Then, we can substitute and cast the result back.
        # Substituting directly into ``a``, does not work for polynomials
        # since substitute works differently.
        expr = SR(a).substitute(l)
        try:
            return a.base_ring()(expr)
        except TypeError:
            return expr
    else:
        return a


# TODO: Should ``reduce_factor`` reduce by e.g. ``a`` if we have ``assume(a == 0)``?
# For ``a == 0``, ``[4*a, 6*a]`` could be reduced to either ``[2, 3]`` or ``[2 a, 3 a]``.
# What do we expect?
def reduce_factor(v):
    r"""
    Cancel a common factor of each entry of this vector. Also works for lists.

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
        sage: reduce_factor([-4, -6]) # TODO: Do we expect this?
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
    g = gcd(v)
    try:
        if isinstance(v, list):
            return type(v)(vi // g for vi in v)
        else:
            return (v/g).change_ring(v.base_ring())
    except ZeroDivisionError:
        return v


def reduce_vector_using_equalities(v, eq):
    r"""Use a list of equalities ``eq`` to simplify expressions in a vector ``v``."""
    if eq:
        return vector(v.base_ring(), [simplify_using_equalities(vi, eq=eq) for vi in v])
    else:
        return v


def reduce_vector(v, eq=None, cancel_factor=True):
    r"""
    Reduces this vector.

    INPUT:

    - ``v`` -- a vector

    - ``eq`` -- a list of equalities (default: ``None``)

    - ``cancel_factor`` -- a boolean (default: ``True``). If true, cancels a common factor.

    EXAMPLES::

        sage: from elementary_vectors.reductions import reduce_vector
        sage: var('a')
        a
        sage: assume(a == 0)
        sage: v = vector([5*a, 10*a])
        sage: reduce_vector(v, eq=assumptions())
        (0, 0)
        sage: reduce_vector(v)
        (1, 2)
        sage: reduce_vector(v, eq=[a == 0], cancel_factor=False)
        (0, 0)
        sage: reduce_vector(v, cancel_factor=False)
        (5*a, 10*a)
    """
    if eq:
        v = reduce_vector_using_equalities(v, eq=eq)
    if cancel_factor:
        v = reduce_factor(v)
    return v


def reduce_vectors_support(L):
    r"""
    Return a sublist of vectors where each vector has distinct support.

    INPUT:

    - ``L`` -- a list of vectors

    OUTPUT:

    Keeps only one vector for each support.

    EXAMPLES::

        sage: from elementary_vectors.reductions import reduce_vectors_support
        sage: l = [vector([1,3,2]), vector([0,0,1]), vector([2,2,0]), vector([0,0,-5])]
        sage: reduce_vectors_support(l)
        [(1, 3, 2), (0, 0, 1), (2, 2, 0)]
        sage: reduce_vectors_support([zero_vector(5)])
        [(0, 0, 0, 0, 0)]
    """
    supp = []
    out = []
    for v in L:
        s = v.support()
        if s not in supp:
            supp.append(s)
            out.append(v)
    return out


def remove_zero_vectors(L):
    r"""
    Remove all zero vectors from this list.

    EXAMPLES::

        sage: from elementary_vectors.reductions import remove_zero_vectors
        sage: var('a')
        a
        sage: remove_zero_vectors([vector([5,2,3]), zero_vector(3), vector([a*0,0,0])])
        [(5, 2, 3)]
    """
    return [v for v in L if v != 0]


def reduce_vectors(L, eq=None, cancel_factors=True, reduce_support=True, remove_zeros=True):
    r"""
    Reduces this list of vectors.

    INPUT:

    - ``L`` -- a list of vectors

    - ``eq`` -- a list of equalities (default: ``None``)

    - ``cancel_factors`` -- a boolean (default: ``True``). If true, cancels common factors of each vector.

    - ``reduce_support`` -- a boolean (default: ``True``). Keeps only the first vector for each different support.

    - ``remove_zeros`` -- a boolean (default: ``False``). If true, removes all zero vectors.

    EXAMPLES::

        sage: from elementary_vectors.reductions import reduce_vectors
        sage: var('a')
        a
        sage: assume(a == 0)
        sage: l = [vector([5*a, 10*a, 0]), vector([5*a, 2*a, a]), vector([4, 6, 0])]
        sage: reduce_vectors(l, eq=assumptions())
        [(2, 3, 0)]
        sage: reduce_vectors(l)
        [(1, 2, 0), (5, 2, 1)]
        sage: reduce_vectors(l, eq=[a == 0], cancel_factors=False)
        [(4, 6, 0)]
        sage: reduce_vectors(l, eq=assumptions(), reduce_support=True)
        [(2, 3, 0)]
        sage: reduce_vectors(l, eq=assumptions(), reduce_support=False, remove_zeros=False)
        [(0, 0, 0), (0, 0, 0), (2, 3, 0)]
    """
    L = [reduce_vector(v, eq=eq, cancel_factor=cancel_factors) for v in L]
    if remove_zeros:
        L = remove_zero_vectors(L)
    if reduce_support:
        L = reduce_vectors_support(L)
    return L
