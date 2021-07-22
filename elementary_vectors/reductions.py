#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.modules.free_module_element import vector, zero_vector
from sage.arith.misc import gcd
from sage.symbolic.ring import SR
from sage.symbolic.assumptions import assumptions
import operator

def simplify_using_equalities(a, eq):
    r"""
    Returns ``0`` if ``a`` is considered zero symbolically (e.g. by assumptions on variabls occuring in ``a``).
    Otherwise ``a`` is returned.
    TODO
    """
    # There might be inequalities in ``eq`` or other expressions like ``a is real``. We get rid of them:
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
        # Substituting directly into ``a``, does not work for polynomials since substitute works differently.
        # For instance, if there are no assumptions, then a segmentation fault would be triggered.
        expr = (SR(a).substitute(l))
        try:
            return a.base_ring()(expr)
        except TypeError:
            return expr
    else:
        return a


def reduce_factor(v):
    r"""Cancels a common factor of each entry of this vector."""
    g = gcd(v)
    try:
#        return vector(v.base_ring(), [vi/g for vi in v])
        return (v/g).change_ring(v.base_ring())
    except ZeroDivisionError:
        return v
        

def reduce_vector_using_equalities(v, eq):
    r"""Use a list of equalities ``eq`` to simplify expressions in a vector ``v``."""
    if eq:
        return vector(v.base_ring(), [simplify_using_equalities(vi, eq=eq) for vi in v])
    else:
        return v


def reduce_vector(v, eq=None, factor=True):
    r"""
    Reduces this vector.
    
    INPUT:
    
    - ``v`` -- a vector
    
    - ``eq`` -- a list of equalities (default: ``None``)

    - ``factor`` -- a boolean (default: ``True``). If true, cancels a common factor.
    """
    if eq:
        v = reduce_vector_using_equalities(v, eq=eq)
    if factor:
        v = reduce_factor(v)
    return v


def reduce_vectors_support(L):
    r"""
    Returns a sublist of vectors where each vector has distinct support.
    
    INPUT:
    
    - ``L`` -- a list of vectors
    
    OUTPUT:
    
    Keeps only one vector for each support.
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
    r"""Removes all zero vectors from this list."""
    return [v for v in L if v != 0]


def reduce_vectors(L, eq=None, factors=True, support=True, remove_zeros=True):
    r"""
    Reduces this list of vectors.
    
    INPUT:
    
    - ``L`` -- a list of vectors
    
    - ``eq`` -- a list of equalities (default: ``None``)

    - ``factors`` -- a boolean (default: ``True``). If true, cancels common factors of each vector.

    - ``support`` -- a boolean (default: ``True``). Keeps only the first vector for each different support.
    
    - ``remove_zeros`` -- a boolean (default: ``False``). If true, removes all zero vectors.
    """
    L = [reduce_vector(v, eq=eq, factor=factors) for v in L]
    if remove_zeros:
        L = remove_zero_vectors(L)
    if support:
        L = reduce_vectors_support(L)
    return L
