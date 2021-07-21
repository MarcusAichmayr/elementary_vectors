#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.combinat.combination import Combinations
from sage.modules.free_module_element import vector, zero_vector
from sage.structure.element import get_coercion_model
from sage.functions.other import binomial
from sage.arith.misc import gcd
from sage.symbolic.ring import SR
from sage.symbolic.assumptions import assume, forget

def reduce_by_support(L):
    r"""
    Returns a sublist of vectors where each vector has distinct support.
    
    INPUT:
    
    - ``L`` -- a list of vectors
    
    OUTPUT:
    
    Returns a sublist of ``L`` such that each vector has distinct support.
    Also removes zero vectors.
    """
    supp = [[]]
    out = []
    for v in L:
        s = v.support()
        if s not in supp:
            supp.append(s)
            out.append(v)
    return out


# TODO: improve name
def has_sign(a):
    r"""Checks whether the sign of ``a`` is determined."""
    if SR(a) > 0 or SR(a) < 0 or SR(a) == 0:
        return True
    else:
        return False

def reduce_if_zero(a):
    r"""
    Returns ``0`` if ``a`` is considered zero symbolically (e.g. by assumptions on variabls occuring in ``a``).
    Otherwise ``a`` is returned.
    """
    if SR(a).is_zero():
        return 0
    else:
        return a
    
def reduce_zero_entries_of_vector(v):
    r"""Replaces symbolic entries of a vector by ``0`` if the corresponding expression is zero."""
    return vector(v.base_ring(), [reduce_if_zero(vi) for vi in v])

def reduce_zero_entries(L):
    r"""Replaces symbolic entries of each vector in the list ``L`` by ``0`` if the corresponding expression is zero."""
    return [reduce_zero_entries_of_vector(v) for v in L]


def conformal_elimination(x, y, S=[]):
    r"""
    Applies conformal elimination to two real vectors to find a new vector.
    
    INPUT:
    
    - ``x`` -- a real vector
    
    - ``y`` -- a real vector
    
    - ``S`` -- a list of indices (default: ``[]``)
    
    OUTPUT:
    
    Returns a new vector ``z = x + a y`` where ``a > 0``, such that ``z[e] == 0``
    for some ``e`` in ``S`` and ``Z_S <= X_S`` and ``Z_f = (X o Y)_f`` for ``f``
    not in ``D(X, Y)``. Here, ``X``, ``Y`` and ``Z`` are the sign vectors
    corresponding to ``x``, ``y`` and ``z``.
        
    .. NOTE::
        
        If ``S`` is the empty list ``[]``, the whole list of separating elements
        will be considered instead. (default)
    """
    assert x.length() == y.length(), 'Vectors have different length.'
    X = sign_vector(x)
    Y = sign_vector(y)
    D = X.separating_elements(Y)
    assert D, 'List of separating elements is empty.'
    if S == []:
        S = D
    else:
        assert all([s in D for s in S]), 'S is not a subset of D.'
    lam = max([x[e]/y[e] for e in S]) # x[e]/y[e] < 0 since e in D
    return x - lam*y
