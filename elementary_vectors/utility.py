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
from sign_vectors import sign_vector

# Todo: change name?
def reduce_by_support(L):
    r"""
    Returns a sublist of vectors where each vector has distinct support.
    
    INPUT:
    
    - ``L`` -- a list of vectors
    
    OUTPUT:
    
    Returns a sublist of ``L`` such that each vector has distinct support.
    """
    supp = []
    out = []
    for v in L:
        s = v.support()
        if s not in supp:
            supp.append(s)
            out.append(v)
    return out


def non_negative_vectors(L):
    r"""
    Returns non-negative vectors.
    
    INPUT:
    
    - ``L`` -- a list of vectors
    
    OUTPUT:
    
    Returns all vectors of ``L`` that are
    - non_negative in each component; or
    - negative in each component. Those will be multiplied by ``-1``; or
    - containing variables such that no opposing signs occur.
    """
    out = []
    for v in L:
        if sign_vector(v) >= 0: # Use ``>=`` instead of ``>``, (0,0,x) -> (000) should be returned
            out.append(v)
        elif sign_vector(v) < 0:
            out.append(-v)
    return out


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
