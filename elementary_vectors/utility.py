#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.modules.free_module_element import vector
from sage.symbolic.ring import SR

# TODO: improve name
def sign_determined(a):
    r"""Checks whether the sign of ``a`` is determined."""
    if SR(a) > 0 or SR(a) < 0 or SR(a) == 0:
        return True
    else:
        return False


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
    if x.length() != y.length():
        raise ValueError('Vectors have different length.')
    X = sign_vector(x)
    Y = sign_vector(y)
    D = X.separating_elements(Y)
    if D == []:
        raise ValueError('List of separating elements is empty.')
    if S == []:
        S = D
    elif not all([s in D for s in S]):
        raise ValueError('S is not a subset of D.')
    lam = max([x[e]/y[e] for e in S]) # x[e]/y[e] < 0 since e in D
    return x - lam*y
