#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sign_vectors import SignVector
from sign_vectors import sign_vector

from sage.modules.free_module_element import vector

def subvector(F, R):
    r"""Returns a function that returns a sign vector or vector consisting of entries not in ``R``."""
    assert F, 'List is empty.'
    n = len(list(F[0]))
    S = [e for e in range(n) if e not in R] # S = E\R

    if isinstance(F[0], SignVector):
        def vec(v):
            return sign_vector(v.list_from_positions(S))
    else:
        def vec(v):
            return vector(v.list_from_positions(S))
    return vec


def contraction(F, R, keep_components=False):
    r"""
    Returns all sign vectors that are zero on ``R``. Also works for real vectors.
    
    INPUT:
    
    - ``F`` -- a list of sign vectors, a list of vectors or a matrix.
    
    - ``R`` -- a list of indices.
    
    - ``keep_components`` -- a boolean (default: ``False``).
        
        - If ``keep_components`` is false (default), remove entries in ``R``.
        
        - If ``keep_components`` is true, keep entries in ``R``.
    """
    if F == []:
        return F


    if keep_components:
        def vec(v):
            return v
    else:
        vec = subvector(F, R)
    
    L = []
    
    for X in F:
        val = True
        for e in X.support():
            if e in R:
                val = False
                break
        if val:
            L.append(vec(X))
    return L


def deletion(F, R):
    r"""
    Removes the components corresponding to ``R`` from a list of sign vectors.
    Also works for real vectors.
    
    INPUT:
    
    - ``F`` -- a list of sign vectors, a list of real vectors or a matrix.
    
    - ``R`` -- a list of indices.
    """
    if F == []:
        return F
    
    vec = subvector(F, R)

    L = []
    for X in F:
        X_R = vec(X)
        if X_R not in L: # naive, might be inefficient 
            L.append(X_R)
    return L
