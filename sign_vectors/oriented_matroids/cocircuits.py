#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sign_vectors import sign_vector
from sign_vectors import zero_sign_vector

from elementary_vectors import elementary_vectors

def cocircuits_from_matrix(A, kernel=False):
    r"""
    Computes a list of cocircuits determined by the matrix ``A``.
    
    INPUT:
    
    - ``A`` -- a matrix with real arguments.

    - ``kernel`` -- a boolean (default: False)
    
    OUTPUT:
    
    - If ``kernel`` is false, returns a list of cocircuits determined by the row
      space of the matrix ``A`` (default).

    - If ``kernel`` is true, returns a list of cocircuits determined by the
      kernel of the matrix ``A``.
    """
    L = elementary_vectors(A, kernel=kernel)
    return [sign_vector(v) for v in L] + [sign_vector(-v) for v in L]


def covectors_from_cocircuits(L):
    r"""
    Uses a list of cocircuits to compute all covectors of the corresponding
    oriented matroid.
    
    INPUT:
    
    - ``L`` -- a list of cocircuits of an oriented matroid.

    OUTPUT:
    
    - a list of all covectors of the oriented matroid.
    """
    assert L, 'List of cocircuits is empty.'
    n = L[0].length()
    F = [zero_sign_vector(n)]
    F_new = [zero_sign_vector(n)]
    while F_new != []:
        Y = F_new.pop()
        for X in L:
            if not X <= Y: # otherwise Z = X.compose(Y) = Y in F
                Z = X.compose(Y)
                if Z not in F:
                    F.append(Z)
                    F_new.append(Z)
    return F
