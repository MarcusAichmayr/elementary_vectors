#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from .topes import *
from .cocircuits import *
from .face_enumeration import face_enumeration

from sage.misc.flatten import flatten

# Todo: could be useful to put these functions in other files

def topes_from_matrix(A, kernel=False):
    r"""
    Return a list of topes of the oriented matroid corresponding to the matrix ``A``.
    
    INPUT:
    
    - ``A`` -- a matrix
    
    - ``kernel`` -- a boolean (default: ``False``)
    
    OUTPUT:
    
    - If ``kernel`` is false, returns a list of topes determined by the row space
      of the matrix ``A``.

    - If ``kernel`` is true, returns a list of topes determined by the kernel of
      the matrix ``A``.
    """
    return topes_from_cocircuits(cocircuits_from_matrix(A, kernel=kernel))


def covectors_from_topes(T, separate=False):
    r"""
    Computes all covectors from a list of topes.
    
    INPUT:
    
    - ``T`` -- a list of topes.
    
    - ``separate`` -- a boolean (default: ``False``)
    
    OUTPUT:
    The list of covectors of the corresponding oriented matroid.
    
    - If ``separate`` is false, returns a list of covectors. The covectors are
      sorted by rank (default).

    - If ``separate`` is true, returns a list of lists of covectors, separated
      by their rank.
    """
    if separate:
        return face_enumeration(T)
    else:
        return flatten(face_enumeration(T))


def cocircuits_from_topes(T):
    r"""
    Computes all cocircuits from a list of topes.
    
    INPUT:
    
    - ``T`` -- a list of topes.
    
    OUTPUT:
    A list of cocircuits of the corresponding oriented matroid.
    """
    return face_enumeration(T)[1]


def covectors_from_matrix(A, kernel=False, algorithm=None, separate=False):
    r"""
    Return a list of covectors of the oriented matroid corresponding to the
    matrix ``A``.
    
    INPUT:
    
    - ``A`` -- a matrix.
    
    - ``kernel`` -- a boolean (default: ``False``)

    - ``algorithm`` -- (optional) either 'face_enumeration' or 'fe'
    
    - ``separate`` -- a boolean (default: ``False``)
    
    OUTPUT:
    
    The list of covectors of an oriented matroid corresponding to the matrix ``A``.
    
    - If ``kernel`` is false, the returned covectors will be determined by the
      row space of the matrix ``A``.

    - If ``kernel`` is true, the returned covectors will be determined by the
      kernel of the matrix ``A``.

    - If ``algorithm`` is 'face_enumeration' or the shortcut 'fe',
      
      - if ``separate`` is false, returns a list of covectors. The covectors are
        sorted by rank (default).
      
      - if ``separate`` is true, returns a list of lists of covectors, separated
        by their rank.
    """
    if algorithm is None:
        return covectors_from_cocircuits(cocircuits_from_matrix(A, kernel=kernel))
    elif algorithm == 'face_enumeration' or algorithm == 'fe':
        return covectors_from_topes(topes_from_matrix(A, kernel=kernel), separate=separate)
    else:
        raise ValueError("no algorithm '%s'"%algorithm)
