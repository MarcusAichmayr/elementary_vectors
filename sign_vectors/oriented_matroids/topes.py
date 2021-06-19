#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sign_vectors import zero_sign_vector
from sign_vectors.utility import loops

def topes_from_cocircuits(D):
    r"""
    Uses the cocircuits of an oriented matroid to compute the topes.
    
    INPUT:
    
    - ``D`` -- a list of cocircuits of an oriented matroid.
    
    OUTPUT:
    
    - a list of topes of the oriented matroid.
    """
    assert D, 'List is empty.'
    n = D[0].length()
    
    F = [zero_sign_vector(n)]
    F_new = [zero_sign_vector(n)]
    T = []
    E0 = loops(D) # intersection of zero-supports of all X in D
    
    while F_new != []:
        Y = F_new.pop()
        for X in D:
            if not X <= Y: # otherwise Z = X.compose(Y) = Y in F
                Z = X.compose(Y)
                if Z not in F:
                    F.append(Z)
                    if Z.zero_support() == E0:
                        T.append(Z)
                    else:
                        F_new.append(Z)
    return T
