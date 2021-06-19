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
from sage.modules.free_module_element import zero_vector

def elementary_vectors(M, kernel=False):
    r"""
    Compute elementary vectors of a subspace determined by the matrix ``M``.
    
    INPUT:
    
    - ``M`` -- a matrix
    
    - ``kernel`` -- a boolean (default: ``False``)
    
    OUTPUT:
    
    - If ``kernel`` is true, returns a list of elementary vectors lying in
      the kernel of ``M``.
    
    - If ``kernel`` is false, returns a list of elementary vectors lying in
      the row space of ``M``
    """
    n = M.ncols()
    if kernel:
        ind = M.pivot_rows()
        M1 = M.matrix_from_rows(ind)
    else:
        M1 = M.right_kernel_matrix()
    r = M1.nrows()
    L = []
    if r >= n:
        return L
    C = Combinations(n,r+1)
    dets = Combinations(n,r)
    minM = M1.minors(r)
    Base = M.base_ring()
    S = []
    for I in C:
        if not any([all([k in I for k in vs]) for vs in S]):
            v = zero_vector(Base, n)
            for k in I:
                pos = I.index(k)
                Jk = list(I)
                Jk.pop(pos)
                v[k] = (-1)**pos * minM[dets.rank(Jk)]
            s = v.support()
            if s != [] and len(s) <= r:
                S.append(s)
            if v != zero_vector(n):
                L.append(v)
    return L
