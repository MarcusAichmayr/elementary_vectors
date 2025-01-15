r"""Finding vectors with components in intervals"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from __future__ import absolute_import

from .setup_intervals import (
    intervals_from_bounds,
    is_vector_in_intervals,
    random_intervals,
    intervals_from_sign_vector,
)
from .existence import exists_vector, exists_orthogonal_vector
from .construction import (
    construct_vector,
    construct_orthogonal_vector,
    vector_from_sign_vector,
    vector_between_sign_vectors,
    sign_vectors_in_intervals,
)
from .linear_inequality_systems import (
    LinearInequalitySystem,
    InhomogeneousSystem,
    HomogeneousSystem,
    HomogeneousSystemCocircuits,
    homogeneous_from_general,
    inhomogeneous_from_general,
    homogeneous_from_inhomogeneous,
)
from .certifying_inequalities import (
    AlternativesGeneral,
    AlternativesInhomogeneous,
    AlternativesHomogeneous,
)
