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

from .intervals import Interval, Intervals
from .existence import exists_vector, exists_orthogonal_vector
from .linear_inequality_systems import (
    LinearInequalitySystem,
    InhomogeneousSystem,
    HomogeneousSystem,
    # HomogeneousSystemCocircuits,
)
from .certifying_inequalities import (
    AlternativesGeneral,
    AlternativesInhomogeneous,
    AlternativesHomogeneous,
)
