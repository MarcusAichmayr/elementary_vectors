#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from __future__ import absolute_import

from .sign_vectors import *

from .closure import closure

from .contraction_deletion import contraction, deletion

from .oriented_matroids import cocircuits_from_matrix, covectors_from_cocircuits, topes_from_cocircuits, lower_faces, face_enumeration, topes_from_matrix, covectors_from_topes, cocircuits_from_topes, covectors_from_matrix

# from .utility import loops, is_parallel, parallel_classes, classes_same_support

