"""Utility functions"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################


def is_symbolic(expression):
    r"""
    Return whether this element is a symbolic expression.

    If it belongs to the symbolic ring but doesn't contain any variables it does not count as "symbolic".

    EXAMPLES::

        sage: from elementary_vectors.utility import is_symbolic
        sage: is_symbolic(5)
        False
        sage: var('a, b')
        (a, b)
        sage: is_symbolic(a)
        True
        sage: is_symbolic(-a)
        True
        sage: is_symbolic(b^2 - a)
        True
        sage: is_symbolic(SR(5))
        False
    """
    if hasattr(expression, "variables"):
        return bool(expression.variables())
    return False
