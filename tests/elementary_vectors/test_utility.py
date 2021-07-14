from sage.all import *
from elementary_vectors.utility import reduce_by_support, non_negative_vectors, conformal_elimination
import unittest

class Tests(unittest.TestCase):
    def test_reduce_by_support(self):
        l = [vector([1,3,2]), vector([0,0,1]), vector([2,2,0]), vector([0,0,-5])]
        
        self.assertEqual(reduce_by_support(l), l[0:3])

# Todo test non_negative_vectors

# Todo test conformal_elimination

if __name__ == '__main__':
    unittest.main()
