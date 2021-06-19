from sage.all import *
from elementary_vectors.exists_vector import exists_vector
import unittest

# Todo:
# * tests with output vector
# * tests with given elementary vectors
# * tests with kernel=False
# * cover special input cases (e.g. l = True, r = None)

class Tests(unittest.TestCase):
    def test_exists_vector_closed(self):
        M = matrix([1,1,0])
        L = [2,5,-1] # lower halves of the intervals
        R = [5,6,1] # upper halves of the intervals
        
        self.assertEqual(exists_vector(M, L, R), True)
        self.assertEqual(exists_vector(M, L, R, r=True), True)
        self.assertEqual(exists_vector(M, L, R, l=True, r=True), True)

    def test_exists_vector_open(self):
        M = matrix([1,1,0])
        L = [2,5,-1] # lower halves of the intervals
        R = [5,6,1] # upper halves of the intervals
        
        self.assertEqual(exists_vector(M, L, R, l=False, r=False), False)

    def test_exists_vector_mixed(self):
        M = matrix([1,1,0])
        L = [2,5,-1] # lower halves of the intervals
        R = [5,6,1] # upper halves of the intervals
        l = [True,True,False]
        r = [False,True,True]
        
        self.assertEqual(exists_vector(M, L, R, l=l, r=r), False)

    def test_exists_vector_infinity(self):
        M = matrix([[1,0,1,0],[0,1,1,1]])
        L = [2,5,0,-oo]
        R = [5,oo,8,5]
        l = [True,True,False,False]
        r = [False,False,False,True]
        
        self.assertEqual(exists_vector(M, L, R, l=l, r=r), True)
        

if __name__ == '__main__':
    unittest.main()
