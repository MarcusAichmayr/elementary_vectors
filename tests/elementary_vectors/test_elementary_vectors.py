from sage.all import *
from elementary_vectors.elementary_vectors import elementary_vectors, elementary_vectors_from_minors
import unittest

class Tests(unittest.TestCase):
    def test_elementary_vectors(self):
        M = matrix([[0,0,1,-1,0],[2,0,0,0,2],[1,1,1,1,1]])
        evs1 = elementary_vectors(M, kernel=True)
        evs2 = elementary_vectors(M, kernel=True, reduce=False)
        evs3 = elementary_vectors(M, kernel=False)
        
        self.assertEqual(len(evs1), 2)
        self.assertEqual(len(evs2), 5)
        self.assertEqual(len(evs3), 4)
        self.assertTrue(vector([0,2,-1,-1,0]) in evs1)

        self.assertEqual(type(evs1[0]), type(evs2[0])) # both should be integer vectors

        A = identity_matrix(3)
        self.assertEqual(elementary_vectors(A),[])

        mM = M.minors(3)

        self.assertEqual(elementary_vectors_from_minors(mM,[3,5]), evs1)
        self.assertEqual(elementary_vectors_from_minors(mM,[3,5], reduce=False), evs2)

    def test_elementary_vectors_polynomial_matrix(self):
        R = PolynomialRing(ZZ, 'x')
        x = R.gen()
        A = matrix([[x, 1, 0],[0,0,2]])
        
        self.assertEqual(elementary_vectors(A), [vector([1,-x,0])])
        
        # elementary_vectors(A, kernel=False) # not implemented for matrices of ZZ[x]
        
        R = PolynomialRing(QQ, 'x')
        x = R.gen()
        B = matrix([[x, 1, 0],[0,0,2]])
        
        elementary_vectors(B, kernel=False) # works for QQ[x]

if __name__ == '__main__':
    unittest.main()
