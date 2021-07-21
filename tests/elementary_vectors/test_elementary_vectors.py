from sage.all import *
from elementary_vectors.elementary_vectors import elementary_vectors, elementary_vectors_from_matrix, elementary_vectors_from_minors
from elementary_vectors.elementary_vectors import non_negative_vectors, positive_elementary_vectors
import unittest

class Tests(unittest.TestCase):
    def test_elementary_vectors_from_matrix(self):
        M = matrix([[0,0,1,-1,0],[2,0,0,0,2],[1,1,1,1,1]])
        evs1 = elementary_vectors_from_matrix(M, kernel=True)
        evs2 = elementary_vectors_from_matrix(M, kernel=True, reduce=False)
        evs3 = elementary_vectors_from_matrix(M, kernel=False)
        
        self.assertEqual(len(evs1), 2)
        self.assertEqual(len(evs2), 5)
        self.assertEqual(len(evs3), 4)
        self.assertTrue(vector([0,2,-1,-1,0]) in evs1)

        self.assertEqual(type(evs1[0]), type(evs2[0])) # both should be integer vectors

        A = identity_matrix(3)
        self.assertEqual(elementary_vectors_from_matrix(A),[])

        mM = M.minors(3)

        self.assertEqual(elementary_vectors_from_minors(mM,[3,5]), evs1)
        self.assertEqual(elementary_vectors_from_minors(mM,[3,5], reduce=False), evs2)

    def test_elementary_vectors_polynomial_matrix(self):
        R = PolynomialRing(ZZ, 'x')
        x = R.gen()
        A = matrix([[x, 1, 0],[0,0,2]])
        
        self.assertEqual(elementary_vectors_from_matrix(A), [vector([1,-x,0])])
        
        # elementary_vectors_from_matrix(A, kernel=False) # not implemented for matrices of ZZ[x]
        
        R = PolynomialRing(QQ, 'x')
        x = R.gen()
        B = matrix([[x, 1, 0],[0,0,2]])
        
        elementary_vectors_from_matrix(B, kernel=False) # works for QQ[x]
    
    def test_elementary_vectors_from_matrix(self):
        M = matrix([[0,0,1,-1,0],[2,0,0,0,2],[1,1,1,1,1]])
        evs1 = elementary_vectors_from_matrix(M, kernel=True)
        evs2 = elementary_vectors_from_matrix(M, kernel=True, reduce=False)
        
        mM = M.minors(3)

        self.assertEqual(elementary_vectors(mM,[3,5]), evs1)
        self.assertEqual(elementary_vectors(mM,[3,5], reduce=False), evs2)
        
        args = {
            "kernel":True,
            "reduce":False,
            "ring":RR,
            "return_minors":True
        }
        self.assertEqual(elementary_vectors(M, **args), elementary_vectors(mM, M.dimensions(), **args))

    def test_positive_elementary_vectors(self):
        A = random_matrix(ZZ,3,5)
        positive_elementary_vectors(A)
        
        # TODO

# TODO test remaining functions

if __name__ == '__main__':
    unittest.main()
