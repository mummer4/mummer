import unittest
import random
import string
import mummer

def gen_seq(size, chars="ACGT"):
    return ''.join(random.choice(chars) for _ in range(size))

def find(array, compare_function):
    for i, v in enumerate(array):
        if compare_function(v):
            return v
    return None

class TestNucmer(unittest.TestCase):
    def test_exact_pair(self):
        s1 = gen_seq(100)
        s2 = s1[80:] + gen_seq(80)
        o = mummer.Options()
        o.minmatch(10)
        o.mincluster(10)
        a = mummer.align_sequences(s1, s2, o)

        self.assertTrue(1 <= len(a))
        al = find(a, lambda x: x.sA == 81 and x.eA == 100)
        self.assertIsNotNone(al)
        self.assertEqual(1, al.dirB)
        self.assertEqual(81, al.sA)
        self.assertEqual(100, al.eA)
        self.assertEqual(1, al.sB)
        self.assertEqual(20, al.eB)
        self.assertEqual(0, al.Errors)
        self.assertEqual(0, al.SimErrors)
        self.assertEqual(0, len(al.delta))
        self.assertEqual(1.0, al.identity())
        self.assertEqual(1.0, al.similarity())

    def test_thread(self):
        nt = mummer.get_num_threads()
        self.assertTrue(nt >= 1)
        mummer.set_num_threads(nt + 1)
        nnt = mummer.get_num_threads()
        self.assertTrue(nnt == 1 or nnt == nt + 1)


if __name__ == '__main__':
    unittest.main()

