import unittest
import random
import string
import mummer

def gen_seq(size, chars="ACGT"):
    return ''.join(random.choice(chars) for _ in range(size))

class TestNucmer(unittest.TestCase):
    def test_exact_pair(self):
        s1 = gen_seq(100)
        s2 = s1[80:] + gen_seq(80)
        o = mummer.Options()
        o.minmatch(10)
        o.mincluster(10)
        a = mummer.align_sequences(s1, s2, o)
        self.assertEqual(1, len(a))
        self.assertEqual(1, a[0].dirB)
        self.assertEqual( 81, a[0].sA)
        self.assertEqual( 100, a[0].eA)
        self.assertEqual( 1, a[0].sB)
        self.assertEqual( 20, a[0].eB)
        self.assertEqual( 0, a[0].Errors)
        self.assertEqual( 0, a[0].SimErrors)
        self.assertEqual( 0, len(a[0].delta))
        self.assertEqual( 1.0, a[0].identity())
        self.assertEqual( 1.0, a[0].similarity())


if __name__ == '__main__':
    unittest.main()

