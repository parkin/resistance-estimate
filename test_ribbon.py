import unittest
from ribbon import parallel
from ribbon import series


class TestSeries(unittest.TestCase):
        def test_normal_2(self):
                r1 = 10
                r2 = 10

                r_out = series(r1, r2)
                self.assertEqual(r_out, 20)

                r1 = 5
                r2 = 10

                r_out = series(r1, r2)
                self.assertEqual(r_out, 15)

        def test_normal_3(self):
                r1 = 10
                r2 = 10
                r3 = 10

                r_out = series(r1, r2, r3)
                self.assertEqual(r_out, 30)

                r1 = 5
                r2 = 10
                r3 = 15

                r_out = series(r1, r2, r3)
                self.assertEqual(r_out, 30)




class TestParallel(unittest.TestCase):
        def test_normal_2(self):
                r1 = 10
                r2 = 10
                
                r_out = parallel(r1, r2)

                self.assertEqual(r_out, 5)

                r1 = 5
                r2 = 20

                r_out = parallel(r1, r2)
                self.assertEqual(r_out, 4)

        def test_normal_3(self):
                r1 = 30
                r2 = 30
                r3 = 30

                r_out = parallel(r1, r2, r3)
                self.assertEqual(r_out, 10)

                r1 = 2
                r2 = 5
                r3 = 10

                r_out = parallel(r1, r2, r3)
                self.assertEqual(r_out, 1.25)
