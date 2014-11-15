import unittest
from estimate import parallel
from estimate import series
from estimate import calculate_dv
from estimate import calculate_r_out
from estimate import calculate_r_pre


class TestCalculateRPre(unittest.TestCase):
	def test_normal_usage_1(self):
		resistivity = 983.
		length = 100.
		width = 10.
		active_radius = 7.
		
		r_pre = calculate_r_pre(length, width, active_radius,
				resistivity)

		should_be = 3019.21428571428571428571

		self.assertAlmostEqual(r_pre, should_be)


class TestCalculateROut(unittest.TestCase):
	def test_normal_usage_1(self):
		resistivity = 983.
		length = 10.
		width = 100.
		active_radius = 7.
		
		r_out = calculate_r_out(length, width, active_radius,
				resistivity)

		should_be = 228.604651162790697

		self.assertAlmostEqual(r_out, should_be)


class TestCalculateDV(unittest.TestCase):
	def test_normal_usage_1(self):
		r = 10.
		v_macro = 20.
		molecule_diameter = 2.
		pore_height = 25.
		pore_diameter = 3.
		concentration_ratio = 100.

		should_be = 0.55500909645

		dv = calculate_dv(r, v_macro, molecule_diameter, pore_height,
				pore_diameter, concentration_ratio)

		self.assertAlmostEqual(dv, should_be)

	def test_normal_usage_2(self):
		r = 99.
		v_macro = 15.
		molecule_diameter = 6.
		pore_height = 50.
		pore_diameter = 10.
		concentration_ratio = 1000.

		should_be = 0.51962267

		dv = calculate_dv(r, v_macro, molecule_diameter, pore_height,
				pore_diameter, concentration_ratio)

		self.assertAlmostEqual(dv, should_be)


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
