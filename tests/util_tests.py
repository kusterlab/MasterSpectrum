from mgf_filter.util import calculateRelativeIntensity
from mgf_filter.util import calculate_Delta_Fixed
from mgf_filter.util import calculate_Delta_based_MZ
from mgf_filter.util import calculate_Delta_by_ppm
from mgf_filter.util import calculatePrecursor
from nose.tools import *
import math


class TestUtil(object):
    def test_fixed_delta(self):
        fixed_0_5 = calculate_Delta_Fixed(0.5)
        assert_equal(0.5, fixed_0_5(4))
        assert_equal(0.5, fixed_0_5(7))

    def test_rel_delta(self):
        assert_equal(100 / 56885.29308438425, calculate_Delta_based_MZ(100))

    def test_rel_double_delta(self):
        assert_equal(2 * 100 / 56885.29308438425, calculate_Delta_based_MZ(100, numSigma=2))

    def test_rel_Intensity(self):
        assert_equal([1, 0.5, 0.1], calculateRelativeIntensity([10, 5, 1]))

    def test_ppm_delta(self):
        delta_ppm_20 = calculate_Delta_by_ppm(20)
        assert_equal(20 * (100 / (math.pow(10, 6))), delta_ppm_20(100))

    def test_calculatePrecursor(self):
        precursor = calculatePrecursor(100, 1)
        assert_equal(precursor, 100)
        precursor = calculatePrecursor(100, 2)
        assert_equal(precursor, 200 - 1.00782503)
