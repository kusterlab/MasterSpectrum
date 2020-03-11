
from nose.tools import *
from mgf_filter.masterPeak import MasterPeak
from mgf_filter.peak import Peak
import math
from mgf_filter.util import calculate_Delta_based_MZ


class TestMasterPeak(object):
    def test_init(self):
        p = Peak(100, 0.5, calculate_Delta_based_MZ)
        mp = MasterPeak(p)
        assert_equal(mp.mz, 100)
        assert_equal(mp.intensity, 0.5)
        assert_not_equal(mp.left, 0)
        assert_not_equal(mp.right, 0)
        assert_equal(mp.counts, 1)

    def test_add(self):
        p1 = Peak(100, 0.5, calculate_Delta_based_MZ)
        p2 = Peak(110, 0.4, calculate_Delta_based_MZ)

        mp = MasterPeak(p1)
        mp.add(p2)
        assert_equal(mp.intensity, 0.9)
        assert_equal(mp.mz, (100 * 0.5 + 110 * 0.4) / (0.5 + 0.4))
        assert_equal(mp.counts, 2)

        p3 = Peak(120, 0.6, calculate_Delta_based_MZ)
        mp.add(p3)
        assert_equal(mp.mz, (100 * 0.5 + 110 * 0.4 + 120 * 0.6) / (0.5 + 0.4 + 0.6))
        assert_equal(mp.intensity, 0.9 + 0.6)
        assert_equal(mp.counts, 3)

    def test_update(self):
        p1 = Peak(100, 0.5, calculate_Delta_based_MZ)
        mp = MasterPeak(p1)

        assert_equal(mp.left, mp.mz - mp.delta)
        assert_equal(mp.right, mp.mz + mp.delta)
        assert_equal(mp.delta, 100 / 56885.29308438425)

        p2 = Peak(110, 0.4, calculate_Delta_based_MZ)
        mp.add(p2)

        mz2 = (100 * 0.5 + 110 * 0.4) / 0.9
        delta = mz2 / (math.pow(10, 5.847 + math.log10(mz2) * (-0.546)))
        assert_equal(mp.delta, delta)
        assert_equal(mp.left, mz2 - delta)
        assert_equal(mp.right, mz2 + delta)

    def test_isInside(self):
        p1 = Peak(100, 0.5, calculate_Delta_based_MZ)
        p2 = Peak(101, 0.3, calculate_Delta_based_MZ)
        p3 = Peak(100 + 0.001, 0.2, calculate_Delta_based_MZ)
        mp = MasterPeak(p1)

        assert_false(mp.isInside(p2))
        assert_true(mp.isInside(p3))

    def test_isInsideMz(self):
        p1 = Peak(100, 0.5, calculate_Delta_based_MZ)
        mp = MasterPeak(p1)

        assert_false(mp.isInsideMz(102))
        assert_true(mp.isInsideMz(100))

    def test_smaller(self):
        #  tests if right border is smaller than peak
        p1 = Peak(100, 0.5, calculate_Delta_based_MZ)
        mp = MasterPeak(p1)

        p2 = Peak(90, 0.5, calculate_Delta_based_MZ)
        p3 = Peak(110, 0.5, calculate_Delta_based_MZ)

        assert_true(mp.smaller(p3))
        assert_false(mp.smaller(p2))

    def test_greater(self):
        #  tests if left border is greater than peak
        p1 = Peak(100, 0.5, calculate_Delta_based_MZ)
        mp = MasterPeak(p1)

        p2 = Peak(90, 0.5, calculate_Delta_based_MZ)
        p3 = Peak(110, 0.5, calculate_Delta_based_MZ)

        assert_true(mp.greater(p2))
        assert_false(mp.greater(p3))

    def test_ratio(self):

        p1 = Peak(100, 0.5, calculate_Delta_based_MZ)
        mp = MasterPeak(p1)

        p2 = Peak(100, 0.5, calculate_Delta_based_MZ)
        mp2 = MasterPeak(p2)

        mp.recalculate_ratio(mp2)
        assert_equal(mp.rel_intensity_ratio, 1)
        assert_equal(mp.counts_ratio, 1)
