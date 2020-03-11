from mgf_filter.perfectSpectra import calculateDoubleCharged
from nose.tools import *


class TestPerfectSpectra(object):
    def test_calculateDoubleCharged(self):
        assert_equal(calculateDoubleCharged(100), 50.503912516035)
