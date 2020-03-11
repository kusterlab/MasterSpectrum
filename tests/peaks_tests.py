
from nose.tools import *
from mgf_filter.peak import Peak
from mgf_filter.util import calculate_Delta_based_MZ

class TestPeak(object):
    def test_init(self):
        p=Peak(100,0.5, calculate_Delta_based_MZ)
        assert_equal(p.mz,100)
        assert_equal(p.intensity,0.5)
        assert_not_equal(p.left,0)
    
    def test_delta(self):
        p=Peak(100,0.5, calculate_Delta_based_MZ)
        assert_equal(p.delta,100/56885.29308438425)

    def test_update(self):
        p=Peak(100,0.5, calculate_Delta_based_MZ)
        assert_equal(p.left,100-100/56885.29308438425)
    
    def test_key(self):
        p=Peak(100,0.5, calculate_Delta_based_MZ)
        assert_not_equal(p.key(),99)
        assert_equal(p.key(),100)
