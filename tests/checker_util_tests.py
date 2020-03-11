
from nose.tools import *
from mgf_checker.util import IonSeriesUsed


class TestChecker_util(object):
    _multiprocess_shared_ = True

    def test_IonSeriesUsed(self):
        r = IonSeriesUsed('0000021020000000000')
        assert_equal(r['b'][1], False)
        assert_equal(r['b'][2], True)
        assert_equal(r['b'][1], False)
        assert_equal(r['b'][2], True)
