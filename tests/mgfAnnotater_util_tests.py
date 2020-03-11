
from mgf_annotation.util import gen_allowed_mass_diff_with_sign
from mgf_annotation.util import parse_scan_id
from nose.tools import *


class TestAnnotationUtil(object):
    def test_allowed_mass_sign(self):
        data = [x for x in gen_allowed_mass_diff_with_sign(3, 3)]
        assert_equal(len(data), 18)

    def test_splitting(self):
        test = "\"controllerType=0 controllerNumber=1 scan=1316\""
        assert_equal(parse_scan_id(test), '1316')
