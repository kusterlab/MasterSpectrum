
from mgf_annotation.mgfAnnotaterDistiller import MgfAnnotaterDistiller
from nose.tools import *


class TestMgfAnnotaterDistiller(object):
    def test_loading_improve(self):
        path = "tests/data/pick_pairs_R4.mgf"
        output = "tests/data/temp/test.csv"
        mad = MgfAnnotaterDistiller(path, output)
        mad.load_improved_csv("tests/data/tests_20_improved.csv")
        assert_equal(mad.ids_to_be_referenced[39], 38)
        assert_equal(mad.ids_to_be_referenced[1199], 1198)
        assert_equal(mad.ids_to_be_referenced[1189], 1188)

    def test_test(self):

        path = "tests/data/distiller.mgf"
        output = "tests/data/temp/test.csv"
        mad = MgfAnnotaterDistiller(path, output)
        mad.load_improved_csv("tests/data/tests_20_improved.csv")
        mad.load_distiller_mgf()
        assert_equal(mad.ids_to_be_referenced[39], 38)
        assert_equal(mad.ids_to_be_referenced[1199], 1198)
