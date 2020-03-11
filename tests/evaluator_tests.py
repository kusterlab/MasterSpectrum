
from nose.tools import *
import os
import shutil
from mgf_recalibration.evaluator import Evaluator


class TestEvaluator(object):
    def test_parsemzid(self):
        e = Evaluator("tests/data/F082488.mzid")
        e.parse_mz_id()
        assert_equal(e.scan1.data['4872'][1][0].indice[0], 1)
        assert_equal(e.scan1.data['4872'][1][0].indice[1], 2)
        assert_equal(e.scan1.data['4872'][1][0].indice[2], 3)

    @classmethod
    def teardown_class(cls):
        shutil.rmtree("tests/data/temp")

    @classmethod
    def setup_class(cls):
        if not os.path.exists("tests/data/temp"):
            os.makedirs("tests/data/temp")
