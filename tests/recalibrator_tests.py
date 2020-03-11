
from nose.tools import *
from mgf_filter.masterSpectrum import MasterSpectrum
from mgf_recalibration.recalibrator import Recalibrator
import os
import shutil
from mgf_recalibration.util import calculate_tag_tmt10
from mgf_recalibration.util import calculate_ppm_shift
from mgf_recalibration.util import calculate_da_shift


class TestRecalibrator(object):
    def test_tmt(self):
        mass = calculate_tag_tmt10()
        assert_equal(round(mass, 5), 230.17021)

    def test_recalibration(self):
        path = "tests/data/weird.mgf"
        re = Recalibrator(path=path, ppm=20)
        re.load_recalibrate()

    def test_ppm_shift(self):
        assert_equal(calculate_ppm_shift(1, 1000000), 1)
        assert_equal(calculate_ppm_shift(-1, 1000000), -1)

    def test_da_shift(self):
        assert_equal(calculate_da_shift(1, 1000000), 1)
        assert_equal(calculate_da_shift(-1, 1000000), -1)

    def test_recalibration_export(self):
        path = "tests/data/weird.mgf"
        re = Recalibrator(path=path, ppm=20, file_out="tests/data/temp/out.mgf")
        re.load_recalibrate()

        ms = MasterSpectrum()
        ms.load_from_mgf("tests/data/temp/out.mgf", ignoreCharges=True)
        assert_equal(round(ms.spectrum[0][231][1].mz, 7), round(calculate_tag_tmt10(), 7))

    @classmethod
    def teardown_class(cls):
        shutil.rmtree("tests/data/temp")

    @classmethod
    def setup_class(cls):
        if not os.path.exists("tests/data/temp"):
            os.makedirs("tests/data/temp")
