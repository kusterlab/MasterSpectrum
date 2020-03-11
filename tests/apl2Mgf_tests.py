
from nose.tools import *
from mgf_filter.apl2Mgf import Apl2Mgf
import shutil
import os


class TestApl2Mgf(object):
    @classmethod
    def teardown_class(cls):
        shutil.rmtree("tests/data/temp")

    @classmethod
    def setup_class(cls):
        if not os.path.exists("tests/data/temp"):
            os.makedirs("tests/data/temp")

    def test_finding_apl(self):
        apl_folder = "tests/data/"
        output_mgf = "tests/data/temp/tadaaa.mgf"
        a2m = Apl2Mgf(apl_folder,
                      output_mgf)
        a2m.find_all_apl()
        assert_equal(a2m.apl_files, ['allSpectra.HCD.FTMS.iso_3_shorted.apl', 'allSpectra.HCD.FTMS.iso_4._shorted.apl'])

    def test_read_apl(self):
        apl_folder = "tests/data/"
        output_mgf = "tests/data/temp/tadaaa.mgf"
        a2m = Apl2Mgf(apl_folder,
                      output_mgf)
        a2m.find_all_apl()
        a2m.load_apls()
        print(a2m.apl_data.keys())
        assert_equal(a2m.apl_data["RawFile: 01509_C01_P015226_S00_A00_CE30 Index: 18670 Precursor: 0 _multi___tobi__3"].mz[112.08698], ["2358.949"])
        a2m.writeMgf()
