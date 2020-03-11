from mgf_annotation.mgfAnnotater import MgfAnnotater
from mgf_filter.masterSpectrum import MasterSpectrum
from nose.tools import *
import os
import shutil


class TestMgfAnnotater(object):
    def test_loading(self):
        path = "tests/data/pick_pairs_R4.mgf"
        output = "tests/data/temp/test.csv"
        ma = MgfAnnotater(path, output)
        ma.load_msconvert_mgf()
        assert_equal(ma.references[0].peak_list_1[0], 100.5521622)
        assert_equal(ma.references[0].peak_list_2[0], float(103.2448959))

    def test_export_compared_spectra(self):
        path = "tests/data/pick_pairs_R4.mgf"
        output = "tests/data/temp/test.csv"
        ma = MgfAnnotater(path, output)
        ma.load_msconvert_mgf()
        ma.export_annotated_spectra_to_csv()

    def test_weird_decharging(self):
        path = "tests/data/weird.mgf"
        output = "tests/data/temp/test.csv"
        ma = MgfAnnotater(path, output)
        ma.load_msconvert_mgf()
        ma.export_annotated_spectra_to_csv()

    def test_mgf(self):
        path = "tests/data/weird.mgf"
        output = "tests/data/temp/test.csv"
        output_mgf = "tests/data/temp/test.mgf"
        ma = MgfAnnotater(path, output)
        ma.load_msconvert_mgf()
        ma.export_annotated_spectra_to_mgf(output_mgf)

    def test_weird_summing_mgf(self):
        path = "tests/data/weird_sum.mgf"
        output = "/home/tobiass/weird_sum.csv"
        output_mgf = "tests/data/temp/weird_mgf.mgf"
        ma = MgfAnnotater(path, output)
        ma.load_msconvert_mgf()
        ma.export_annotated_spectra_to_mgf(output_mgf)

        ms = MasterSpectrum()
        ms.load_from_mgf("tests/data/temp/weird_mgf.mgf", ignoreCharges=True)
        assert_equal(round(ms.spectrum[0][179][0].intensity, 4), round(0.0025, 4))

    @classmethod
    def teardown_class(cls):
        shutil.rmtree("tests/data/temp")

    @classmethod
    def setup_class(cls):
        if not os.path.exists("tests/data/temp"):
            os.makedirs("tests/data/temp")
