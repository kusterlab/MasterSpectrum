from mgf_filter.deltaExtractor import DeltaExtractor
from nose.tools import *


class TestDeltaExtractor(object):
    def test_append_different_Bins(self):
        dEx = DeltaExtractor("tests/data/cetsa_101287_A01_P013190_S00_N01_R1_TMT10_3spectra_shorted.mgf")
        dEx.createDeltaMasterSpectrum(0)
        print(dEx.masterSpectrum.spectrum)
        assert_true(359.0 in dEx.masterSpectrum.spectrum[0])
        assert_true(15 in dEx.masterSpectrum.spectrum[0])
