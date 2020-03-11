from nose.tools import *
from mgf_filter.masterSpectrum import MasterSpectrum
import os
import shutil
from mgf_filter.mgfPatcher import MgfPatcher
from mgf_filter.util import calculate_Delta_by_ppm


class TestMgfPatcher(object):

    @classmethod
    def teardown_class(cls):
        shutil.rmtree("tests/data/temp")

    @classmethod
    def setup_class(cls):
        if not os.path.exists("tests/data/temp"):
            os.makedirs("tests/data/temp")

    def test_init(self):
        mgfP = MgfPatcher(delta_func=calculate_Delta_by_ppm(20))
        a = mgfP.delta_function(1000000)
        assert_equal(a, 20)

    def test_load_exclusion_list(self):
        mgfP = MgfPatcher(delta_func=calculate_Delta_by_ppm(20))
        mgfP.readExclusionList('tests/data/exclusionList.txt')
        assert_equal(mgfP.exclusionSpectrum.spectrum[0][177][0].mz, 176.309174)
        assert_equal(mgfP.precursorDeltas, [229])

    def test_patch_mgf(self):
        ms = MasterSpectrum()
        ms.load_from_mgf('tests/data/cetsa_101287_A01_P013190_S00_N01_R1_TMT10_3spectra_shorted.mgf', ignoreCharges=True)
        assert_equal(177 in ms.spectrum[0].keys(), True)

        mgfP = MgfPatcher(delta_func=calculate_Delta_by_ppm(20))
        mgfP.readExclusionList('tests/data/exclusionList.txt')
        mgfP.patchMgf('tests/data/cetsa_101287_A01_P013190_S00_N01_R1_TMT10_3spectra_shorted.mgf', 'tests/data/temp/cetsa_101287_A01_P013190_S00_N01_R1_TMT10_3spectra_shorted_patched.mgf')

        ms = MasterSpectrum()
        ms.load_from_mgf('tests/data/temp/cetsa_101287_A01_P013190_S00_N01_R1_TMT10_3spectra_shorted_patched.mgf', ignoreCharges=True)
        assert_equal(177 in ms.spectrum[0].keys(), False)
        assert_equal(433 in ms.spectrum[0].keys(), False)
