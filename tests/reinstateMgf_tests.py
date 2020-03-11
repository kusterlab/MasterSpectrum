
from nose.tools import *
from mgf_filter.reinstateMgf import ReinstateMgf


class TestReinstateMgf(object):
    def test_init(self):
        patchedV = "tests/data/cetsa_101287_A01_P013190_S00_N01_R1_TMT10_3spectra_shorted_filtered.mgf"
        originalV = "tests/data/cetsa_101287_A01_P013190_S00_N01_R1_TMT10_3spectra_shorted.mgf"
        rM = ReinstateMgf(patchedV,
                          originalV)
        rM.loadPatched()
        title = "msmsid:F000019,rt:0.099,survey:S000018,parent:418.89,AnalTime:110.00,Activation:HCD"
        d = rM.patched_Data[title].mz
        assert_equal(d[129.102417][0], '661.950378')
        assert_equal(d[134.956573][0], '2347.313965')
        assert_equal(len(d), 2)
        assert_true(134.614151 not in d)
