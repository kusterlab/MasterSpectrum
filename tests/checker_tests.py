

from nose.tools import *
import os
import shutil
from mgf_checker.checker import Checker


class TestChecker(object):
    _multiprocess_shared_ = True

    def test_parsemzid(self):
        c = Checker("tests/data/F082488.mzid")
        c.parse_mz_id()
        assert_equal(c.identifications[0][1].ion_series_ary[0].fragtype, "frag: b ion")
        assert_equal(c.identifications[0][1].ion_series_ary[0].charge, 1)
        assert_equal(c.identifications[0][1].ion_series_ary[0].mz_ary, 129.10202)
        """
              <IonType index="1" charge="1">
                <FragmentArray values="129.10202" measure_ref="m_mz" />
                <FragmentArray values="3745" measure_ref="m_intensity" />
                <FragmentArray values="0.0362" measure_ref="m_error" />
                <cvParam cvRef="PSI-MS" accession="MS:1001224" name="frag: b ion" />
              </IonType>
        """

    def test_readmgf(self):
        c = Checker("tests/data/F082488.mzid")
        c.read_enhanced_spectrum("/home/tobiass/goto/master/data/HeLa/2_nd_measurenment/01533_B12_P015940_B00_A00_R3_20_improved.mgf")
        assert_equal(c.mgf_reads['39'].mz_ary[0], 136.06149289999999)

        ms = c.mgf_reads['39'].request_ms()

        assert_equal(len(ms.spectrum[0][351]), 1)

    def test_compare_mzid_mgf(self):
        # c = Checker("/home/tobiass/Desktop/F083063.mzid")
        # c.parse_mz_id()
        # c.read_enhanced_spectrum("/home/tobiass/goto/master/data/HeLa/2_nd_measurenment/01533_B12_P015940_B00_A00_R3_20_improved.mgf")
        # c.analyse_mzid_vs_mgf()
        # assert_equal(1, 2)
        pass

    def test_get_peptide_info(self):
        c = Checker("tests/data/F082488.mzid")
        c.parse_mz_id_peptide_ref()
        assert_equal(c.peptide_evidence['peptide_6190_1'].peptide_sequence, 'CCTESLVNRRPCFSALTPDETYVPK')
        assert_equal(c.peptide_evidence['peptide_6190_1'].modification[0].name, 'TMT6plex')
        assert_equal(c.peptide_evidence['peptide_6190_1'].modification[25].name, 'TMT6plex')
        assert_equal(c.peptide_evidence['peptide_6190_1'].modification[2].name, 'Carbamidomethyl')

    def test_info_about_peptide_tag_amount(self):
        """
        peptide_6180_1 has two TMT and is therefore more interesting
        sequence is 23 amino acid long
        """
        c = Checker("tests/data/F082488.mzid")
        c.parse_mz_id()
        c.parse_mz_id_peptide_ref()
        assert_equal(c.peptide_evidence['peptide_6180_1'].peptide_sequence, 'CCTKPESERMPCTEDYLSLILNR')
        b_tmt, y_tmt = c.peptide_evidence['peptide_6180_1'].get_annotated_positions()
        assert_equal(b_tmt[4], 2)
        assert_equal(b_tmt[3], 1)
        assert_equal(b_tmt[1], 1)
        assert_equal(b_tmt[23], 2)
        assert_equal(y_tmt[23], 2)
        assert_equal(y_tmt[21], 1)
        assert_equal(y_tmt[20], 1)
        assert_equal(y_tmt[19], 0)
        assert_equal(y_tmt[3], 0)

    def test_read_score_file(self):
        c = Checker("tests/data/F082488.mzid")
        allowed_ids = [8358]
        spectra = c.read_score_file("/home/tobiass/goto/master/data/HeLa/2_nd_measurenment/01533_B12_P015940_B00_A00_R3_20_improved.csv", allowed_ids)
        assert_equal(8358 in spectra, True)
        assert_equal(len(spectra[8358]), 2)
        assert_equal(len(spectra[8358]['mz']), 129)

    def test_compare_mz_id_score_file(self):
        # c = Checker("/home/tobiass/Desktop/F083270.mzid")
        # c.parse_mz_id()
        # c.analyze_mzid_id_vs_score_file(path_score="/home/tobiass/goto/master/data/HeLa/2_nd_measurenment/01533_B12_P015940_B00_A00_R3_20_improved.csv",
        #                                output_path="/home/tobiass/Desktop/out.csv")
        assert_equal(1, 1)

    @classmethod
    def teardown_class(cls):
        shutil.rmtree("tests/data/temp")

    @classmethod
    def setup_class(cls):
        if not os.path.exists("tests/data/temp"):
            os.makedirs("tests/data/temp")
