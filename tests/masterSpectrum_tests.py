from nose.tools import *
from mgf_filter.peak import Peak
from mgf_filter.masterSpectrum import MasterSpectrum
from mgf_filter.util import calculate_Delta_based_MZ
from mgf_filter.util import calculate_Delta_Fixed
import os
import shutil


class TestSpectrum(object):
    def test_init(self):
        ms = MasterSpectrum()
        assert_equal(ms.appended, 0)
        assert_equal(ms.merged, 0)

    def test_append_different_Bins(self):
        ms = MasterSpectrum()

        # create first MP
        p1 = Peak(100, 0.5, calculate_Delta_based_MZ)
        ms.add(p1)

        print(p1)
        assert_equal(ms.appended, 1)
        assert_equal(ms.merged, 0)

        p2 = Peak(100.1, 0.5, calculate_Delta_based_MZ)
        print(p2)
        ms.add(p2)

        assert_equal(ms.appended, 2)
        assert_equal(ms.merged, 0)

        assert_equal(len(ms.spectrum[0][100]), 1)
        assert_equal(len(ms.spectrum[0][101]), 1)
        assert_equal(len(ms.spectrum[0]), 2)

    def test_append_same_bin(self):
        ms = MasterSpectrum()

        # create first MP
        p1 = Peak(100.2, 0.5, calculate_Delta_based_MZ)
        ms.add(p1)

        assert_equal(ms.appended, 1)
        assert_equal(ms.merged, 0)

        p2 = Peak(100.7, 0.5, calculate_Delta_based_MZ)
        ms.add(p2)

        assert_equal(ms.appended, 2)
        assert_equal(ms.merged, 0)

        assert_equal(len(ms.spectrum[0][101]), 2)
        assert_equal(len(ms.spectrum[0]), 1)

    def test_binary_search_1_exist_1add_not_overlapping(self):

        ms = MasterSpectrum()

        # create first MP
        p1 = Peak(100.5, 0.5, calculate_Delta_based_MZ)
        ms.add(p1)

        p2 = Peak(100.8, 0.5, calculate_Delta_based_MZ)
        assert_equal(p2.key(), p1.key())
        assert_false(ms.spectrum[0][101][0].isInside(p2))

        size_bin = len(ms.spectrum[0][p1.key()])
        assert_equal(size_bin, 1)

    def test_binary_search_1_exist_1add_overlapping(self):

        ms = MasterSpectrum()

        # create first MP
        p1 = Peak(100.5, 0.5, calculate_Delta_based_MZ)
        ms.add(p1)

        p2 = Peak(100.50000000001, 0.5, calculate_Delta_based_MZ)
        assert_true(ms.spectrum[0][101][0].isInside(p2))
        assert_equal(p2.key(), p1.key())

        size_bin = len(ms.spectrum[0][p1.key()])
        assert_equal(size_bin, 1)

    def test_binary_search_2_exist_not_overlapping_1add_in_upper_window(self):
        ms = MasterSpectrum()
        # create first MP
        p1 = Peak(100.5, 0.5, calculate_Delta_based_MZ)
        ms.add(p1)

        p2 = Peak(100.9, 0.5, calculate_Delta_based_MZ)
        ms.add(p2)

        p3 = Peak(100.90000001, 0.5, calculate_Delta_based_MZ)

        assert_true(ms.spectrum[0][101][1].isInside(p3))

        size_bin = len(ms.spectrum[0][p1.key()])
        assert_equal(size_bin, 2)

    def test_binary_search_2_exist_not_overlapping_1add_in_lower_window(self):
        ms = MasterSpectrum()
        # create first MP
        p1 = Peak(100.5, 0.5, calculate_Delta_based_MZ)
        ms.add(p1)

        p2 = Peak(100.9, 0.5, calculate_Delta_based_MZ)
        ms.add(p2)

        p3 = Peak(100.4999999, 0.5, calculate_Delta_based_MZ)
        print(ms)
        assert_true(ms.spectrum[0][101][0].isInside(p3))

        size_bin = len(ms.spectrum[0][p1.key()])
        assert_equal(size_bin, 2)

    def test_binary_search_2_exist_not_overlapping_1add_in_both_windows(self):
        ms = MasterSpectrum()

        p1 = Peak(100.78, 0.5, calculate_Delta_Fixed(0.1))
        ms.add(p1)

        p2 = Peak(100.9, 0.7, calculate_Delta_Fixed(0.1))
        ms.add(p2)

        assert_equal(len(ms.spectrum[0][101]), 2)
        assert_equal(ms.spectrum[0][101][0].intensity, 0.5)
        assert_equal(ms.spectrum[0][101][1].intensity, 0.7)
        p3 = Peak(100.85, 0.5, calculate_Delta_Fixed(0.1))
        ms.add(p3)
        print(ms)
        assert_equal(len(ms.spectrum[0][101]), 1)
        assert_equal(ms.spectrum[0][101][0].intensity, 1.7)
        assert_equal(ms.multimerged, 1)
        print(ms)
        assert_equal(ms.spectrum[0][101][0].mz_origin, 100.9)
        assert_equal(ms.spectrum[0][101][0].counts, 3)
        print(ms)
        p_pull_left = Peak(100.76, 1, calculate_Delta_Fixed(0.1))
        ms.add(p_pull_left)
        p_pull_left = Peak(100.76, 4, calculate_Delta_Fixed(0.1))
        ms.add(p_pull_left)
        p4 = Peak(100.89, 0.5, calculate_Delta_Fixed(0.1))
        ms.add(p4)
        print(ms)

    def test_2_peaks_made_problems_above_boarder_bigger_bin_first_merged_moves_left(self):
        """
        p1 in bin1
        p2 is inside p1
        p2 should move bin1
        checks
        p1 in p2
        p3 = p1 + p2
        p3 moves bin0
        """
        ms = MasterSpectrum()
        p1 = Peak(2823.396729, 0.3, calculate_Delta_based_MZ)
        ms.add(p1)
        p2 = Peak(2823.174072, 0.5, calculate_Delta_based_MZ)

        print(p1)
        print(p2)
        assert_true(ms.spectrum[0][2824][0].isInside(p2))
        ms.add(p2)
        assert_equal(len(ms.spectrum[0][2823]), 1)
        assert_equal(ms.spectrum[0][2823][0].intensity, 0.8)

    def test_2_peaks_made_problems_above_boarder_bigger_bin_first_merged_stays_right(self):

        """
        p1 in bin1
        p2 is inside p1
        p2 should move to bin
        checks
        p1 in p2
        p3 = p1 + p2
        p3 moves bin1
        """
        ms = MasterSpectrum()
        p1 = Peak(2823.396729, 0.9, calculate_Delta_based_MZ)
        ms.add(p1)
        p2 = Peak(2823.174072, 0.1, calculate_Delta_based_MZ)

        assert_true(ms.spectrum[0][2824][0].isInside(p2))
        ms.add(p2)
        assert_equal(len(ms.spectrum[0][2824]), 1)
        assert_equal(ms.spectrum[0][2824][0].intensity, 1)
        assert_true(not(2823 in ms.spectrum[0]))
        print(ms)

    def test_3_peaks_2_in_a_bin_3rd_moves_first_left(self):
        """
        p1 in bin1
        p2 in bin1
        p2 > p1
        p3 is inside p1
        p3 should move to bin0
        checks
        p4 = p1+p3
        p4 moves in b0
        """
        ms = MasterSpectrum()
        p1 = Peak(2823.396729, 0.3, calculate_Delta_based_MZ)
        ms.add(p1)
        p2 = Peak(2823.774072, 0.5, calculate_Delta_based_MZ)
        ms.add(p2)

        assert_false(ms.spectrum[0][2824][0].isInside(p2))

        assert_equal(len(ms.spectrum[0][2824]), 2)

        p3 = Peak(2823.174072, 0.5, calculate_Delta_based_MZ)
        ms.add(p3)

        assert_equal(ms.spectrum[0][2823][0].intensity, 0.8)
        assert_equal(len(ms.spectrum[0][2824]), 1)

    def test_3_peaks_2_in_a_bin_3rd_lets_stay(self):
        """
        p1 in bin1
        p2 in bin1
        p2 > p1
        p3 is inside p1
        p3 should move to bin0
        checks
        p4 = p1+p3
        p4 stays in b1
        """

        ms = MasterSpectrum()
        p1 = Peak(2823.396729, 0.9, calculate_Delta_based_MZ)
        ms.add(p1)
        p2 = Peak(2823.774072, 0.5, calculate_Delta_based_MZ)
        ms.add(p2)

        assert_false(ms.spectrum[0][2824][0].isInside(p2))

        assert_equal(len(ms.spectrum[0][2824]), 2)

        p3 = Peak(2823.174072, 0.1, calculate_Delta_based_MZ)
        ms.add(p3)

        assert_equal(ms.spectrum[0][2824][0].intensity, 1)
        assert_equal(len(ms.spectrum[0][2824]), 2)
        assert_true(not(2823 in ms.spectrum[0]))

    def test_2_peaks_merged_moved_left_insert_3rd_peak_right(self):
        """
        p1 in bin1
        p2 is inside p1
        p2 should move bin1
        checks
        p1 in p2
        p3 = p1 + p2
        p3 moves bin0
        p4 in bin1
        """

        ms = MasterSpectrum()
        p1 = Peak(2823.396729, 0.3, calculate_Delta_based_MZ)
        ms.add(p1)
        p2 = Peak(2823.174072, 0.5, calculate_Delta_based_MZ)

        assert_true(ms.spectrum[0][2824][0].isInside(p2))
        ms.add(p2)

        p3 = Peak(2823.596729, 0.3, calculate_Delta_based_MZ)
        ms.add(p3)
        print(ms)
        assert_equal(len(ms.spectrum[0][2823]), 1)
        assert_equal(len(ms.spectrum[0]), 2)
        assert_equal(len(ms.spectrum[0][2824]), 1)

    def test_for_refactored_append_class(self):
        ms = MasterSpectrum()

        p3 = Peak(2566.308, 0.02237079, calculate_Delta_based_MZ)
        ms.add(p3)

        p1 = Peak(2567.371190, 0.0525852066, calculate_Delta_based_MZ)
        print(p1)
        ms.add(p1)

        p2 = Peak(2567.199755, 0.0577020220, calculate_Delta_based_MZ)
        print(p2)

        assert_true(ms.spectrum[0][2568][0].isInside(p2))

        ms.add(p2)
        p3 = Peak(2566.308, 0.02237079, calculate_Delta_based_MZ)
        ms.add(p3)
        assert_equal(ms.spectrum[0][2568][0].intensity, 0.0525852066 + 0.0577020220)
        assert_equal(len(ms.spectrum[0][2567]), 1)
        assert_equal(len(ms.spectrum[0][2568]), 1)

    def test_strange_missing_left_merges(self):
        ms = MasterSpectrum()
        p1 = Peak(128.002106, 0.11, calculate_Delta_based_MZ)
        ms.add(p1)

        p2 = Peak(128.00296464, 0.39, calculate_Delta_based_MZ)
        ms.add(p2)

        assert_equal(ms.spectrum[0][129][0].counts, 2)
        ms = MasterSpectrum()
        p1 = Peak(128.002106, 0.11, calculate_Delta_based_MZ)

        p2 = Peak(128.00296464, 0.39, calculate_Delta_based_MZ)
        ms.add(p2)
        ms.add(p1)

        assert_equal(ms.spectrum[0][129][0].counts, 2)
        assert_equal(len(ms.spectrum[0][129]), 1)

    def test_loading_mgf_same_charges(self):
        ms = MasterSpectrum()
        ms.load_from_mgf("tests/data/cetsa_101287_A01_P013190_S00_N01_R1_TMT10_2_charge_states.mgf", ignoreCharges=False)
        assert_equal(len(ms.spectrum), 3)

    @classmethod
    def teardown_class(cls):
        shutil.rmtree("tests/data/temp")

    @classmethod
    def setup_class(cls):
        if not os.path.exists("tests/data/temp"):
            os.makedirs("tests/data/temp")

    def test_load_mgf_binning_reloading(self):
        ms = MasterSpectrum()
        ms.load_from_mgf("tests/data/cetsa_101287_A01_P013190_S00_N01_R1_TMT10_3spectra_shorted.mgf", ignoreCharges=False)
        print("file is loaded")
        ms.export_to_csv("tests/data/temp/cetsa_101287_A01_P013190_S00_N01_R1_TMT10_3spectra_shorted.csv")
        print("file is written")
        ms2 = MasterSpectrum()
        ms2.load_from_csv("tests/data/temp/cetsa_101287_A01_P013190_S00_N01_R1_TMT10_3spectra_shorted.csv")
        print("file is loaded from csv")
        ms2.multimerged = ms.multimerged
        ms2.merged = ms.merged
        ms2.appended = ms.appended
        ms2.export_to_csv("tests/data/temp/cetsa_101287_A01_P013190_S00_N01_R1_TMT10_3spectra_shorted_2.csv")
        assert_true(ms == ms2)

    def test_compare_two_spectra(self):
        ms = MasterSpectrum()
        ms.load_from_csv("tests/data/df_charged.csv")
        ms2 = MasterSpectrum()
        ms2.load_from_csv("tests/data/df_1_brain_charged.csv")

        ms3 = ms.compare_other_ms(ms2)
        assert_equal(ms3.spectrum['0'][90][0].counts_ratio, -0.5)
        assert_equal(ms3.spectrum['0'][90][0].counts_ratio, -0.5)

    def test_loading_mgf_to_one_charge(self):
        ms = MasterSpectrum()
        ms.load_from_mgf("tests/data/cetsa_101287_A01_P013190_S00_N01_R1_TMT10_2_charge_states.mgf", ignoreCharges=True)
        assert_equal(len(ms.spectrum), 1)

    def test_reactToPeaksMovingLeftBinConnectingRightBin(self):
        ms = MasterSpectrum()

        p1 = Peak(100.5, 0.005, calculate_Delta_Fixed(0.2))
        ms.add(p1)
        p1 = Peak(100.31, 1, calculate_Delta_Fixed(0.2))
        ms.add(p1)
        print(ms)
        p1 = Peak(100.1, 1, calculate_Delta_Fixed(0.3))
        ms.add(p1)
        print(ms)
        p1 = Peak(100.2, 1, calculate_Delta_Fixed(0.3))
        ms.add(p1)

        print(ms)
        assert_equal(len(ms.spectrum[0][100]), 1)
        assert_equal(ms.spectrum[0][100][0].counts, 4)
