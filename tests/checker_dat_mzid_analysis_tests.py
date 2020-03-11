
from nose.tools import *
from mgf_checker.checker_dat_mzid_analysis import MS_parser_result
from mgf_checker.checker_dat_mzid_analysis import MZID_Comparison_result
from mgf_checker.checker_dat_mzid_analysis import Comparer


class TestChecker_Data_Loader(object):
    _multiprocess_shared_ = True

    def test_readMZID_data_result_file(self):
        m = MZID_Comparison_result("tests/data/out.csv")
        m.read_file()
        assert_equal(m.data[24647][1].has_Series('B'), True)
        assert_equal(m.data[24647][1].has_Series('Y'), True)

        assert_equal(m.data[9852][1].BSerie.byCharge[2][1]['Peak'], 172.153061)

    def test_MS_parser_result_file(self):
        m = MS_parser_result("tests/data/F083270_score.txt")
        m.read_file()
        assert_equal(m.data[24643][1].dScore, 9.06)

    def test_ms_parser_vs_mzid_out_result(self):

        mzid_comparison = MZID_Comparison_result("tests/data/out.csv")
        mzid_comparison.read_file()
        ms_parser_result = MS_parser_result("tests/data/F083270_score.txt")
        ms_parser_result.read_file()

    def test_more_stuff(self):
        c = Comparer(stPath_mzid_comparison_file="tests/data/out.csv", stMs_parser_out_file="tests/data/F083270_score.txt")
        c.parse()
        c.compare()
        assert_equal(1, 2)
