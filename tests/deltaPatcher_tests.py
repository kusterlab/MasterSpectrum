
from nose.tools import *
import os
import shutil
from mgf_filter.deltaPatcher import DeltaPatcher
from mgf_filter.util import calculate_Delta_by_ppm
from mgf_filter.perfectSpectra import GeneratorForPerfectSpectra
from mgf_filter.deltaExtractor import DeltaExtractor
import csv

class TestDeltaPatcher(object):

    @classmethod
    def teardown_class(cls):
        shutil.rmtree("tests/data/temp")

    @classmethod
    def setup_class(cls):
        if not os.path.exists("tests/data/temp"):
            os.makedirs("tests/data/temp")

    def test_init(self):
        deltaP = DeltaPatcher(delta_func=calculate_Delta_by_ppm(20))
        gps = GeneratorForPerfectSpectra()
        gps.generateAminoAcidDeltaList("tests/data/temp/", 1, 0)
        deltaP.readExclusionList("tests/data/temp/exclusionListDelta_1_0.csv")
        dEx = DeltaExtractor("tests/data/cetsa_101287_A01_P013190_S00_N01_R1_TMT10_3spectra_shorted.mgf")
        dEx.createDeltaMasterSpectrum(0)
        dEx.exportCsv("tests/data/temp/a.csv")
        out = "tests/data/temp/b.csv"
        out = "/home/tobiass/b.csv"
        deltaP.patchDelta(input_path="tests/data/temp/a.csv", output_path=out)
        with open(out, 'r') as csvfile:
            readr = csv.reader(csvfile)
            header = True
            for row in readr:
                if header:
                    header = False
                else:
                    if '71.' in row[0]:
                        raise ValueError('A very specific bad thing happened')
