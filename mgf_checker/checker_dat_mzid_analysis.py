import os
import csv
from mgf_checker.util import IonSeriesUsed
from mgf_checker.reader import MZID_Comparison_result
from mgf_checker.reader import MS_parser_result


class Comparer(object):
    """
        reads an out file and an mid file
    """

    def __init__(self, stPath_mzid_comparison_file, stMs_parser_out_file):
        self.stPath_mzid_comparison_file = stPath_mzid_comparison_file
        self.stMs_parser_out_file = stMs_parser_out_file
        self.mzid_info = None
        self.msparser_info = None

    def parse(self):
        self.mzid_info = MZID_Comparison_result(self.stPath_mzid_comparison_file)
        self.mzid_info.read_file()

        self.msparser_info = MS_parser_result(self.stMs_parser_out_file)
        self.msparser_info.read_file()

    def compare(self, output_path="/home/tobiass/Desktop/gaga.csv"):
        with open(output_path, "wt") as csvfile:
            writr = csv.writer(csvfile, lineterminator=os.linesep)
            writr.writerow(("msmsid", "series", "charge", "ratio_not_found", "num_all_peaks", "query", "score", "rank", "numPeaksMatched"))
            msg = "msmsid:\t\t\t{0}\n"
            msg += "series:\t\t\t{1}\n"
            msg += "charges:\t\t\t{2}\n"
            msg += "ratio of not found:\t\t\t{3}\n"
            msg += "num_all_peak:\t\t\t{4}\n"
            msg += "query:\t\t\t{5}\n"
            msg += "score:\t\t\t{6}\n"
            for msmsid in self.msparser_info.data:
                for nRank in self.msparser_info.data[msmsid]:
                    # print(self.msparser_info.data[msmsid].stSeriesUsedStr)
                    ion_info_sig = IonSeriesUsed(self.msparser_info.data[msmsid][nRank].stSeriesUsedStr)
                    for stIon_serie_name in ion_info_sig:
                        for nCharge in ion_info_sig[stIon_serie_name]:
                            if ion_info_sig[stIon_serie_name][nCharge]:

                                query = self.msparser_info.data[msmsid][nRank].nQuery
                                score = self.msparser_info.data[msmsid][nRank].dScore

                                if stIon_serie_name == 'b':
                                    if not self.mzid_info.data[msmsid][nRank].has_Series('B'):
                                        raise ValueError("gn1")
                                    else:
                                        diff = self.mzid_info.data[msmsid][nRank].BSerie.get_ratio(nCharge)
                                        dPeaks_matched = self.mzid_info.data[msmsid][nRank].BSerie.get_num_Peaks_matched_in_serie()
                                        numPeaks = self.mzid_info.data[msmsid][nRank].BSerie.same[nCharge]
                                        # print(msg.format(msmsid, stIon_serie_name, nCharge, diff, numPeaks, query, score))
                                        writr.writerow((msmsid, stIon_serie_name, nCharge, diff, numPeaks, query, score, nRank, dPeaks_matched))
                                elif stIon_serie_name == 'y':
                                    if not self.mzid_info.data[msmsid][nRank].has_Series('Y'):
                                        raise ValueError("gn2")
                                    else:
                                        diff = self.mzid_info.data[msmsid][nRank].YSerie.get_ratio(nCharge)
                                        dPeaks_matched = self.mzid_info.data[msmsid][nRank].YSerie.get_num_Peaks_matched_in_serie()
                                        numPeaks = self.mzid_info.data[msmsid][nRank].YSerie.same[nCharge]
                                        # print(msg.format(msmsid, stIon_serie_name, nCharge, diff, numPeaks, query, score))
                                        writr.writerow((msmsid, stIon_serie_name, nCharge, diff, numPeaks, query, score, nRank, dPeaks_matched))
                                else:
                                    raise ValueError("doooohhh")
