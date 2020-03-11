
from mgf_checker.msParserInfo import MS_parser_info
from abc import ABCMeta, abstractmethod
import csv
from mgf_checker.mzidInfo import MZID_info


class Reader(metaclass=ABCMeta):
    def __init__(self, path):
        self.path = path
        self.data = {}

    def read_file(self):
        with open(self.path) as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',')
            for row in csvreader:
                self.parse(row)

    @abstractmethod
    def parse(self):
        pass


class MZID_Comparison_result(Reader):
    """
    mzid comparison reults reader
    everything is a rank1 result
    expected is calculted by mzid
    found is value of my mgf comparison
    """
    def __init__(self, path):
        super().__init__(path)

    def parse(self, row):
        if "scanid" not in row:
            scanid = int(row[0])
            rank = int(row[1])
            peak = float(row[2])
            position = int(row[3])
            frag = row[4]
            if 'b' in frag:
                frag = 'B'
            elif 'y' in frag:
                frag = 'Y'
            else:
                raise ValueError("strange frag")
            expected = float(row[5])
            found = float(row[6])
            charge = int(row[7])
            dPeaks_matched = int(row[8])

            if scanid in self.data:
                if rank in self.data[scanid]:
                    self.data[scanid][rank].add(stFrag=frag, nPosition=position, dPeak=peak, dExpected=expected, dFound=found, nCharge=charge, dPeaks_matched=dPeaks_matched)
                else:
                    self.data[scanid][rank] = MZID_info(nScanid=scanid, dPeak=peak, nPosition=position, stFrag=frag, dExpected=expected, dFound=found, nCharge=charge, dPeaks_matched=dPeaks_matched)
            else:
                self.data[scanid] = {}
                self.data[scanid][rank] = MZID_info(nScanid=scanid, dPeak=peak, nPosition=position, stFrag=frag, dExpected=expected, dFound=found, nCharge=charge, dPeaks_matched=dPeaks_matched)


class MS_parser_result(Reader):
    """
    store more than first rank result
    """
    def parse(self, row):
        rank_max = 2
        if "query" not in row:

            nQuery = int(row[0])
            dScore = float(row[1])
            if dScore > 0:
                stPeptide = row[2]
                stVarModsStr = row[3]
                stReadableVarMods = row[4]
                nMsmsid = int(row[5])
                stProtein_match = row[6]
                stFilename = row[7]
                nRank = int(row[8])
                stSeriesUsedStr = row[9]

                if nRank <= rank_max:

                    if nMsmsid in self.data:
                        if nRank in self.data[nMsmsid]:
                            pass
                        else:
                            self.data[nMsmsid][nRank] = MS_parser_info(nQuery=nQuery, dScore=dScore, stPeptide=stPeptide, stVarModsStr=stVarModsStr, stReadableVarMods=stReadableVarMods,
                                                                       nMsmsid=nMsmsid,
                                                                       stProtein_match=stProtein_match, stFilename=stFilename, nRank=nRank, stSeriesUsedStr=stSeriesUsedStr)
                    else:
                        self.data[nMsmsid] = {}
                        self.data[nMsmsid][nRank] = MS_parser_info(nQuery=nQuery, dScore=dScore, stPeptide=stPeptide, stVarModsStr=stVarModsStr, stReadableVarMods=stReadableVarMods,
                                                                   nMsmsid=nMsmsid,
                                                                   stProtein_match=stProtein_match, stFilename=stFilename, nRank=nRank, stSeriesUsedStr=stSeriesUsedStr)
