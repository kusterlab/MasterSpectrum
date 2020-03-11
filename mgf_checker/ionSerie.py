
from mgf_checker.dataContainer import DataContainer


class IonSerie(DataContainer):
    def __init__(self):
        self.byCharge = {}
        self.diff = {}
        self.same = {}
        self.dPeaks_matched = -1  # number of peaks matched for scan id

    def add(self, nPosition, dPeak, dExpected, dFound, nCharge, dPeaks_matched):
        self.dPeaks_matched = dPeaks_matched
        if nCharge not in self.byCharge:
            self.byCharge[nCharge] = {}
            self.diff[nCharge] = 0
            self.same[nCharge] = 0

        self.byCharge[nCharge][nPosition] = {'Peak': dPeak, 'dExpected': dExpected, 'dFound': dFound}
        if dExpected != dFound:
            self.diff[nCharge] += 1
        self.same[nCharge] += 1

    def get_by_charge(self, nCharge):
        if nCharge not in self.byCharge:
            return {}, False
        else:
            return self.byCharge[nCharge], True

    def get_ratio(self, nCharge):
        """
        returns how many of explained peaks (by mascot) can be matched to peak with same diff (to lighter spectrum) as expected
        e.g. 8 peaks were matched and all of them were part of 21 matched peaks --> ratio 1

        """
        return self.diff[nCharge] / self.same[nCharge]

    def get_num_Peaks_matched_in_serie(self):
        """
        returns number of peaks matched for this scanid
        not at all depending on series, charge, ...
        """
        return self.dPeaks_matched
