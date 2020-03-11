
from mgf_checker.dataContainer import DataContainer
from mgf_checker.ionSerie import IonSerie


class MZID_info(DataContainer):
    def __init__(self, nScanid, dPeak, nPosition, stFrag, dExpected, dFound, nCharge, dPeaks_matched):
        self.nScanid = nScanid
        self.BSerie = IonSerie()
        self.YSerie = IonSerie()
        self.hasB = False
        self.hasY = False

        self.add(stFrag=stFrag, nPosition=nPosition, dPeak=dPeak, dExpected=dExpected, dFound=dFound, nCharge=nCharge, dPeaks_matched=dPeaks_matched)

    def add(self, stFrag, nPosition, dPeak, dExpected, dFound, nCharge, dPeaks_matched):
        if stFrag == 'Y':
            self.YSerie.add(nPosition, dPeak, dExpected, dFound, nCharge, dPeaks_matched)
            self.hasY = True
        elif stFrag == 'B':
            self.BSerie.add(nPosition, dPeak, dExpected, dFound, nCharge, dPeaks_matched)
            self.hasB = True
        else:
            raise ValueError("i cannot land here")

    def has_Series(self, stSerie):
        if stSerie == 'B':
            return self.hasB
        elif stSerie == 'Y':
            return self.hasY
        else:
            raise ValueError("nowhere land")

    def ratio_explained(self):
        if self.hasB:
            self.BSerie
        pass
