from mgf_filter.masterPeak import MasterPeak
from mgf_filter.masterSpectrum import MasterSpectrum
from mgf_filter.peak import Peak
import csv
from mgf_filter.util import calculate_Delta_by_ppm
import os


class DeltaPatcher(object):
    def __init__(self, delta_func=calculate_Delta_by_ppm(20)):
        self.exclusionSpectrum = MasterSpectrum()
        self.delta_function = delta_func
        self.precursorDeltas = []

    def readExclusionList(self, path):
        '''
        exclusionList: 3 columns
        m/z , comments
        '''
        with open(path, 'r') as csvfile:
            readr = csv.reader(csvfile)
            header = True
            for row in readr:
                if header:
                    header = False
                else:
                    if row[2] == 'absolute':
                        p = Peak(float(row[0]), 1.0, self.delta_function)
                        mp = MasterPeak(p)
                        # 0 for no differentiation of charge states
                        self.exclusionSpectrum.add(mp, 0)
                    elif row[2] == 'precursor':
                        self.precursorDeltas.append(float(row[0]))
                    else:
                        raise ValueError('A very specific bad thing happened')

    def patchDelta(self, input_path, output_path):
        '''
        '''
        with open(output_path, "wt") as csvfile:
            writr = csv.writer(csvfile, lineterminator=os.linesep)
            with open(input_path, 'r') as r_csvfile:
                readr = csv.reader(r_csvfile)
                header = True
                for row in readr:
                    if header:
                        header = False
                        writr.writerow(row)
                    else:
                        peak = Peak(float(row[0]), 0, self.delta_function)
                        if peak.key() in self.exclusionSpectrum.spectrum[0]:
                            idx, bin_to_ack, should_merge_left_peak, should_merge_right_peak = self.exclusionSpectrum.binary(peak, 0, len(self.exclusionSpectrum.spectrum[0][peak.key()]) - 1, 0)
                            if idx != -1:  # found
                                print("found it")
                            else:  # not found
                                writr.writerow(row)
                        else:
                            writr.writerow(row)
