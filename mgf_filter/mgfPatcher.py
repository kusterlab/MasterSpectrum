from pyteomics import mgf
import numpy as np
from mgf_filter.masterPeak import MasterPeak
from mgf_filter.masterSpectrum import MasterSpectrum
from mgf_filter.peak import Peak
import csv
from mgf_filter.util import calculate_Delta_by_ppm
from mgf_filter.util import calculatePrecursor


class MgfPatcher(object):
    def __init__(self, delta_func=calculate_Delta_by_ppm(20)):
        self.exclusionSpectrum = MasterSpectrum()
        self.delta_function = delta_func
        self.precursorDeltas = []

    def readExclusionList(self, path):
        '''
        exclusionList: 2 columns
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

    def patchMgf(self, input_path, output_path):
        '''
        '''

        with mgf.read(input_path) as spectra:
            spectra_out = []
            for spectrum in spectra:
                int_dic = spectrum['intensity array']
                mz_dic = spectrum['m/z array']
                param_dic = spectrum['params']

                chrg_spec = spectrum['params']['charge'][0]
                precursor = calculatePrecursor(mz=spectrum['params']['pepmass'][0], charge=chrg_spec)
                pos = 0
                del_array = []
                for m in mz_dic:
                    peak = Peak(m, 0, self.delta_function)
                    if peak.key() in self.exclusionSpectrum.spectrum[0]:
                        idx, bin_to_ack, should_merge_left_peak, should_merge_right_peak = self.exclusionSpectrum.binary(peak, 0, len(self.exclusionSpectrum.spectrum[0][peak.key()]) - 1, 0)
                        if idx != -1:  # found
                            del_array.append(pos)
                    else:
                        mp = MasterPeak(peak)
                        for precursorDelta in self.precursorDeltas:
                            if mp.isInsideMz(precursor - precursorDelta):
                                del_array.append(pos)
                            else:
                                pass
                    pos += 1

                int_dic = np.delete(int_dic, del_array, 0)
                mz_dic = np.delete(mz_dic, del_array, 0)

                spectra_out.append({'m/z array': mz_dic, 'intensity array': int_dic, 'params': param_dic})

        mgf.write(spectra=spectra_out, output=output_path)
