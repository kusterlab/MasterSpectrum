from mgf_filter.peak import Peak
from pyteomics import mgf
from mgf_filter.util import calculate_Delta_by_ppm
from mgf_filter.masterSpectrum import MasterSpectrum
from mgf_recalibration.util import calculate_tag_tmt10
from mgf_recalibration.util import calculate_ppm_shift
from mgf_recalibration.util import calculate_da_shift


class Recalibrator(object):
    def __init__(self, path, ppm=20, file_out="/dev/null", type='absolute'):
        self.ppm = ppm
        self.path = path
        self.file_out = file_out
        self.type = type

    def load_recalibrate(self):
        fc = calculate_Delta_by_ppm(self.ppm)
        tmt_mass = calculate_tag_tmt10()
        with mgf.read(self.path) as spectra:
            for spectrum in spectra:
                ms = MasterSpectrum()
                params = spectrum['params']
                for mass, intensity in zip(spectrum['m/z array'], spectrum['intensity array']):
                    ms.add(Peak(mass, intensity, fc))

                peak = Peak(tmt_mass, 0.5, fc)
                if peak.key() not in ms.spectrum[0]:
                    recalibrate = False
                else:
                    idx, bin_to_ack, a, b = ms.binary(peak, 0, len(ms.spectrum[0][peak.key()]) - 1, 0)
                    if idx == -1:
                        recalibrate = False
                    else:
                        recalibrate = True
                        recalibration_mass = ms.spectrum[0][peak.key()][idx].mz
                        diff = tmt_mass - recalibration_mass
                        print(params['title'])
                        print("original={0}\tdiff={1}".format(recalibration_mass, diff))

                mass_list = []
                int_list = []
                if recalibrate:
                    ppm_shift = calculate_ppm_shift(diff, tmt_mass)

                for key in ms.spectrum[0].keys():
                    for mp in ms.spectrum[0][key]:
                        if recalibrate:
                            if self.type == 'ppm':
                                diff = calculate_da_shift(mp.mz, ppm_shift)
                                mass_list.append(mp.mz + diff)
                            elif self.type == 'absolute':
                                diff = diff
                                mass_list.append(mp.mz + diff)
                            else:
                                print(self.type)
                                raise ValueError("what did you dooooo")
                        else:
                            mass_list.append(mp.mz)
                        int_list.append(mp.intensity)
                print("len is:\t{0}".format(len(mass_list)))
                mgf.write(spectra=[{'m/z array': mass_list, 'intensity array': int_list, 'params': params}], output=self.file_out)
