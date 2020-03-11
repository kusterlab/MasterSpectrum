from mgf_filter.util import calculatePrecursor
from mgf_filter.util import calculate_Delta_by_ppm
from pyteomics import mgf
from mgf_filter.util import calculateRelativeIntensity
from mgf_filter.peak import Peak
from mgf_filter.masterSpectrum import MasterSpectrum


class PrecursorDeltaPeaksExtractor(object):
    def __init__(self, path):
        self.path = path
        self.masterSpectrum = MasterSpectrum()

    def createDeltaPrecursorMasterSpectrum(self, delta_func=calculate_Delta_by_ppm(20)):
        with mgf.read(self.path) as spectra:
            for spectrum in spectra:
                int_dic = spectrum['intensity array']
                mz_dic = spectrum['m/z array']
                chrg_spec = spectrum['params']['charge'][0]
                precursor = calculatePrecursor(mz=spectrum['params']['pepmass'][0], charge=chrg_spec)

                rel_int = calculateRelativeIntensity(int_dic)
                for m, i in zip(mz_dic, rel_int):
                    p = Peak(precursor - float(m), float(i), delta_func)
                    self.masterSpectrum.add(p, 0)

    def exportCsv(self, output_path):
        self.masterSpectrum.export_to_csv(output_path)
