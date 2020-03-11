
from mgf_filter.util import calculate_Delta_by_ppm
from pyteomics import mgf
from mgf_filter.util import calculateRelativeIntensity
from mgf_filter.peak import Peak
from mgf_filter.masterSpectrum import MasterSpectrum


class DeltaExtractor(object):
    def __init__(self, path):
        self.path = path
        self.masterSpectrum = MasterSpectrum()

    def createDeltaMasterSpectrum(self, min_rel_intensity, delta_func=calculate_Delta_by_ppm(20)):
        up = 0
        with mgf.read(self.path) as spectra:
            for spectrum in spectra:
                print(up)
                up += 1
                int_dic = spectrum['intensity array']
                mz_dic = spectrum['m/z array']

                rel_int = calculateRelativeIntensity(int_dic)
                smallerArea = [(i, j) for i, j in zip(mz_dic, rel_int) if j >= min_rel_intensity]
                mz_dic = [i for i, j in smallerArea]
                rel_int = [j for i, j in smallerArea]
                for i in range(len(mz_dic) - 1, -1, -1):
                    for j in range(i - 1, - 1, -1):
                        diff = mz_dic[i] - mz_dic[j]
                        p = Peak(diff, rel_int[j], delta_func)  # intensities are from lower peak
                        self.masterSpectrum.add(p, 0)

    def exportCsv(self, output_path):
        self.masterSpectrum.export_to_csv(output_path)


if __name__ == "__main__":
    t = [3, 5, 10, 15, 20]
    for i in range(len(t) - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            print(t[j] - t[i])
