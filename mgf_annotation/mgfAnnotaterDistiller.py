from mgf_annotation.util import gen_allowed_mass_diff_with_sign
from mgf_filter.peak import Peak
from mgf_filter.masterSpectrum import MasterSpectrum
from pyteomics import mgf
from mgf_filter.util import calculate_Delta_by_ppm
from mgf_annotation.util import parse_scan_id
import csv
from mgf_filter.util import calculatePrecursor
from mgf_annotation.mgfAnnotater import MgfAnnotater
from mgf_annotation.mgfAnnotater import Reference


class MgfAnnotaterDistiller(MgfAnnotater):
    """
    def __init__(self, path, output_path, ppm=10):
        self.ppm = ppm
        self.ms = MasterSpectrum()
        self.path = path
        self.out = []
        self.output_path = output_path
        self.references = sortedlist(key=lambda i: i.id_1)
    """

    def load_improved_csv(self, path_score):
        spectra_to_be_referenced = {}
        with open(path_score) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if row['id1'] > row['id2']:
                    spectra_to_be_referenced[int(row['id1'])] = int(row['id2'])
                else:
                    spectra_to_be_referenced[int(row['id2'])] = int(row['id1'])

        self.ids_to_be_referenced = spectra_to_be_referenced

    def load_distiller_mgf(self):
        """
        creates references based on improved csv
        a missing scanid means an ms1 event
        """
        data = {}

        alm = [i for i in gen_allowed_mass_diff_with_sign(n=4, z=1)]

        with mgf.read(self.path) as spectra:
            for spectrum in spectra:
                mass = spectrum['params']['pepmass'][0]
                precursor_chrg = int(spectrum['params']['charge'][0])
                mass = calculatePrecursor(mass, precursor_chrg)

                scanid = int(parse_scan_id(spectrum['params']['title']))

                if scanid in self.ids_to_be_referenced:
                    if self.ids_to_be_referenced[scanid] in data:
                        mass1 = data[self.ids_to_be_referenced[scanid]]['params']['pepmass'][0]
                        precursor_chrg1 = int(data[self.ids_to_be_referenced[scanid]]['params']['charge'][0])
                        mass1 = calculatePrecursor(mass1, precursor_chrg1)
                        diff = abs(mass1 - mass)
                        diff2 = [abs(diff - abs(i)) for i in alm]
                        pos = diff2.index(min(diff2))
                        p = "mass1:\t {0}\n"
                        p += "mass:\t {1}\n"
                        p += "scanid:\t {2}\n"
                        p += "charge:\t {3}\n"
                        p += "charge2:\t {4}\n"
                        p += "scanid2:\t {5}\n"
                        if diff > 21:  # distiller changes precursor charge therefore precurosr mass calculation is wrong
                            print(p.format(mass1, mass, scanid, precursor_chrg1, spectrum['params']['charge'][0], self.ids_to_be_referenced[scanid]))
                            print(diff)
                            print(diff2)
                            print("----------------")
                        else:
                            self.references.add(Reference(ppm=self.ppm,
                                                          id_2=scanid,
                                                          id_1=self.ids_to_be_referenced[scanid],  # also scanid
                                                          peak_list_2=spectrum['m/z array'],
                                                          peak_list_1=data[self.ids_to_be_referenced[scanid]]['m/z array'],
                                                          mass_2=mass,
                                                          mass_1=mass1,
                                                          charge=spectrum['params']['charge'][0],
                                                          extra_mass=alm[pos],
                                                          int_list_2=spectrum['intensity array'],
                                                          int_list_1=data[self.ids_to_be_referenced[scanid]]['intensity array'],
                                                          params2=spectrum['params'],
                                                          params1=data[self.ids_to_be_referenced[scanid]]['params']))
                        del(data[self.ids_to_be_referenced[scanid]])
                        del(self.ids_to_be_referenced[scanid])
                else:
                    data[scanid] = spectrum

    def load_distiller_mgf2(self):
        """
        creates references based on precursor mass
        a missing scanid means an ms1 event
        by default referencing works just within one ms2 block
        """
        fc = calculate_Delta_by_ppm(self.ppm)
        error = 0
        self.ms = MasterSpectrum()
        with mgf.read(self.path) as spectra:
            for spectrum in spectra:
                mass = spectrum['params']['pepmass'][0]
                precursor_chrg = int(spectrum['params']['charge'][0])
                mass = calculatePrecursor(mass, precursor_chrg)

                scanid = int(parse_scan_id(spectrum['params']['title']))
                found = False
                if len(self.ms.spectrum) == 0:
                    peak = Peak(mass, scanid, fc)
                    self.ms.add(Peak(mass, scanid, fc, meta={'ms': spectrum['m/z array'], 'int': spectrum['intensity array'], 'params': spectrum['params']}),
                                charge=precursor_chrg)
                    found = True
                else:
                    if (precursor_chrg in self.ms.spectrum.keys()):  # react to charge !!!!!!
                        if len(self.ms.spectrum[precursor_chrg]) == 0:
                            peak = Peak(mass, scanid, fc)
                            self.ms.add(Peak(mass, scanid, fc, meta={'ms': spectrum['m/z array'], 'int': spectrum['intensity array'], 'params': spectrum['params']}),
                                        charge=precursor_chrg)
                            found = True
                        else:
                            for extra_mass in gen_allowed_mass_diff_with_sign(n=4, z=1):
                                if found is False:
                                    peak = Peak(mass + extra_mass, 0.5, fc)
                                    if peak.key() in self.ms.spectrum[precursor_chrg]:
                                        print(precursor_chrg)
                                        idx, bin_to_ack, a, b = self.ms.binary(peak, 0, len(self.ms.spectrum[precursor_chrg][peak.key()]) - 1, precursor_chrg)
                                        if idx != -1:
                                            self.references.add(Reference(ppm=self.ppm,
                                                                          id_2=scanid,
                                                                          id_1=self.ms.spectrum[precursor_chrg][peak.key()][idx].intensity,  # also scanid
                                                                          peak_list_2=spectrum['m/z array'],
                                                                          peak_list_1=self.ms.spectrum[precursor_chrg][peak.key()][idx].meta['ms'],
                                                                          mass_2=mass,
                                                                          mass_1=self.ms.spectrum[precursor_chrg][peak.key()][idx].mz,
                                                                          charge=spectrum['params']['charge'][0],
                                                                          extra_mass=extra_mass,
                                                                          int_list_2=spectrum['intensity array'],
                                                                          int_list_1=self.ms.spectrum[precursor_chrg][peak.key()][idx].meta['int'],
                                                                          params2=spectrum['params'],
                                                                          params1=self.ms.spectrum[precursor_chrg][peak.key()][idx].meta['params']))
                                            found = True
                                            del(self.ms.spectrum[precursor_chrg][peak.key()][idx])
                                            if len(self.ms.spectrum[precursor_chrg][peak.key()]) == 0:
                                                del(self.ms.spectrum[precursor_chrg][peak.key()])
                                                if len(self.ms.spectrum[precursor_chrg]) == 0:
                                                    del(self.ms.spectrum[precursor_chrg])

                if found is False:
                    limit_scan_id = scanid - 20  # could start at -19
                    ms_bac = MasterSpectrum()
                    for chrg in self.ms.spectrum:
                        for key in self.ms.spectrum[chrg].keys():
                            for mp in self.ms.spectrum[chrg][key]:
                                if mp.intensity >= limit_scan_id:
                                    ms_bac.add(mp, charge=chrg)
                    self.ms = ms_bac
                    self.ms.add(Peak(mass, scanid, fc, meta={'ms': spectrum['m/z array'], 'int': spectrum['intensity array'], 'params': spectrum['params']}),
                                charge=precursor_chrg)
            if error > 0:
                print(" delete valid information {0}".format(error))
                # raise ValueError("taddaaa")
