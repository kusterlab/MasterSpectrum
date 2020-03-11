from mgf_filter.util import calculateRelativeIntensity
import os
import csv
from mgf_annotation.util import gen_allowed_mass_diff_with_sign
from mgf_annotation.util import calculate_Mass_Diff
from mgf_annotation.util import parse_scan_id
from mgf_filter.peak import Peak
from mgf_filter.masterSpectrum import MasterSpectrum
from pyteomics import mgf
from mgf_filter.util import calculate_Delta_by_ppm
from blist import sortedlist
from mgf_filter.util import calculatePrecursor
from mgf_annotation.util import mass_diff_decharging_stuff


class Reference(object):
    """
    Every Reference:
    Spectrum1 <---> Spectrum2
    Therefore
    peak list1 = all peaks of 1 spectrum
    Therefore a complete MGF consists of N references
    """
    def __init__(self, ppm, id_2, id_1, peak_list_2, peak_list_1, mass_2, mass_1, charge, extra_mass, int_list_2, int_list_1, params2, params1):
        self.ppm = ppm
        self.charge = int(charge)
        self.extra_mass = extra_mass
        self.fc = calculate_Delta_by_ppm(20)

        if float(mass_2) > float(mass_1):
            self.mass_1 = mass_2
            self.mass_2 = mass_1

            self.id_1 = int(id_2)
            self.id_2 = int(id_1)

            self.peak_list_1 = peak_list_2
            self.peak_list_2 = peak_list_1

            self.int_list_1 = int_list_2
            self.int_list_2 = int_list_1

            self.params = params2

        else:
            self.mass_1 = mass_1
            self.mass_2 = mass_2

            self.id_1 = int(id_1)
            self.id_2 = int(id_2)

            self.peak_list_1 = peak_list_1
            self.peak_list_2 = peak_list_2

            self.int_list_1 = int_list_1
            self.int_list_2 = int_list_2

            self.params = params1

            self.params["#numtags"] = "{0}".format("-1")

    def __str__(self):
        return "a={0} \t b={1}\n{2}\n{3}".format(self.id_1, self.id_2, self.peak_list_1, self.peak_list_2)

    def just_create_heavy_ms(self):
        ms = MasterSpectrum()
        for mass, i in zip(self.peak_list_1, self.int_list_1):
            ms.add(Peak(mass, i, self.fc, meta={'delta': 511, 'mass': -1, 'decharged': False}))
        return ms

    def create_ms(self, iMin_similarity):
        """
        int list 1 and mz list 1 are the ones with a higher precursor mass
        int/mz list1 will create a master spectrum

        for every peak in list2
        check all allowed mass diffs (depending on charge and number of tags)
        -> find peak + delta
        --> for all save a similarity value (calculating ratio of int1/int2) depending which is bigger

        -> sort for highest similarity score
        -> check if masterspectrum peak already has an refering peak2
        --> if new similarity is higher: replace

        a peak must have at least int_similarity = 0.5

        """
        req_min_similarity = iMin_similarity
        ms = MasterSpectrum()
        rel_int = calculateRelativeIntensity(self.int_list_1)
        for mass, rel_int, i in zip(self.peak_list_1, rel_int, self.int_list_1):
            ms.add(Peak(mass, i, self.fc, meta={'delta': 511, 'mass': -1, 'decharged': False, 'similarity': -1, 'originated_mz': mass, 'rel_int': rel_int}))

        num = round((abs(self.extra_mass) / calculate_Mass_Diff()))  # num of tags possible
        self.params["#numtags"] = "{0}".format(num)
        deltas = mass_diff_decharging_stuff(n=int(num), z=self.charge)
        dd = deltas.keys()
        dd = list(dd)
        dd.sort()
        rel_int = calculateRelativeIntensity(self.int_list_2)
        for mass, rel_int, i in zip(self.peak_list_2, rel_int, self.int_list_2):
            similarity_most_similar = -1
            idx_most_similar = -1
            peak_key_most_similar = -1
            delta_most_similar = -1
            mass_most_similar = -1

            for delta in dd:
                peak = Peak(mass + delta, i, self.fc)
                if peak.key() in ms.spectrum[0]:
                    idx, bin_to_ack, should_merge_left_peak, should_merge_right_peak = ms.binary(peak, 0, len(ms.spectrum[0][peak.key()]) - 1, 0)
                    if idx != -1:  # found
                        rel_int1 = ms.spectrum[0][peak.key()][idx].meta['rel_int']
                        ratio = rel_int1 / rel_int
                        if ratio > 1:  # i want ratio to be between 0 - 1
                            ratio = 1 / ratio

                        if ratio > similarity_most_similar:
                            similarity_most_similar = ratio
                            idx_most_similar = idx
                            peak_key_most_similar = peak.key()
                            mass_most_similar = mass
                            delta_most_similar = delta

            if similarity_most_similar > req_min_similarity:
                if similarity_most_similar > ms.spectrum[0][peak_key_most_similar][idx_most_similar].meta['similarity']:
                    ms.spectrum[0][peak_key_most_similar][idx_most_similar].meta = {'delta': delta_most_similar,
                                                                                    'mass': mass_most_similar,
                                                                                    'decharged': deltas[delta_most_similar]['decharge']['state'],
                                                                                    'similarity': similarity_most_similar,
                                                                                    'originated_mz': ms.spectrum[0][peak_key_most_similar][idx_most_similar].mz,
                                                                                    'rel_int': ratio}

                    if deltas[delta_most_similar]['decharge']['state']:
                        ms.add(Peak(calculatePrecursor(ms.spectrum[0][peak_key_most_similar][idx_most_similar].mz, deltas[delta_most_similar]['decharge']['z']),
                                    ms.spectrum[0][peak_key_most_similar][idx_most_similar].intensity,
                                    self.fc,
                                    meta=ms.spectrum[0][peak_key_most_similar][idx_most_similar].meta))
                        del(ms.spectrum[0][peak_key_most_similar][idx_most_similar])
                        if len(ms.spectrum[0][peak_key_most_similar]) == 0:
                            del(ms.spectrum[0][peak_key_most_similar])

        return ms


class MgfAnnotater(object):
    def __init__(self, path, output_path, min_rel_similarity, ppm=10):
        self.ppm = ppm
        self.ms = MasterSpectrum()
        self.path = path
        self.out = []
        self.output_path = output_path
        self.references = sortedlist(key=lambda i: i.id_1)
        self.min_rel_similarity = min_rel_similarity

    def load_msconvert_mgf(self):
        """
        creates references based on precursor mass
        a missing scanid means an ms1 event
        by default referencing works just within one ms2 block
        """
        fc = calculate_Delta_by_ppm(self.ppm)
        scan_id_ary = []
        problems = []
        error = 0
        with mgf.read(self.path) as spectra:
            for spectrum in spectra:
                mass = spectrum['params']['pepmass'][0]
                precursor_chrg = int(spectrum['params']['charge'][0])
                mass = calculatePrecursor(mass, precursor_chrg)

                scanid = int(parse_scan_id(spectrum['params']['title']))
                if len(scan_id_ary) == 0:
                    scan_id_ary.append(scanid)
                else:
                    if scanid != scan_id_ary[-1] + 1:
                        if len(scan_id_ary) % 2 == 1:
                            problems.append(scan_id_ary[0])
                            error += 1
                            scan_id_ary = []
                            scan_id_ary.append(scanid)
                        else:
                            scan_id_ary = []
                            scan_id_ary.append(scanid)
                        self.ms = MasterSpectrum()  # new MS if scan_id group (seperated by ms1) is completed
                    else:
                        scan_id_ary.append(scanid)

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
                    self.ms.add(Peak(mass, scanid, fc, meta={'ms': spectrum['m/z array'], 'int': spectrum['intensity array'], 'params': spectrum['params']}),
                                charge=precursor_chrg)
            if error > 0:
                print(" delete valid information {0}".format(error))
                # raise ValueError("taddaaa")

    def export_annotated_spectra_to_csv(self):
        """
        id1
        id2
        mz
        intensityy
        annotation:     delta mass
        mz2:
        decharged:      True/False
        similarity:     intensity similaritz of mz1 vs mz2
        originated mz:  if decharged the original mz is saved
        """
        print("start exporting information")
        with open(self.output_path, "wt") as csvfile:
            writr = csv.writer(csvfile, lineterminator=os.linesep)
            writr.writerow(("id1", "id2", "mz", "intensity", "annotation", "mz2", "decharged", 'similarity', 'originated mz', "extra_mass"))
            for ref in self.references:
                ms = ref.create_ms(iMin_similarity=self.min_rel_similarity)
                for chrg in ms.spectrum:
                    for key in ms.spectrum[chrg].keys():
                        for mp in ms.spectrum[chrg][key]:
                            if mp.meta['mass'] != -1:
                                writr.writerow((ref.id_1, ref.id_2, mp.mz, mp.intensity, mp.meta['delta'], mp.meta['mass'], mp.meta['decharged'], mp.meta['similarity'], mp.meta['originated_mz'], ref.extra_mass))

    def export_annotated_spectra_to_mgf(self, mgf_path, report_just_heavy=False):
        spectra_out = []
        for ref in self.references:
            if report_just_heavy:
                ms = ref.just_create_heavy_ms()
            else:
                ms = ref.create_ms(iMin_similarity=self.min_rel_similarity)
            buf_peaks = []
            buf_int = []
            for chrg in ms.spectrum:
                for key in ms.spectrum[chrg].keys():
                    for mp in ms.spectrum[chrg][key]:
                        if report_just_heavy:
                            buf_peaks.append(mp.mz)
                            buf_int.append(mp.intensity)
                        else:
                            if mp.meta['mass'] != -1:
                                buf_peaks.append(mp.mz)
                                buf_int.append(mp.intensity)
                if len(buf_peaks) != 0:
                    spectra_out.append({'m/z array': buf_peaks, 'intensity array': buf_int, 'params': ref.params})
        mgf.write(spectra=spectra_out, output=mgf_path)
