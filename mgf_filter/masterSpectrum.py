import os
from pyteomics import mgf
from mgf_filter.util import calculateRelativeIntensity
from mgf_filter.util import calculate_Delta_by_ppm
from blist import sortedlist
from mgf_filter.masterPeak import MasterPeak
from mgf_filter.peak import Peak
import csv
import pyximport; pyximport.install()
from mgf_filter.cython.mat import ceil


class MasterSpectrum:
    """
    Master spectrum takes peaks and converts them to master peaks
    if there are no master peaks to be inserted
    spectrum is {charge of spectrum: {bins: MP}}
    """
    def __init__(self, logger=None):
        self.spectrum = {}
        self.merged = 0
        self.appended = 0
        self.multimerged = 0  # within one insertion 3 peaks are merged

    def binary(self, peak, imin, imax, charge):
        """
        Input values:
        peak
        imin is minimum search position
        imax is max search poistion
        charge defines in which masterspectrum to search
        Return values:
        first argument: position of insert, -1 if can not be added to any peak
        sc argument : add peak from left or right bin
        rd argument : must left peak also be added (merge case)
        4th argument: right peak also be added (think about 3 peaks and a merge between 2 and 3)
        """
        key = peak.key()
        imid = int(ceil((imax + imin) / 2))
        if imax < imin:
            exist_left_bin = key - 1 in self.spectrum[charge].keys()
            exist_right_bin = key + 1 in self.spectrum[charge].keys()

            if exist_left_bin:
                if self.spectrum[charge][key - 1][-1].isInside(peak):
                    return -1, -1, False, False

            if exist_right_bin:
                if self.spectrum[charge][key + 1][0].isInside(peak):
                    return -1, 1, False, False

            return -1, 0, False, False

        # search must go on
        mspeak_in_bin = self.spectrum[charge][key][imid]
        if mspeak_in_bin.greater(peak):
            return self.binary(peak, imin, imid - 1, charge)
        elif mspeak_in_bin.smaller(peak):
            return self.binary(peak, imid + 1, imax, charge)
        # search results in peak that should be added
        else:
            if imid == 0:
                exist_left_bin = key - 1 in self.spectrum[charge]
                exist_right_bin = key + 1 in self.spectrum[charge]
                if exist_left_bin:
                    if self.spectrum[charge][key - 1][-1].isInside(peak):
                        return 0, -1, False, False
                    else:
                        return 0, 0, False, False  # peak must be added to peak in pos1 but but left bin can be ignored
                elif exist_right_bin:
                    if self.spectrum[charge][key + 1][0].isInside(peak):
                        return 0, 1, False, False
                    else:
                        return 0, 0, False, False  # peak must be added to peak in pos1 but but left bin can be ignored
                else:
                    return 0, 0, False, False

            else:  # peak is somewhere between 1 and last
                # imid -1 exists alway
                is_last_entry = len(self.spectrum[charge][key]) - 1 == imid
                exists_bigger_peak = not(is_last_entry)
                if exists_bigger_peak:
                    if self.spectrum[charge][key][imid - 1].isInside(peak):
                        return imid, 0, True, False
                    if self.spectrum[charge][key][imid + 1].isInside(peak):
                        return imid, 0, False, True
                    else:
                        return imid, 0, False, False
                else:  # is last entry
                    exist_right_bin = key + 1 in self.spectrum[charge]
                    if self.spectrum[charge][key][imid - 1].isInside(peak):
                        return imid, 0, True, False
                    elif exist_right_bin:
                        if self.spectrum[charge][key + 1][0].isInside(peak):
                            return imid, 1, False, False
                        else:
                            return imid, 0, False, False
                    else:
                        return imid, 0, False, False

    def add(self, peak, charge=0):
        """
        charge is by default 0, so if a spectrum should be summed without including charge information
        it is defaulted to 0
        """
        key = peak.key()
        if charge not in self.spectrum:
            self.spectrum[charge] = {}

        if key in self.spectrum[charge]:
            idx, bin_to_ack, should_merge_left_peak, should_merge_right_peak = self.binary(peak, 0, len(self.spectrum[charge][key]) - 1, charge)
            if idx == -1:  # does not have to react to merge cases!
                if bin_to_ack == 0:
                    self.appended += 1
                    if type(peak).__name__ == 'Peak':
                        self.spectrum[charge][key].add(MasterPeak(peak))
                    else:
                        self.spectrum[charge][key].add(peak)
                elif bin_to_ack == -1:
                    if len(self.spectrum[charge][key]) == 0:
                        del(self.spectrum[charge][key])
                    get_masterPeak_left_bin = self.spectrum[charge][key - 1][-1]
                    del(self.spectrum[charge][key - 1][-1])
                    if len(self.spectrum[charge][key - 1]) == 0:
                        del(self.spectrum[charge][key - 1])
                    get_masterPeak_left_bin.add(peak)
                    self.merged += 1
                    if (get_masterPeak_left_bin.key() not in self.spectrum[charge]):
                        self.spectrum[charge][get_masterPeak_left_bin.key()] = sortedlist(key=lambda i: i.left)
                    self.spectrum[charge][get_masterPeak_left_bin.key()].add(get_masterPeak_left_bin)
                else:  # +1
                    if len(self.spectrum[charge][key]) == 0:
                        del(self.spectrum[charge][key])
                    get_masterPeak_right_bin = self.spectrum[charge][key + 1][0]
                    del(self.spectrum[charge][key + 1][0])
                    if len(self.spectrum[charge][key + 1]) == 0:
                        del(self.spectrum[charge][key + 1])
                    get_masterPeak_right_bin.add(peak)
                    self.merged += 1
                    if (get_masterPeak_right_bin.key() not in self.spectrum[charge]):
                        self.spectrum[charge][get_masterPeak_right_bin.key()] = sortedlist(key=lambda i: i.left)
                    self.spectrum[charge][get_masterPeak_right_bin.key()].add(get_masterPeak_right_bin)

            else:  # idx != -1
                if bin_to_ack == 0:
                    if should_merge_left_peak:
                        get_left_masterPeak = self.spectrum[charge][key][idx - 1]
                        get_idx_masterPeak = self.spectrum[charge][key][idx]
                        del(self.spectrum[charge][key][idx - 1:idx + 1])  # delete is slicing, therefore one more
                        if len(self.spectrum[charge][key]) == 0:
                            del(self.spectrum[charge][key])
                        get_idx_masterPeak.add(get_left_masterPeak)
                        get_idx_masterPeak.add(peak)
                        self.multimerged += 1
                        if (get_idx_masterPeak.key() not in self.spectrum[charge]):
                            self.spectrum[charge][get_idx_masterPeak.key()] = sortedlist(key=lambda i: i.left)
                        self.spectrum[charge][get_idx_masterPeak.key()].add(get_idx_masterPeak)
                    elif should_merge_right_peak:
                        get_idx_masterPeak = self.spectrum[charge][key][idx]
                        get_right_masterPeak = self.spectrum[charge][key][idx + 1]
                        del(self.spectrum[charge][key][idx:idx + 2])  # delete is slicing, therefore one more
                        if len(self.spectrum[charge][key]) == 0:
                            del(self.spectrum[charge][key])
                        get_idx_masterPeak.add(get_right_masterPeak)
                        get_idx_masterPeak.add(peak)
                        self.multimerged += 1
                        if (get_idx_masterPeak.key() not in self.spectrum[charge]):
                            self.spectrum[charge][get_idx_masterPeak.key()] = sortedlist(key=lambda i: i.left)
                        self.spectrum[charge][get_idx_masterPeak.key()].add(get_idx_masterPeak)
                    else:  # case idx = 0, 0, False, False
                        get_idx_masterPeak = self.spectrum[charge][key][idx]
                        del(self.spectrum[charge][key][idx])  # delete is slicing, therefore one more
                        if len(self.spectrum[charge][key]) == 0:
                            del(self.spectrum[charge][key])
                        get_idx_masterPeak.add(peak)
                        self.merged += 1
                        if (get_idx_masterPeak.key() not in self.spectrum[charge]):
                            self.spectrum[charge][get_idx_masterPeak.key()] = sortedlist(key=lambda i: i.left)
                        self.spectrum[charge][get_idx_masterPeak.key()].add(get_idx_masterPeak)
                elif bin_to_ack == -1:  # idx != -1
                    # 0, -1, F, F
                    get_idx_masterPeak = self.spectrum[charge][key][idx]
                    get_masterPeak_before = self.spectrum[charge][key - 1][-1]
                    del(self.spectrum[charge][key][idx])
                    if(len(self.spectrum[charge][key]) == 0):
                        del(self.spectrum[charge][key])
                    del(self.spectrum[charge][key - 1][-1])
                    if(len(self.spectrum[charge][key - 1]) == 0):
                        del(self.spectrum[charge][key - 1])
                    get_idx_masterPeak.add(get_masterPeak_before)
                    get_idx_masterPeak.add(peak)
                    self.merged += 1
                    if (get_idx_masterPeak.key() not in self.spectrum[charge]):
                        self.spectrum[charge][get_idx_masterPeak.key()] = sortedlist(key=lambda i: i.left)
                    self.spectrum[charge][get_idx_masterPeak.key()].add(get_idx_masterPeak)
                else:  # bin_to_ack == 1, idx !=  -1
                    get_idx_masterPeak = self.spectrum[charge][key][idx]
                    get_masterPeak_after = self.spectrum[charge][key + 1][0]
                    del(self.spectrum[charge][key][idx])
                    if(len(self.spectrum[charge][key]) == 0):
                        del(self.spectrum[charge][key])
                    del(self.spectrum[charge][key + 1][0])
                    if(len(self.spectrum[charge][key + 1]) == 0):
                        del(self.spectrum[charge][key + 1])
                    get_idx_masterPeak.add(get_masterPeak_after)
                    get_idx_masterPeak.add(peak)
                    self.merged += 1
                    if (get_idx_masterPeak.key() not in self.spectrum[charge]):
                        self.spectrum[charge][get_idx_masterPeak.key()] = sortedlist(key=lambda i: i.left)
                    self.spectrum[charge][get_idx_masterPeak.key()].add(get_idx_masterPeak)
        else:
            self.spectrum[charge][key] = sortedlist(key=lambda i: i.left)
            self.add(peak, charge)

    def printStats(self):
        size = 0
        for key in self.spectrum.keys():
            size += len(self.spectrum[key])
        print("m:", self.merged, "a:", self.appended, "s:", size, "av:",
              size / len(self.spectrum),
              "cc:", ((1.0 + self.merged + self.appended) / size))

    def printSize(self):
        for charges in self.spectrum.keys():
            for key in self.spectrum.keys():
                print(charges, key, len(self.spectrum[key]))

    def __str__(self):
        info = ""
        info += "appended: " + str(self.appended) + "\n merged: " + str(self.merged) + "\n"
        for charges in self.spectrum.keys():
            info += "Charge :" + str(charges) + "\n" + "#######################" + "\n"
            for key in self.spectrum[charges].keys():
                info += "Key : " + str(key) + "\n" + "===================" + "\n"
                for mp in self.spectrum[charges][key]:
                        info += str(mp) + "\n" + "-----------------" + "\n"
        return info + "###########################"

    def load_from_mgf(self, path, ignoreCharges, delta_func=calculate_Delta_by_ppm(20)):
        up = 0
        with mgf.read(path) as spectra:
            for spectrum in spectra:
                up = up + 1
                rel_int = calculateRelativeIntensity(spectrum['intensity array'])
                charge_of_spectrum = str(spectrum['params']['charge'][0])
                for m, i in zip(spectrum['m/z array'], rel_int):
                    p = Peak(float(m), float(i), delta_func)
                    if ignoreCharges:
                        self.add(p, 0)
                    else:
                        self.add(p, charge_of_spectrum)

    def export_to_csv(self, path):
        with open(path, "wt") as csvfile:
            writr = csv.writer(csvfile, lineterminator=os.linesep)
            writr.writerow(("mz", "intensity", "counts", "left border", "right border", "start_mz", "ms1_charge", "rel_intensity_ratio", "counts_ratio"))
            for charges in self.spectrum.keys():
                for key in self.spectrum[charges].keys():
                    for mp in self.spectrum[charges][key]:
                        writr.writerow((mp.mz, mp.intensity, mp.counts, mp.left, mp.right, mp.mz_origin, charges, mp.rel_intensity_ratio, mp.counts_ratio))

    def load_from_csv(self, path, delta_function=calculate_Delta_by_ppm(20)):
        with open(path, 'r') as csvfile:
            readr = csv.reader(csvfile)
            header = True
            for row in readr:
                if header:
                    header = False
                else:
                    p = Peak(float(row[0]), float(row[1]), delta_function)
                    mp = MasterPeak(p)
                    mp.counts = int(row[2])
                    mp.mz_origin = float(row[5])
                    self.add(mp, str(row[6]))

    def __eq__(self, other):
        """
        reports true if both have the same member variables!
        normally you test for the same memory address
        """
        if self.multimerged != other.multimerged:
            print("there1")
            return False
        if self.appended != other.appended:
            print("there2")
            return False
        if self.merged != other.merged:
            print("there3")
            return False
        for charges in self.spectrum.keys():
            for key in self.spectrum[charges].keys():
                i = 0
                for mp in self.spectrum[charges][key]:
                    if mp != other.spectrum[charges][key][i]:
                        print(mp)
                        print(other.spectrum[charges][key][i])
                        print("there")
                        return False
                    i += 1
        return True

    def compare_other_ms(self, ms2):
        """
        compare master spectrum to a second masterspectrum
        returns a new masterSpectrum
        """
        # load self in a new MasterSpectrum
        ms3 = MasterSpectrum()
        for charges in self.spectrum.keys():
            for key in self.spectrum[charges].keys():
                for mp in self.spectrum[charges][key]:
                    mp.rel_intensity_ratio = 0.5
                    mp.counts_ratio = 0.5
                    ms3.add(mp, charges)

        # check for all peaks in ms2, if nothing comparable is in self --> add to ms3
        often = 0
        for charges in ms2.spectrum.keys():
            for key in ms2.spectrum[charges].keys():
                for mp in ms2.spectrum[charges][key]:
                    if charges not in ms3.spectrum.keys():
                        mp.rel_intensity_ratio = -0.5
                        mp.counts_ratio = -0.5
                        ms3.add(mp, charges)
                    elif mp.key() not in ms3.spectrum[charges]:
                        mp.rel_intensity_ratio = -0.5
                        mp.counts_ratio = -0.5
                        ms3.add(mp, charges)
                    else:
                        idx, bin_to_ack, should_merge_left_peak, should_merge_right_peak = ms3.binary(mp, 0, len(ms3.spectrum[charges][mp.key()]) - 1, charges)
                        if idx == -1:
                            if bin_to_ack == 0:
                                mp.rel_intensity_ratio = -0.5
                                mp.counts_ratio = -0.5
                                ms3.add(mp, charges)
                            elif bin_to_ack == -1:
                                ms3.spectrum[charges][mp.key() - 1][-1].recalculate_ratio(mp)
                            else:  # +1
                                ms3.spectrum[charges][mp.key() + 1][0].recalculate_ratio(mp)
                        else:  # idx != -1
                            if bin_to_ack == 0:
                                if should_merge_left_peak:
                                    get_masterPeak_before = ms3.spectrum[charges][mp.key()][idx - 1]
                                    get_masterPeak = ms3.spectrum[charges][mp.key()][idx]
                                    get_masterPeak.add(get_masterPeak_before)
                                    del(ms3.spectrum[charges][mp.key()][idx - 1: idx + 1])
                                    if(len(ms3.spectrum[charges][mp.key()]) == 0):
                                        del(ms3.spectrum[charges][mp.key()])
                                    get_masterPeak.recalculate_ratio(mp)
                                    ms3.add(get_masterPeak, charges)
                                    often += 1
                                elif should_merge_right_peak:
                                    get_masterPeak_after = ms3.spectrum[charges][mp.key()][idx + 1]
                                    get_masterPeak = ms3.spectrum[charges][mp.key()][idx]
                                    get_masterPeak.add(get_masterPeak_after)
                                    del(ms3.spectrum[charges][mp.key()][idx: idx + 2])
                                    if(len(ms3.spectrum[charges][mp.key()]) == 0):
                                        del(ms3.spectrum[charges][mp.key()])
                                    get_masterPeak.recalculate_ratio(mp)
                                    ms3.add(get_masterPeak, charges)
                                    often += 1
                                else:
                                    ms3.spectrum[charges][mp.key()][idx].recalculate_ratio(mp)
                            elif bin_to_ack == -1:
                                get_masterPeak_binBefore = ms3.spectrum[charges][mp.key() - 1][-1]
                                get_masterPeak = ms3.spectrum[charges][mp.key()][idx]
                                get_masterPeak.add(get_masterPeak_binBefore)
                                del(ms3.spectrum[charges][mp.key() - 1][-1])
                                if(len(ms3.spectrum[charges][mp.key() - 1]) == 0):
                                    del(ms3.spectrum[charges][mp.key() - 1])

                                del(ms3.spectrum[charges][mp.key()][idx])
                                if(len(ms3.spectrum[charges][mp.key()]) == 0):
                                    del(ms3.spectrum[charges][mp.key()])

                                get_masterPeak.recalculate_ratio(mp)
                                ms3.add(get_masterPeak, charges)

                                often += 1
                            else:
                                get_masterPeak_binAfter = ms3.spectrum[charges][mp.key() + 1][0]
                                get_masterPeak = ms3.spectrum[charges][mp.key()][idx]
                                get_masterPeak.add(get_masterPeak_binAfter)
                                del(ms3.spectrum[charges][mp.key() + 1][0])
                                if(len(ms3.spectrum[charges][mp.key() + 1]) == 0):
                                    del(ms3.spectrum[charges][mp.key() + 1])

                                del(ms3.spectrum[charges][mp.key()][idx])
                                if(len(ms3.spectrum[charges][mp.key()]) == 0):
                                    del(ms3.spectrum[charges][mp.key()])

                                get_masterPeak.recalculate_ratio(mp)
                                ms3.add(get_masterPeak, charges)
                                often += 1
        print("how often: ", often)
        return ms3
