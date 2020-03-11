import os
import csv
from mgf_annotation.util import calculate_allowed_Mass_Diff
from mgf_filter.util import calculate_Delta_by_ppm
from mgf_filter.peak import Peak
from pyteomics import mgf
from mgf_annotation.util import parse_scan_id
from pyteomics import mzid
from mgf_checker.peptide_evidence import PeptideEvidence
from mgf_checker.identification import Identification
from mgf_checker.reference import Reference
from mgf_checker.util import calculate_max_tmt
from mgf_checker.util import generateMS_by_score_file


class Checker(object):
    def __init__(self, path):
        self.path = path
        self.identifications = []
        self.mgf_reads = {}
        self.peptide_evidence = {}

    def read_enhanced_spectrum(self, path):
        """
        saving a masterSpectrum for every spectrum is not a memory efficient idea
        4 GB for 8247 spectra
        instead creating on request (saving reference and spectra object)

        """
        with mgf.read(path) as spectra:
            for spectrum in spectra:
                # charge_of_spectrum = str(spectrum['params']['charge'][0])
                scan_id = parse_scan_id(spectrum['params']['title'])
                self.mgf_reads[scan_id] = Reference(scan_id, spectrum['m/z array'])

    def read_score_file(self, path_score, allowed_ids):
        """
        input:
        path
        allowed_ids

        reads improved csv
        returns hashset: scan_id : list of mz and diffs

        info:
        originated_mz is read because after decharging mz is shifted and can not be used for my comparison
        """
        read_spectra = {}
        with open(path_score) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if int(row['id1']) in allowed_ids:
                    id = int(row['id1'])
                    if id not in read_spectra:
                        read_spectra[id] = {'mz': [float(row['originated mz'])], 'diff': [float(row['annotation'])]}
                    else:
                        read_spectra[id]['mz'].append(float(row['originated mz']))
                        read_spectra[id]['diff'].append(float(row['annotation']))

        return read_spectra

    def parse_mz_id(self):
        """
        reading mzid
        saving every spectrum identification (but just rank 1)
        returns:
        None
        """
        data = mzid.read(self.path)

        for d in data:
            title = parse_scan_id(d['spectrum title'])
            ident = {}
            len_ranks = len(d['SpectrumIdentificationItem'])
            if len_ranks > 1:
                for i in [0, 1]:
                    identification = d['SpectrumIdentificationItem'][i]  # 0 because just first rank
                    peptide_ref = identification['peptide_ref']
                    ident[i + 1] = Identification(mzid_info_lvl_fragmentation=identification['IonType'],
                                                  peptide_ref=peptide_ref,
                                                  title=title)
            else:
                identification = d['SpectrumIdentificationItem'][0]  # 0 because just first rank
                peptide_ref = identification['peptide_ref']
                ident[1] = Identification(mzid_info_lvl_fragmentation=identification['IonType'],
                                          peptide_ref=peptide_ref,
                                          title=title)

            self.identifications.append(ident)

    def parse_mz_id_peptide_ref(self):
        """
        reads mzid and reports peptide evidence
        following a read of peptide information and combining these two
        returns:
        None
        """
        data = mzid.iterfind(self.path, "PeptideEvidence")
        for d in data:
            self.peptide_evidence[d['peptide_ref']] = PeptideEvidence(isDecoy=d['isDecoy'],
                                                                      start=d['start'],
                                                                      end=d['end'],
                                                                      peptide_ref=d['peptide_ref'],
                                                                      id=d['id'],
                                                                      pre=d['pre'],
                                                                      dBSequence_ref=d['dBSequence_ref'],
                                                                      post=d['post'])
        data = mzid.iterfind(self.path, "Peptide")
        for d in data:
            if d['id'] not in self.peptide_evidence:
                raise ValueError("{0} is not in peptide_evidence".format(d['id']))
            else:
                self.peptide_evidence[d['id']].set_peptide_modification(d)

    def analyze_mzid_id_vs_score_file(self, path_score, output_path="/home/tobiass/Desktop/out.csv"):
        """
        Prerequirenment:
        loading mzid(spectrum and peptide info)

        input:
        path to score file (created by control.py improve_csv)

        Description:
        - finding all important scan ids
        - loading all scans as spectra from score file
            - saving delta as meta information for every peak
        -
        cheching all identifiacation (--> peptide_ref)
            -   for peptide ref:
                    read series (b/y, index and charge, mz)
                    generate tag amount for every position
                    tag series of b/y with tmt diff

                    for series:
                        check in scans for delta

        Output:
        none
        """
        with open(output_path, "wt") as csvfile:
            writr = csv.writer(csvfile, lineterminator=os.linesep)
            writr.writerow(("scanid", "rank", "peak", "position", "frag", "expected by mascot", "found", "charge", "PeaksMatchedInSerie", "max_number_of_tmt_tag"))
            valid_scan_ids = {}
            delta_func = calculate_Delta_by_ppm(20)

            for identification_hashobject in self.identifications:
                valid_scan_ids[identification_hashobject[1].scan_id] = True

            valid_scan_ids = valid_scan_ids.keys()  # quick and dirty - smarter with some kind of unique tree
            valid_scan_ids = [int(i) for i in valid_scan_ids]
            valid_spectra_info = self.read_score_file(path_score, valid_scan_ids)

            self.parse_mz_id_peptide_ref()

            for identification_hashobject in self.identifications:
                for i in identification_hashobject:
                    peptide_ref = identification_hashobject[i].peptide_ref
                    scan_id = int(identification_hashobject[i].scan_id)

                    b_tmt, y_tmt = self.peptide_evidence[peptide_ref].get_annotated_positions()
                    max_tmt_tag = calculate_max_tmt(b_tmt)  # same num for b and y
                    ms, dPeaks_matched = generateMS_by_score_file(valid_spectra_info[scan_id])
                    for ion_serie in identification_hashobject[i].ion_series_ary:
                        tmt_masses = calculate_allowed_Mass_Diff(n=max_tmt_tag, z=ion_serie.charge)
                        if 'b' in ion_serie.fragtype:
                            tmt_pos_ary = b_tmt
                        elif 'y' in ion_serie.fragtype:
                            tmt_pos_ary = y_tmt
                        else:
                            raise ValueError("{0}\tis not valid fragtype".format(ion_serie.fragtype))

                        for pos, mz in zip(ion_serie.ions_index, ion_serie.mz_ary):
                            num_tag_at_pos = tmt_pos_ary[pos]
                            if num_tag_at_pos == 0:
                                expected_mass_delta = 0
                            else:
                                expected_mass_delta = tmt_masses[num_tag_at_pos][ion_serie.charge]
                            peak = Peak(mz, 1, delta_func)
                            if peak.key() in ms.spectrum[0]:
                                idx, bin_to_ack, should_merge_left_peak, should_merge_right_peak = ms.binary(peak, 0, len(ms.spectrum[0][peak.key()]) - 1, 0)
                                if idx == -1:
                                    # just heavy spectra are just all peaks of the heavy labeled partner of a mix (there is no garantee that both
                                    # have the same number of peaks
                                    # BUT the number of spectra is the same
                                    msg = "Peak:\t{0}\n"
                                    msg += "frag:\t{4}\n"
                                    msg += "position:\t{5}\n"
                                    msg += "scanid:\t{1}\n"
                                    msg += "expected:\t{2}\n"
                                    msg += "found:\t{3}\n"
                                    msg += "peak was just part of heavy spectra(but same number before \".\"\n"
                                    # print(msg.format(mz, scan_id, expected_mass_delta, -2,
                                    #                 ion_serie.fragtype, pos))
                                    writr.writerow((scan_id, i, mz, pos, ion_serie.fragtype, expected_mass_delta, -2, ion_serie.charge, dPeaks_matched, max_tmt_tag))
                                    # raise ValueError("identified Peak:\t{0}\nscanid:\t{1}\ncouldnt be found in spectra".format(mz, scan_id, expected_mass_delta, -1,
                                    #                                                                                           ion_serie.fragtype, pos))
                                else:
                                    found_mass_diff = ms.spectrum[0][peak.key()][idx].meta
                                    msg = "Peak:\t{0}\n"
                                    msg += "frag:\t{4}\n"
                                    msg += "position:\t{5}\n"
                                    msg += "scanid:\t{1}\n"
                                    msg += "expected:\t{2}\n"
                                    msg += "found:\t{3}\n"
                                    if found_mass_diff != expected_mass_delta:
                                        writr.writerow((scan_id, i, mz, pos, ion_serie.fragtype, expected_mass_delta, found_mass_diff, ion_serie.charge, dPeaks_matched, max_tmt_tag))
                                        # raise ValueError(msg.format(mz, scan_id, expected_mass_delta, found_mass_diff, ion_serie.fragtype, pos))
                                    else:
                                        # here is my goal
                                        writr.writerow((scan_id, i, mz, pos, ion_serie.fragtype, expected_mass_delta, found_mass_diff, ion_serie.charge, dPeaks_matched, max_tmt_tag))
                            else:
                                # just heavy spectra are just all peaks of the heavy labeled partner of a mix (there is no garantee that both
                                # have the same number of peaks
                                # BUT the number of spectra is the same
                                msg = "Peak:\t{0}\n"
                                msg += "frag:\t{4}\n"
                                msg += "position:\t{5}\n"
                                msg += "scanid:\t{1}\n"
                                msg += "expected:\t{2}\n"
                                msg += "found:\t{3}\n"
                                msg += "peak was just part of heavy spectra\n"
                                # print(msg.format(mz, scan_id, expected_mass_delta, -1,
                                #                 ion_serie.fragtype, pos))
                                writr.writerow((scan_id, i, mz, pos, ion_serie.fragtype, expected_mass_delta, -1, ion_serie.charge, dPeaks_matched, max_tmt_tag))
                                # raise ValueError("identified Peak:\t{0}\nscanid:\t{1}\ncouldnt be found in spectra".format(mz, scan_id, expected_mass_delta, -1,
                                #                                                                                           ion_serie.fragtype, pos))

    def analyse_mzid_vs_mgf(self):
        delta_func = calculate_Delta_by_ppm(20)
        if self.mgf_reads == {}:
            raise ValueError("need the mgf beforehand read_enhanced_spectrum")
        else:
            for ids in self.identifications:
                mzs = ids.report_all_mzs()
                ms = self.mgf_reads[ids.scan_id].request_ms()

                for mz in mzs:

                    peak = Peak(mz, 2, delta_func)
                    if peak.key() in ms.spectrum[0]:
                        idx, bin_to_ack, should_merge_left_peak, should_merge_right_peak = ms.binary(peak, 0, len(ms.spectrum[0][peak.key()]) - 1, 0)
                        if idx == -1:
                            error = "mz:\t{0}\nscan_id:\t{1}".format(mz, ids.scan_id)
                            raise ValueError(error)
                        else:
                            pass
                    else:
                        error = "mz:\t{0}\nscan_id:\t{1}".format(mz, ids.scan_id)
                        raise ValueError(error)
