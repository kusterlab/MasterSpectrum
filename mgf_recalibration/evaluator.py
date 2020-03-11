import csv
from statistics import median
import os
from pyteomics import mzid
from mgf_annotation.util import parse_scan_id


class Fragment(object):
    def __init__(self, name, indice, charge, mz, error, scanid):
        self.name = name
        self.indice = indice
        self.charge = charge
        self.mz = mz
        self.error = error
        self.scanid = scanid


class Scans(object):
    def __init__(self):
        self.data = {}


class Evaluator(object):
    def __init__(self, path):
        self.path = path
        self.scan1 = Scans()

    def parse_mz_id(self):
        data = mzid.read(self.path)
        max_rank = 2
        for d in data:
            title = d['spectrum title']
            scan_id = parse_scan_id(title)
            fragments = []
            len_frag = len(d['SpectrumIdentificationItem'])
            pos = 1
            while(pos <= min(max_rank, len_frag)):
                for fragmentation in d['SpectrumIdentificationItem'][pos - 1]['IonType']:  # 0 because just first rank
                    for f in fragmentation['FragmentArray']:
                        if f['measure_ref'] == 'm_mz':
                            mz = f['values']
                        elif f['measure_ref'] == 'm_error':
                            error = f['values']
                        else:
                            pass
                    fragments.append(Fragment(name=fragmentation['name'],
                                     indice=fragmentation['index'],
                                     charge=fragmentation['charge'],
                                     mz=mz,
                                     error=error,
                                     scanid=parse_scan_id(title)))
                if scan_id in self.scan1.data:
                    self.scan1.data[scan_id][pos] = fragments
                else:
                    self.scan1.data[scan_id] = {}
                    self.scan1.data[scan_id][pos] = fragments
                pos += 1

    def report_average_error(self, behaviour):
        n = 0
        data = 0
        for scan_id in self.scan1.data:
            for fragment in self.scan1.data[scan_id]:
                for err in fragment.error:
                    if behaviour == 'all':
                        n += 1
                        data += abs(err)
                    elif behaviour == 'ignoreNH3':
                        if "NH3" in fragment.name:
                            pass
                        else:
                            n += 1
                            data += abs(err)
                    else:
                        raise ValueError("Problem")

        print("error (av.)\t{0}".format(data / n))

    def export_error(self, stOutputPath, bIgnore_non_standard=True):
        """
        bIgnore_non_standard controls if error for NH3 and CO2 loss is reported
        """
        ppm = pow(10, 6)
        with open(stOutputPath, "wt") as csvfile:
            writr = csv.writer(csvfile, lineterminator=os.linesep)
            writr.writerow(("msmsid", "rank", "series", "charge", "error_avg", "ppm_err_avg", "max_ppm_err", "median_ppm_err"))
            for scan_id in self.scan1.data:
                for nRank in self.scan1.data[scan_id]:
                    for fragment in self.scan1.data[scan_id][nRank]:
                        aError = [abs(i) for i in fragment.error]
                        ppmError = [abs(error) * ppm / mz for error, mz in zip(fragment.error, fragment.mz)]
                        err = float(sum(aError)) / len(aError) if len(aError) > 0 else float('nan')
                        ppm_err = float(sum(ppmError)) / len(ppmError) if len(ppmError) > 0 else float('nan')
                        max_ppm_err = max(ppmError)
                        median_ppm_err = median(ppmError)
                        if bIgnore_non_standard:
                            if "-" in fragment.name:
                                pass
                            else:
                                name = 'b' if 'b' in fragment.name else 'y'
                                writr.writerow((scan_id,
                                                nRank,
                                                name,
                                                fragment.charge,
                                                err,
                                                ppm_err,
                                                max_ppm_err,
                                                median_ppm_err))
                        else:
                            writr.writerow((scan_id,
                                            nRank,
                                            fragment.name,
                                            fragment.charge,
                                            err,
                                            ppm_err,
                                            max_ppm_err,
                                            median_ppm_err))
