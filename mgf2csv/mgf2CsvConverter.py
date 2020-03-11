
from mgf_filter.util import calculateRelativeIntensity
from pyteomics import mgf
import csv
from mgf_annotation.util import parse_scan_id
import os


class Mgf2CsvConverter(object):
    def __init__(self, mgf_path, output_csv):
        self.mgf_path = mgf_path
        self.output_csv = output_csv

    def load_export(self):
        with mgf.read(self.mgf_path) as spectra, open(self.output_csv, "wt") as csvfile:
            writr = csv.writer(csvfile, lineterminator=os.linesep)
            writr.writerow(("scanid", "peak", "rel_int"))
            for spectrum in spectra:
                scan_id = parse_scan_id(spectrum['params']['title'])
                st = "scanid:\t{0}\n"
                print(st.format(scan_id))
                rel_int = calculateRelativeIntensity(spectrum['intensity array'])
                for m, i in zip(spectrum['m/z array'], rel_int):
                    writr.writerow((scan_id,
                                    m,
                                    i))
