
import csv
from mgf_annotation.util import parse_scan_id
from pyteomics import mgf


class SimpleMgfFilter(object):
    def __init__(self, stPathSelection, path_mgf_in, path_mgf_out):
        self.path_mgf_in = path_mgf_in
        self.path_mgf_out = path_mgf_out
        self.stPathSelection = stPathSelection
        self.list_chosen = []

    def read_selection_file(self):
        with open(self.stPathSelection) as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',')
            for row in csvreader:
                if "x" not in row:
                    print(row[0])
                    self.list_chosen.append(int(row[0]))

    def select_mgf(self):
        spectra_out = []
        with mgf.read(self.path_mgf_in) as spectra:
            for spectrum in spectra:

                scanid = int(parse_scan_id(spectrum['params']['title']))
                if scanid in self.list_chosen:
                    spectra_out.append({'m/z array': spectrum['m/z array'],
                                        'intensity array': spectrum['intensity array'],
                                        'params': spectrum['params']})

        mgf.write(spectra=spectra_out, output=self.path_mgf_out)
