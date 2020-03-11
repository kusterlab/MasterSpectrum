from pyteomics import pepxml
import csv
import os


class Extractor(object):
    def __init__(self, path):
        self.result = []
        self.path = path

    def digest_pepxml(self):
        with pepxml.read(self.path) as psms:
            for psm in psms:
                psm_result = []
                if 'search_hit' in psm:
                    psm_result.append(psm['spectrum'])
                    psm_result.append(psm['search_hit'][0]['peptide'])
                    psm_result.append(psm['search_hit'][0]['search_score']['ionscore'])
                    psm_result.append(len(psm['search_hit'][0]['proteins']))
                    psm_result.append(psm['search_hit'][0]['num_matched_ions'])
                else:
                    psm_result.append(psm['spectrum'])
                    psm_result.append('')
                    psm_result.append(0)
                    psm_result.append(0)
                    psm_result.append(0)
                self.result.append(psm_result)

    def write_csv_output(self, output_path):
        with open(output_path, "wt") as csvfile:
            writr = csv.writer(csvfile, lineterminator=os.linesep)
            writr.writerow(("spectrum", "peptide", "ionscore_rank1", "num_tot_proteins", "num_matched_ions"))
            for psm in self.result:
                writr.writerow((psm[0], psm[1], psm[2], psm[3], psm[4]))
