import os
from pyteomics import pepxml
from pyteomics import mgf


class MgfSeperator(object):
    def __init__(self, pathXml, pathMgf, scoreCutoff):
        self.pathXml = pathXml
        self.pathMgf = pathMgf
        self.file_name = os.path.splitext(os.path.basename(pathMgf))[0]
        self.goodBadUgly = {'good': [], 'bad': []}
        self.goodSpectraCriteria = scoreCutoff

    def getInfoFromPepXml(self):
        count_good = 0
        count_bad = 0
        with pepxml.read(self.pathXml) as psms:
            for psm in psms:
                if 'search_hit' in psm:
                    score = psm['search_hit'][0]['search_score']['ionscore']
                    if score > self.goodSpectraCriteria:
                        count_good += 1
                        self.goodBadUgly['good'].append(psm['spectrum'])
                    else:
                        count_bad += 1
                        self.goodBadUgly['bad'].append(psm['spectrum'])
                else:
                    count_bad += 1
                    self.goodBadUgly['bad'].append(psm['spectrum'])
        print("good: ", count_good)
        print("bad: ", count_bad)

    def seperateMgf(self, folder_out):
        '''
        '''
        for idea in self.goodBadUgly.keys():
            spectra_out = []
            with mgf.read(self.pathMgf) as spectra:
                for spectrum in spectra:
                    int_dic = spectrum['intensity array']
                    mz_dic = spectrum['m/z array']
                    param_dic = spectrum['params']
                    if spectrum['params']['title'] in self.goodBadUgly[idea]:
                        spectra_out.append({'m/z array': mz_dic, 'intensity array': int_dic, 'params': param_dic})
            output_path = folder_out + self.file_name + "_" + str(idea) + ".mgf"
            mgf.write(spectra=spectra_out, output=output_path)
