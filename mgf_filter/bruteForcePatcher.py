from pyteomics import mgf
import numpy as np


class BruteForcePatcher(object):
    def __init__(self):
        pass

    def patchMgf(self, input_path, output_path):
        maxWindowDiff = 156.10112 + 2 * 1.00782503 + 15.9949146
        with mgf.read(input_path) as spectra:
            spectra_out = []
            for spectrum in spectra:
                int_dic = spectrum['intensity array']
                mz_dic = spectrum['m/z array']
                param_dic = spectrum['params']
                chrg_spec = spectrum['params']['charge'][0]

                pos = 0
                del_array = []
                for m in mz_dic:
                    if m < 175: # smallest y ion - arginin
                        del_array.append(pos)
                    elif m > spectrum['params']['pepmass'][0] * chrg_spec - (chrg_spec - 1) * 1.00782503 - maxWindowDiff:
                        del_array.append(pos)
                    pos += 1

                int_dic = np.delete(int_dic, del_array, 0)
                mz_dic = np.delete(mz_dic, del_array, 0)

                spectra_out.append({'m/z array': mz_dic, 'intensity array': int_dic, 'params': param_dic})

        mgf.write(spectra=spectra_out, output=output_path)
