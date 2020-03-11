
class Ion_series(object):
    def __init__(self, mzid_info):
        self.fragtype = mzid_info['name']
        self.charge = mzid_info['charge']
        self.ions_index = mzid_info['index']
        for peak_info in mzid_info['FragmentArray']:
            if peak_info['measure_ref'] == 'm_mz':
                self.mz_ary = peak_info['values']
            elif peak_info['measure_ref'] == 'm_intensity':
                self.int_ary = peak_info['values']
            elif peak_info['measure_ref'] == 'm_error':
                pass
            else:
                raise ValueError("not known identifier")
