
from mgf_filter.masterSpectrum import MasterSpectrum
from mgf_filter.util import calculate_Delta_by_ppm
from mgf_filter.peak import Peak


class Reference(object):
    def __init__(self, scan_id, mz_ary):
        self.scan_id = scan_id
        self.mz_ary = mz_ary

    def request_ms(self):
        ms = MasterSpectrum()

        delta_func = calculate_Delta_by_ppm(20)
        for m in self.mz_ary:
            p = Peak(float(m), 1, delta_func)
            ms.add(p, 0)

        return ms
