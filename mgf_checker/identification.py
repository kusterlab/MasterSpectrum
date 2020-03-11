from mgf_checker.ion_series import Ion_series


class Identification(object):
    def __init__(self, mzid_info_lvl_fragmentation, peptide_ref, title):
        self.peptide_ref = peptide_ref
        self.scan_id = title
        self.ion_series_ary = []
        for ion_serie in mzid_info_lvl_fragmentation:
            if "-" not in ion_serie['name']:  # ignoring all NH3 or CO2 series
                self.ion_series_ary.append(Ion_series(ion_serie))

    def report_all_mzs(self):
        """
        returns all mz of all fragements
        """
        mzs = []
        for ser in self.ion_series_ary:
            for mz in ser.mz_ary:
                mzs.append(mz)
        return mzs

    def __str__(self):
        print(self.peptide_ref)
        return {0}.format(self.peptide_ref)
