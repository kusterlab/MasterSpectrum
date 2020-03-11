
from mgf_filter.masterSpectrum import MasterSpectrum
from mgf_filter.util import calculate_Delta_by_ppm
from mgf_filter.peak import Peak


def calculate_max_tmt(tmt_set):
    i = -1
    for pos in tmt_set:
        i = max(i, tmt_set[pos])
    return i


def generateMS_by_score_file(score_info_object):
    """
    every object has an mz array
    every object has an diff array
    """
    ms = MasterSpectrum()

    delta_func = calculate_Delta_by_ppm(20)
    dPeaks_matched = 0
    for m, diff in zip(score_info_object['mz'], score_info_object['diff']):
        p = Peak(float(m), 1, delta_func, meta=diff)
        ms.add(p, 0)
        dPeaks_matched += 1

    return ms, dPeaks_matched


def IonSeriesUsed(IonSeries):
    """
    01 - A series
    02 - placeholder (was A - NH3 series in older versions of Mascot)
    03 - A++ series
    04 - B series
    05 - placeholder (was B - NH3 series in older versions of Mascot)
    06 - B++ series
    07 - Y series
    08 - placeholder (was Y - NH3 series in older versions of Mascot)
    09 - Y++ series
    10 - C series
    11 - C++ series
    12 - X series
    13 - X++ series
    14 - Z series
    15 - Z++ series
    16 - Z+H series
    17 - Z+H++ series
    18 - Z+2H series
    19 - Z+2H++ series
    """
    b1_sig = False
    b2_sig = False
    y1_sig = False
    y2_sig = False

    if IonSeries[3] == '2':
        b1_sig = True
    if IonSeries[5] == '2':
        b2_sig = True
    if IonSeries[6] == '2':
        y1_sig = True
    if IonSeries[8] == '2':
        y2_sig = True

    return {'y': {1: y1_sig, 2: y2_sig}, 'b': {1: b1_sig, 2: b2_sig}}
