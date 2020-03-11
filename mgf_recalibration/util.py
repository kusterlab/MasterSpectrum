
import pyximport
pyximport.install()
from mgf_filter.cython.mat import pow
from pyteomics import mass


def calculate_ppm_shift(da_diff, tmt_mass):
    return da_diff * pow(10, 6) / tmt_mass


def calculate_da_shift(mass, ppm_shift):
    return mass * ppm_shift / pow(10, 6)


def calculate_tag_tmt10():
    base = mass.calculate_mass('C8NO2H20')
    base += 4 * mass.nist_mass['C'][13][0]
    base += mass.nist_mass['N'][15][0]
    base += mass.nist_mass['H+'][1][0]
    return base
