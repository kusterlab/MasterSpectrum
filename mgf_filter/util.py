import math
import datetime
import pyximport; pyximport.install()
from mgf_filter.cython.mat import pow


def calculate_Resolution_based_MZ(mz):
    """
    based on MZ calculates resolution
    provided by a LR
    """
    return math.pow(10, 5.847 + math.log10(mz) * (-0.546))


def calculate_Delta_based_MZ(mz, numSigma=1):
    """
    based on mz a resolution is calculated and multiplied by sec. (not necessary)
    Multiplier is 1 by default
    """
    return numSigma * mz / calculate_Resolution_based_MZ(mz)


def calculate_Delta_Fixed(delta):
    """
    returns a function that can use
    """
    def fixed(mz):
        return delta
    return fixed


def calculateRelativeIntensity(aIntensity):
    maxIntensity = max(aIntensity)
    return [x / maxIntensity for x in aIntensity]


def calculate_Delta_by_ppm(ppm):
    def fix_ppm(mz):
        return ppm * float(mz) / (pow(10, 6))
    return fix_ppm


def timeStamped(fname, fmt='%Y-%m-%d-%H-%M-%S_{fname}'):
    return datetime.datetime.now().strftime(fmt).format(fname=fname)


def calculatePrecursor(mz, charge):
    return mz * charge - (charge - 1) * 1.00782503
