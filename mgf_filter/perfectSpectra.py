import os
import csv
from pyteomics import mass
import itertools


def Nist_mass(element):
    return mass.nist_mass[element][0][0]


def generatePeptides(length):
    allAmino = mass.std_aa_mass.keys()
    for le in range(1, length + 1):
        yield itertools.combinations(allAmino, le)


def calculateDoubleCharged(mz):
    h_atom_mass = Nist_mass('H')
    return (mz + h_atom_mass) / 2.0


def generateMascotIons(length, starter_mass):
    """
    all mascotions used for scoring
    """
    water_mass = 2.0 * Nist_mass('H') + Nist_mass('O')
    amin_mass = 3.0 * Nist_mass('H') + Nist_mass('N')

    for generatorPeptideCombinations in generatePeptides(length):
        for peptides in generatorPeptideCombinations:
            for ion_type in ('b', 'y'):
                ion_type_1 = mass.fast_mass(sequence=peptides,
                                            charge=1,
                                            ion_type=ion_type) + starter_mass
                ion_type_1_star = ion_type_1 - amin_mass
                ion_type_1_o = ion_type_1 - water_mass
                ion_type_2 = calculateDoubleCharged(ion_type_1)
                ion_type_2_star = calculateDoubleCharged(ion_type_1_star)
                ion_type_2_o = calculateDoubleCharged(ion_type_1_o)

                yield([ion_type_1, "".join(peptides) + "_ion_" + str(ion_type) + "_1"],
                      [ion_type_1_star, "".join(peptides) + "_ion_" + str(ion_type) + "_1_star"],
                      [ion_type_1_o, "".join(peptides) + "_ion_" + str(ion_type) + "_1_o"],
                      [ion_type_2, "".join(peptides) + "_ion_" + str(ion_type) + "_2"],
                      [ion_type_2_star, "".join(peptides) + "_ion_" + str(ion_type) + "_2_star"],
                      [ion_type_2_o, "".join(peptides) + "_ion_" + str(ion_type) + "_2_o"])


class GeneratorForPerfectSpectra(object):

    def generateCrazyExclusionList(self, path_dir, length, starter_mass=0.0):
        path = path_dir + "exclusionList" + "_" + str(length) + "_" + str(starter_mass) + ".csv"
        with open(path, "a") as csvfile:
            writr = csv.writer(csvfile, lineterminator=os.linesep)
            writr.writerow(("mz", "comment", "position"))
            for i in generateMascotIons(length, starter_mass):
                for j in i:
                    writr.writerow((j[0], j[1] + "_" + str(starter_mass), "absolute"))

    def addInternalFragments(self, path_dir, length, starter_mass=0.0):
        with open(path_dir + "exclusionList" + "_" + str(length) + "_" + str(starter_mass) + ".csv", "a") as csvfile:
            writr = csv.writer(csvfile, lineterminator=os.linesep)
            writr.writerow(("mz", "comment", "position"))
            for i in generateMascotIons(length, starter_mass):
                for j in i:
                    writr.writerow((j[0], "", "internal"))

    def generateAminoAcidDeltaList(self, path_dir, length, starter_mass=0.0):
        path = path_dir + "exclusionListDelta" + "_" + str(length) + "_" + str(starter_mass) + ".csv"
        with open(path, "a") as csvfile:
            writr = csv.writer(csvfile, lineterminator=os.linesep)
            writr.writerow(("mz", "comment", "position"))
            for i in generatePeptides(length):
                for j in i:
                    mass_pep = mass.fast_mass(j, charge=0, ion_type='b')
                    writr.writerow((mass_pep, "", "absolute"))
