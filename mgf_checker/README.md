# Checker
This program should read mzid files and create a list of all identified fragments
Mapping starts from SpectrumIdentificationResult --> Fragmentation + peptide_ref
peptide_ref shows peptide sequence in normal annotation
Fragmentions are in order for b ions (1, 2, 3) and in reverse order for y ions (5, 3, 2)

Every fragment ion (with and without losses - meaning NH2 and CO2) should be found in a dual
(after enhancing) spectrum.

Validation of this idea:
P015940 - B12 - R3 
comparing results of just_heavy (given by msconvert) [mascot id : 6244]
vs
improved spectra (given by msconvert)  [mascot id : 6286]

Aim:
in the end perhaps unvalidating rank1 results
