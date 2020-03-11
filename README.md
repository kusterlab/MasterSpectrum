
# MasterSpectrum
This software was used for my master thesis in 2016. The aim was to find diagnostic ions and eventually remove them to improve identification of MS/MS spectra. It finally shined in 2018 as a tool to quickly find systematic side products of 21 post-translational protein modifications https://doi.org/10.1074/mcp.TIR118.000783.

The original filtering option was also the reason for repository name *mgf_filter*

## What is does?
**MasterSpectrum** is a representation of one or more spectra and consists of MasterPeaks.
The aggregation has to be fuzzy in m/z space because same measurenment of the same peptide can result in different decimal places. After inserting a spectrum into a Masterspectrum every peak is normalized to a relative value between [0, 1] based on the most
intense peak. This is necessary to provide an easier way to compare intensity between different measurenments.


Insertion of a peak triggers a search for similar MasterPeaks. A peak is defined as similar
to a MasterPeak if the peak lies within a predefined m/z window of a found MasterPeak.
These windows can be defined relative to the m/z value of the actual peak (ppm) or by
a fix distance (Da). As a standard setting for all analysis a windows size of 20 ppm was
used.

### A search can result in three different actions:
1. If no similar peak is found, a peak is converted to MasterPeak and is placed into MasterSpectrum according to its m/z value. Every new MasterPeak saves also the number of already aggregated peaks and its starting m/z value
2. If a similar peak is found, the search has to be further refined.  
  - If the next MasterPeak above or below the found m/z value is not similar to the searched peak, the first found MasterPeak and the inserted peak merge 
  - If the next MasterPeak above or below the found m/z is also similar to the searched peak, all three peaks are merged Merging of one or two MasterPeaks with a new peak results alwawys in a new MasterPeak. The M/z values are adjusted to the weighted average of all found MasterPeaks and the new peak. 


To speed up insertion, all function calls that consist of numeric calculation (i) (re)calculation of m/z windows, (ii) weighted averages and (iii) normalisation were replaced with C functions. 

Results of MasterSpectrum can be saved as a CSV file for analysis or reimported if more aggregation with additional MGF files is necessary.
Because MasterSpectrum allows a fast search in m/z space it was also used for different tasks e.g. creating a filter for peaks.
To ensure correctness of MasterSpectrum within its specification and because it works on loosly defined files as MGF files a test suit of 40 test cases for the basic implementation of MasterSpectrum was added. All test cases can be run with minimised real data and
takes about 20 seconds. After all a code coverage of about 82 % was achieved.

# Installation
The package is not (yet) pushed to Pypi, but it is easily installable via pip (requires a up-to-date pip version `pip install -U pip`)
```
pip install -e git+https://github.com/kusterlab/MasterSpectrum.git@v1.0#egg=mgf_filter
```

## Platform Dependency
MasterSpectrum is available on Windows, UNIX/Linux and Mac OSX and the only necessary dependency is a working Python installation and a C compiler.

# How to run it

```
MasterSpectrum create /MGF/PATH/lorem.mgf /tmp/ True
```
First argument => path to mgf
Second argument => output folder
Third argument => aggregation take **not** care of precursor charge

## Output

| mz               | intensity          | counts | left border      | right border     | start_mz  | ms1_charge | rel_intensity_ratio | counts_ratio | 
|------------------|--------------------|--------|------------------|------------------|-----------|------------|---------------------|--------------| 
| 175.10254        | 0.0107771314851379 | 1      | 175.0990379492   | 175.1060420508   | 175.10254 | 0          | 0                   | 0            | 
| 175.106629101904 | 1.60120116043193   | 27     | 175.103126969322 | 175.110131234486 | 175.10474 | 0          | 0                   | 0            | 
| 175.1121535397   | 0.119175527922278  | 3      | 175.108651296629 | 175.115655782771 | 175.11229 | 0          | 0                   | 0            | 
| 175.119032460438 | 430.916364238763   | 4338   | 175.115530079789 | 175.122534841087 | 175.11354 | 0          | 0                   | 0            | 
| 175.133620190728 | 0.25654267792363   | 9      | 175.130117518324 | 175.137122863132 | 175.13344 | 0          | 0                   | 0            | 

### Explanation

| Output variable     	| explanation                                                                         	|
|---------------------	|-------------------------------------------------------------------------------------	|
| mz                  	| mz of the aggregated MasterPeak                                                     	|
| intensity           	| intensity weighted m/z average of all peaks                                         	|
|                     	| number of aggregated peaks in a MasterPeak                                          	|
| left border         	| left border of MasterPeak                                                           	|
| right border        	| right border of MasterPeak                                                          	|
| start_mz            	| original mz value of the first peak in a MasterPeak (to investigate drift of peaks) 	|
| ms1_charge          	| deprecated                                                                          	|
| rel_intensity_ratio 	| feature used for folded MasterSpectra                                               	|
| counts_ratio        	| feature used for folded MasterSpectra                                               	|


## Help
**MasterSpectrum** provides a lot of different functionalities and the help manual should provide a good starting point.

```
>>> MasterSpectrum -h 


usage: MasterSpectrum [-h]
                      {create,fold,patch,bruteforce_patch,extract,precursor,mascotExclusion,seperateMgfbyScore,allDelta,reinstate,deleteAminoAcidDelta,apl}
                      ...

mgf_filter

positional arguments:
  {create,fold,patch,bruteforce_patch,extract,precursor,mascotExclusion,seperateMgfbyScore,allDelta,reinstate,deleteAminoAcidDelta,apl}
                        commands
    create              takes an MGF files and creates a MasterSpectrum based
                        on a 20 ppm window, it can ignore charges
    fold                it takes two MGF files and creates a relative spectrum
                        between both runs based on a 20 ppm window. it cant
                        ignore charges
    patch               creates a MasterSpectrum based on 1-N exclusion lists
    bruteforce_patch    creates a MasterSpectrum similar to the already
                        published results
    extract             extract important information from a pepxml
    precursor           generates a MasterSpectrum with m/z values relative to
                        precursor value
    mascotExclusion     delete peaks based on a provided peptide length and/or
                        an additional mass
    seperateMgfbyScore  generate two spectra seperated by a provided score
                        cuttoff
    allDelta            calculate all possible deltas given a minimal relative
                        intensity
    reinstate           takes two mgfs and reports difference in intensity
                        world as a new MGF file
    deleteAminoAcidDelta
                        takes a list of possible amino acid deltas and deletes
                        such peaks in a MasterSpectrum
    apl                 convert an folder full of apl files to a single MGF
                        file

optional arguments:
  -h, --help            show this help message and exit

```
