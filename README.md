# Prey Size Analysis

This depository contains the code & datasets used in **Santin *et. al*,  Modulation of prey size reveals adaptability and robustness in the cell cycle of an intracellular predator., 2022**. In this work, time lapse microscopy was used to monitor the *Bdellovibrio* cell cycle at the single cell level. The depository can be divided in two parts: the *analysis pipeline* used to obtain the data, and the *experiments*, where the datasets from different experiments are gathered. The experiments folder also contains the scripts used for the visualization of each dataset.

## Pipeline

The pipeline is written as an R markdown file `Pipeline_final.Rmd`, where the user can choose to run the complete script or run the analysis chunk-by-chunk. The pipeline loads the miscroscopy & segmentation data (see below), detects the popping time, and (depending on the data type) detects & counts foci and/or performs growth analysis based on cytoplasmic fluorescence.

### Data curation: before running the pipeline

The pipeline contains code to process movies of bdellovibrio containing fluorescent foci and/or cytoplasmic signal & uses phase-contrast to detect the bdelloplasts. As input it needs folders containing

* A .tiff or .tiffs of the full phase contrast movie(s); aligned so there is no drift (for instance using [Image Stabalizer](https://imagej.net/plugins/image-stabilizer))
* A .tiff or .tiffs of the full fluorescent movie(s); aligned so there is no drift
* A .tiff or .tiffs of only the first frame of the phase movie(s)
* The output of either [Oufti](www.oufti.org) or [MicrobeJ](https://microbej.com) segmentation of the bdelloplasts in the first frame of the phase contrast movies. This pipeline was tested with Oufti.

### Software requirements

The pipeline uses the following software:

* [R => 4.2.2](https://cran.r-project.org) 
* ImageJ2 with the [Trackmate Plugin](https://imagej.net/plugins/trackmate/) (standard in the [FIJI distribution of ImageJ](https://imagej.net/software/fiji/))
* BFconvert from BioFormats bftools. (Download & unzip bftools.zip [here](https://downloads.openmicroscopy.org/bio-formats/5.5.2/artifacts/bftools.zip))

### Other notes/requirements

* It is convenient but not required to run R chunks from [Rstudio](https://posit.co/products/open-source/rstudio/)
* The pipeline was used & tested in Windows 10

### The parameter file

This file `parameters.yml` contains all the settings the pipeline needs to know to run properly. This is changed for each experiment. You can find a `parameters.yml` file in each Experiments subfolder. This file contains the paths to the required softwares & datasets and the parameters used to compute the output data. When working with new data, make sure to at least change the folder paths & check the Microscopy settings (pixel to micron conversion & minutes per frame). 

## Experiments

The *experiments* folder contains the data output & visualization code used in **Santin *et al.* ** This is divided in experiments exploring the *Bdellovibrio* **growth rate** (folder *Growth_rates*) and the number of **Progenies** and **Origins of Replication** (folder *Progenies_&_cell_cycle*). The code files containing the downstream statistical analysis and visualization are saved in the main Experiments folder.



- Pipeline
- Experiments
  - GR
    - LB E. coli
    - M9 E. coli
  - Progenies & cell cycle
    - LB E. coli
    - M9 E. coli
    - LB Citrobacter
      - Parameters.yml
      - *original pipeline outputs*.rmd
      - *compiled and summarized data*.csv
      - Plots.Rmd
    
    
Test_push
