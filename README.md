# Prey Size Analysis

This depository contains the code & datasets used in **Santin *et. al*,  Modulation of prey size reveals adaptability and robustness in the cell cycle of an intracellular predator., 2022**. In this work, time lapse microscopy was used to monitor the *Bdellovibrio* cell cycle at the single cell level. The depository can be divided in two parts: the *analysis pipeline* used to obtain the data, and the *experiments*, where the datasets from different experiments are gathered. The experiments folder also contains the scripts used for the visualization of each dataset.

## Pipeline

The pipeline is written as an R markdown file, where the user can choose to run the complete script or run the analysis chunk-by-chunk.

### Data curation: before running the pipeline

The pipeline contains code to process movies of bdellovibrio containing fluorescent foci and/or cytoplasmic signal & uses phase-contrast to detect the bdelloplasts. As input it needs folders containing

* A .tiff or .tiffs of the full phase contrast movie(s); aligned so there is no drift 
* A .tiff or .tiffsof the full fluorescent movie(s); aligned so there is no drift
* A .tiff or .tiffs of only the first frame of the phase movie(s)
* The output of either Oufti (www.oufti.org) or MicrobeJ (www.microbej.com) segmentation of the bdelloplasts in the first frame of the phase contrast movies.

### Software requirements


### The parameter file

This file contains all the settings the pipeline needs to know to run properly. This is changed for each experiment. 

## Experiments

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
