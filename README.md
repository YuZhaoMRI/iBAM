# iBAM MATLAB Package

This repository provides the core MATLAB implementation of **implicit Brain Activity Mapping (iBAM)**, a framework designed to detect brain neural activations from task fMRI data.

## Overview

iBAM is intended to identify brain activations that may not be fully captured by conventional task timing-based activation mapping. The core program uses task fMRI data as the primary input for brain activation detection. Resting-state fMRI data are also required to construct the null models used for statistical inference.

## Input data

The program requires two types of fMRI data:

1. **Task fMRI data**  
   Used for detecting brain neural activations.

2. **Resting-state fMRI data**  
   Used for constructing null models.

Before running iBAM, all fMRI data should first be preprocessed using either:

- **fMRIPrep**, or
- **HCP minimal preprocessing pipelines**.

After standard preprocessing, the data should further undergo the **ICE-based filtering** procedure described in the associated manuscript.

## Dependencies

This package is implemented in MATLAB and requires the following external toolboxes.

### SPM12

The program calls functions from **SPM12** during execution. Please install SPM12 and add it to the MATLAB path before running the code.

SPM12 is available at:  
https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

### Tools for NIfTI and ANALYZE image

The program also uses Jimmy Shen’s **Tools for NIfTI and ANALYZE image**, version 1.26.0.0, for reading and writing NIfTI-format neuroimaging data.

The toolbox is available from MATLAB File Exchange:  
https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

Please download the toolbox and add it to the MATLAB path before running iBAM.

## Recommended workflow

A typical iBAM workflow consists of the following steps:

1. Preprocess task and resting-state fMRI data using fMRIPrep or the HCP minimal preprocessing pipelines.
2. Apply the ICE-based filtering procedure described in the manuscript.
3. Run `calculateMl_SWrun3.m` to calculate the MI map from the task fMRI data.
4. Run `generateNullModels_SWrun3.m` to generate null models from the resting-state fMRI data.
5. Run `calculate_Zmap.m` to perform statistical inference and generate the final Z map.

## Core scripts

### `calculateMl_SWrun.m`

This script calculates the MI map from the input task fMRI data. The resulting MI map represents the spatial distribution of mode intensity estimated by the iBAM framework and serves as the primary activation-related map for subsequent statistical inference.

### `generateNullModels_SWrun.m`

This script generates null models using resting-state fMRI data. These null models provide the reference distribution required to statistically evaluate whether the observed MI values from task fMRI data exceed the expected background level.

### `calculate_Zmap.m`

This script performs statistical inference by comparing the task-derived MI map with the null models. The output is a Z map that quantifies statistically significant brain activations detected by iBAM.

## Notes

This repository contains the core code used for implementing the iBAM framework. Users should refer to the associated manuscript for the detailed methodological description, including preprocessing, ICE-based filtering, MI map calculation, null model construction, and statistical inference procedures.

