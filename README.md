# iBAM MATLAB Package

This repository provides the core MATLAB implementation of **implicit Brain Activity Mapping (iBAM)**, a framework designed to detect brain neural activations from task fMRI data.

## Overview

The unique methodological contribution of iBAM in brain science resides in its capacity to map brain activations that elude explicit experimental manipulations. Conventional General Linear Model (GLM)-based activation mapping is structurally constrained by two limitations: first, its reliance on strict temporal coupling between stimuli and neural dynamics renders it blind to autonomous, endogenous neural processes that are partially decoupled from external stimuli; second, its dependence on aggregating signal segments across multiple repeated trials to counteract noise in activation detection "washes out" the trial-varied brain dynamics in complex cognition. 

iBAM transcends these constraints through a dual technical breakthrough:

1. **Task Timing-Independent Modeling**: By identifying latent signal modes, iBAM emancipates functional activation detection from traditional task timing-based fMRI signal modeling.
2. **Single-Trial Sensitivity:** Its single-trial sensitivity allows for the capture of transient, trial-specific neural activations without the "washing-out" effect inherent to multi-trial averaging.

The core program uses task fMRI data as the primary input for brain activation detection, while resting-state fMRI data are utilized to construct the null models required for statistical inference.

## Input data

The program requires two types of fMRI data:

1. **Task fMRI data** Used for detecting brain neural activations.

2. **Resting-state fMRI data** Used for constructing null models.

Before running iBAM, all fMRI data should first be preprocessed using either:

- **fMRIPrep**, or
- **HCP minimal preprocessing pipelines**.

After standard preprocessing, the data should further undergo the **ICE-based filtering** procedure described in the associated manuscript.

## Dependencies

This package is implemented in MATLAB and requires the following external toolboxes and folders to be added to the MATLAB path.

### Required Path Configurations

- **Masks Folder** Please ensure that the **`Masks`** folder included in this repository is added to your MATLAB path before running any scripts.

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

1. **For HCP Data Processing:** If you are running the pipeline on HCP (Human Connectome Project) data, please use the code provided in the **`iBAMpipeline_HCP`** directory. Carefully read the instructions within that folder before proceeding. Note that to accommodate the specific file organization of the HCP dataset and optimize computational efficiency, certain parts of this code include hard-coded parameters. Therefore, if you intend to run your own custom (non-HCP) data, you should either modify this code accordingly or refer to the alternative pipeline below.

2. **For Language Study Data & Custom Datasets:** If you are analyzing data related to language studies, please use the code in the **`iBAMpipeline_Language`** directory. Furthermore, if you have your own independently collected datasets that need to be processed, it is **strongly recommended** to use the scripts in this folder as your primary pipeline.

## Notes

This repository contains the core code used for implementing the iBAM framework. Users should refer to the associated manuscript for the detailed methodological description, including preprocessing, ICE-based filtering, MI map calculation, null model construction, and statistical inference procedures.
