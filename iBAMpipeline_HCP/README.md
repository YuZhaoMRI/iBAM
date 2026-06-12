# iBAM Pipeline for HCP Dataset (`iBAMpipeline_HCP`)

This directory contains the pipeline optimized for running the **implicit Brain Activity Mapping (iBAM)** framework on the Human Connectome Project (HCP) dataset. The code is structured and written with a **batch-processing** mindset to handle large-scale neuroimaging datasets efficiently.

---

##  CRITICAL PRE-REQUISITES 

Before executing any scripts in this directory, you **must** fulfill the following requirements:

1. **Read the Root Documentation:** Thoroughly read the `README.md` file in the parent (root) directory of this repository.
2. **Path Configuration:** Ensure all required dependencies, toolboxes (SPM12, NIfTI tools), and the **`Masks`** folder are properly configured and added to your MATLAB path.

---

##  Directory Structure & File Naming Conventions

This pipeline relies on preset, hard-coded file organization and naming rules tailored to the HCP dataset:

* **Directory Structure:** The code expects a specific fMRI file structure that mirrors the directory tree layout of downloaded HCP data.
* **File Naming Warning:** Please note that the HCP data used to develop this pipeline underwent **further downstream processing** after initial download. Consequently, the actual fMRI data filenames inside the scripts **do not match** the raw, default filenames provided directly by the HCP database.

>  **Action Required:** Before running the scripts, you must inspect the filenames specified inside the code. Your local fMRI data filenames and your folder structures **must match exactly** with what is hard-coded in the scripts for the pipeline to run successfully.

---

##  Execution Order (Recommended Workflow)

To process your data, execute the batch scripts strictly in the following sequence:

1. **`generateModeIntensityMap_SW.m`** Calculates the primary Mode Intensity (MI) maps from the task fMRI data.

2. **`generateNullModels_SW.m`** Generates the reference null models using the resting-state fMRI data.

3. **`iBAMstatisticalInference.m`** Performs statistical inference by evaluating the MI maps against the generated null models to output the final activation Z-maps.

### Post-Processing & Downstream Analysis

* After successfully executing the three core steps above, you can run any auxiliary analysis or metric-gathering functions. These are the scripts whose names begin with the prefix **`calculate...`** (e.g., calculating region-of-interest statistics or block averages).
