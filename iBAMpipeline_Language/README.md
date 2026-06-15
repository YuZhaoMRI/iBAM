# iBAM Pipeline for fMRIprep-processed Dataset (`iBAMpipeline_fMRIprep`)

This directory contains the pipeline optimized for running the **implicit Brain Activity Mapping (iBAM)** framework on datasets processed via **fMRIprep**. The code is structured and written with a **batch-processing** mindset to handle large-scale neuroimaging datasets efficiently.

---

##  CRITICAL PRE-REQUISITES

Before executing any scripts in this directory, you **must** fulfill the following requirements:

1. **Read the Root Documentation:** Thoroughly read the `README.md` file in the parent (root) directory of this repository.
2. **Path Configuration:** Ensure all required dependencies, toolboxes (SPM12, NIfTI tools), and the **`Masks`** folder are properly configured and added to your MATLAB path.
3. **Special Note for In-House/Custom Datasets:** * If you are using **your own collected dataset**, it is highly recommended that you use this pipeline for your analysis.
   * Before running this program, you **must** preprocess your custom data strictly following the preprocessing steps described in our paper.

---

##  Directory Structure & File Naming Conventions

This pipeline relies on preset file organization and naming rules tailored to the dataset configuration:

* **Directory Structure:** The code expects a specific fMRI file structure that matches **how it is configured in the code**.
* **File Naming Warning:** The fMRI data filenames inside the scripts match **how they are configured in the code**.

>  **Action Required:** Before running the scripts, you must inspect the filenames specified inside the code. Your local fMRI data filenames and your folder structures **must match our template script paths exactly** for the pipeline to run successfully.

---

##  Execution Order (Recommended Workflow)

To process your data, execute the batch scripts strictly in the following sequence:

1. **`generateModeIntensityMap_SW.m`** Calculates the primary Mode Intensity (MI) maps from the task fMRI data.

2. **`generateNullModels_SW.m`** Generates the reference null models using the resting-state fMRI data.

3. **`performStatisticalInference.m`** Performs statistical inference by evaluating the MI maps against the generated null models to output the final activation Z-maps.

### Post-Processing & Downstream Analysis

* After successfully executing the three core steps above, you can run any auxiliary analysis or metric-gathering functions. These are the scripts whose names begin with the prefix **`calculate...`** (e.g., calculating region-of-interest statistics or block averages).
