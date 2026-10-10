SPM based data import, preprocessing and analyses for data in BIDS format.  
Includes a function to download the data from the ISN DICOM data base to the project directory  

> **A major revision (Version 3) is in preparation on the branch `v3-overhaul`.**  
> It is not finished yet (the analysis module `sbp_analyses.m` still has to be refactored) and it is **not compatible** with this version (file names, `get_study_specs.m`).  
> Please keep using this branch (`main`, tag `Version-2.1`) until Version 3 has been merged.  

THIS IS A MAJOR REVISION, it is very likely that you have adapt your scripts if you have used a previous version

## Organisation  

**The data are organized like this:**  

c:\lcpa - the top data directory  
inside you will get a directory \derivatives  
The data import routine (sbp_import_data.m) will put all images in /derivatives/spm_preprocessing but will also make a  
gzipped copy in the main study directory  
  
You can clone this repository to any directory, but should keep it as is  

**It is best not to edit files inside spm_bids_pipeline**  
  
It is advisable to have a "scripts" directory somewhere where your wrapper functions are located  
In this wrapper function make sure to strat with an   
**addpath(genpath('your_path_2_spm_bids_pipeline'))**

inside derivatives you have spm_preprocessing where all the data goes using \sub-01 etc  
inside sub-01 you have \anat and \func (and \fmap if specified)  

Before you start you need to create a get_study_specs.m in your "scripts" folder. There is an example inside the repository "get_study_specs.txt"  
Copy it to "scripts" edit it and rename it to "get_study_specs.m"  

Then you should create a wrapper file such as "do_em_all.m" that will sequentially call the subfunctions in the repository  
There is a commented example in the repo (with ext .txt)

The new version of the toolbox not only uses a nonlinear coregistration as in CAT12, but now also incorporates correction with fieldmaps using VDM files. 
Both workflows are described in **do_em_all.txt**

This toolbox now also includes quality control routines found in \qc by Alexandra Tinnermann
See full_diagn_qc.txt for how to call these functions. 

## New features

### 10/25

1) Multiple Analysis Definition: ``analysis`` is now a struct array instead of single struct, allowing multiple analyses to be run sequentially. This includes a major overhaul of the get_study_specs structure. All changes are backwards compatible

2) Overhauled PPI analysis: now includes reuse of existing base models, improved consistency checks across subjects, and better user interaction for managing corrupted analyses

3) Critical bug fix to physio regressor creation pipeline. 

### 06/25

1) The new version of the toolbox now also incorporates correction with fieldmaps using VDM files.

2) The new version uses an extension to participants.tsv using a column labeled 'scan_id' (the old import still works)
These contain the individual PRISMA numbers. If there are more sessions, seperate the PRISMA numbers by spaces NOT tabs !
In addition you can specify series that should NOT be imported by adding them in () with NO space in between e.g.

participant_id	sex	age	scan_id
sub-09	F	20	25879 25878(13)

imports sub-09, with sessions 25879 and 25878. Series 13 in 25878 will be omitted

3) all major functions were now renamed to start with sbp_
4) Can now do PPI analyses (Thanks to Lars Keuter for implementing this)
5) Can now incorporate physio (retroicor type) correction with the TAPAS toolbox (Thanks to Leo Asan for implementing this)
    

## Pipeline for physio regressor creation 


IMPORTANT: 
- you need to define path.physioDir, which is baseDir/physio.

- in PhysioDir, you need to have a folder 'mat' or a folder 'smr', depending on which format your file is.
- Your smr files need to be stored as such:
	path.physioDir/physio/smr/sub-XX_ses-YY_physio.smr
- Your mat files need to be stored as such:
	path.physioDir/physio/mat/sub-XX_ses-YY_physio.mat
	
- if you have session-wise files and wish to exclude certain runs, you need to add a json sidecar, which has the following structure:
{
  "includedRuns": [1, NaN, 2, 3, 4]
  }
 where 1-4 are the runs that are included, and NaN is the run you want to exclude.
 save this json sidecar with the same name and in the same location as your session file, with a .json extension.

= OPERATIONS (functions found in utils) =:
1. convert_and_clean_physio: flexibly reads in mat- or smr- session- or run-files, cleans them and converts them to BIDS-compatible TSV files and 
2. create_tapas: Uses tapas toolbox to create txt file of 18 physiological noise regressors in func dirs of subjects.


= OUTPUTS =:
-BIDS-compatible run-wise tsv files of physiological data with jsons in func dir of each sub
-txt file of 18 physiological noise regressors in func dirs of each sub.
-info_physio mat file containing some info coming from tapas toolbox, in func dir of each sub


= REQUIREMENTS =:

1. TAPAS PhysIO Toolbox
	Download:
    **HAS TO BE DOWNLOADED AND INITIALIZED**

    https://github.com/translationalneuromodeling/tapas/releases
    go to tapas_master and type tapas_init() in command line
    you need to do this only once. It will be set up for good in your spm12 after


## PPI Analysis

Psychophysiological interaction (PPI) analysis examines how functional connectivity between brain regions changes depending on psychological context.

### Requirements and Setup

- Every PPI analysis requires a base first-level model defined in `get_study_specs.m`. When first running the PPI analysis, this base model will be estimated. For subsequent PPIs using the same base model, it won't be estimated again.
  - Multiple PPIs (using different ROIs or conditions) can be run sequentially by defining multiple elements in the analysis struct array. These PPIs will use the same base model
  - **Important**: If some subjects don't have the base model yet, you must first estimate the base analysis for these subjects before running the PPI

- PPI Configuration: PPI settings (ROI definition, condition selection etc.) are defined in `get_study_specs.m`

### Best practise

- Verify that condition indices in your PPI specification reference valid conditions in your base model
- For two-pass ROI optimization, make sure the contrast vector matches your model dimensions