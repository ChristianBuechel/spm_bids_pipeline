SPM based data import, preprocessing and analyses for data in BIDS format.  
Includes a function to download the data from the ISN DICOM data base to the project directory  

THIS IS A MAJOR REVISION, it is very likely that you have adapt your scripts if you have used a previous version


**The data are organized like this:**  
# Organisation  
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

Before you start you need to create a **get_study_specs.m** in your "scripts" folder. There is an example inside the repository "get_study_specs.txt"  
Copy it to "scripts" edit it and rename it to "get_study_specs.m"  

Then you should create a wrapper file such as "do_em_all.m" that will sequentially call the subfunctions in the repository  
There is a commented example in the repo (with ext .txt)

The new version of the toolbox not only uses a nonlinear coregistration as in CAT12, but now also incorporates correction with fieldmaps using VDM files. 
Both workflows are described in **do_em_all.txt**

This toolbox now also includes quality control routines found in \qc by Alexandra Tinnermann
See full_diagn_qc.txt for how to call these functions. 

# New features
1) The new version of the toolbox now also incorporates correction with fieldmaps using VDM files.

2) The new version uses an extension to participants.tsv using a column labeled 'scan_id' (the old import still works)
These contain the individual PRISMA numbers. If there are more sessions, seperate the PRISMA numbers by spaces NOT tabs !
In addition you can specify series that should NOT be imported by adding them in () with NO space in between e.g.

participant_id	sex	age	scan_id
sub-09	F	20	25879 25878(13)

imports sub-09, with sessions 25879 and 25878. Series 13 in 25878 will be omitted

3) all major functions were now renamed to start with sbp_
4) Can now do PPI analyses (Thanks to Lars Geuter for implementing this)
5) Can now incorporate physio (retroicor type) correction with the TAPAS toolbox (Thanks to Leo Asan for implementing this)
    

# Physiological Noise correction
- the function sbp_create_physio can be integrated into do_em_all in the same way as the other preprocessing/noise creation functions.

IMPORTANT: Your smr files need to be stored as such:
path.baseDir/physio/smr/sub-XX_ses-YY_physio.smr

= OPERATIONS (functions found in utils) =:
1. convert_physio_smr: converts raw SMR to mat-files that are readable for cleaning and regressor creation. via https://github.com/tjrantal/Spike-smr-reader
2. clean_physio: cleans mats (divide into runs, clean excess scanner pulses / add scanner pulses where necessary) and saves as BIDS-compatible run-wise tsv files with jsons in func dirs
3. create_tapas: Uses tapas toolbox to create txt file of 18 physiological noise regressors in func dirs of subjects.

= OUTPUTS =:
-BIDS-compatible run-wise tsv files of physiological data with jsons in func dir of each sub
-txt file of 18 physiological noise regressors in func dirs of each sub.
-info_physio mat file containing some info coming from tapas toolbox, in func dir of each sub

= REQUIREMENTS =:
1. TAPAS PhysIO Toolbox
	Download:
    **HAS TO BE DOWNLOADED**

    https://github.com/translationalneuromodeling/tapas/releases
    go to tapas_master and type tapas_init() in command line
    you need to do this only once. It will be set up for good in your spm12 after

 