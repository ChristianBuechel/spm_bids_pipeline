SPM based data import, preprocessing and analyses for data in BIDS format  
Includes a function to download the data from the ISN DICOM data base to the project directory  

The data are organized like this:  

# Organisation  
c:\peep - the top data directory  
inside you will get a directory \derivatives  
The data import routine (import_data_ses.m) will put all images in /derivatives/spm_preprocessing but will also make a  
gzipped copy in the main study directory  
  
You can clone this repository to any directory, but should keep it as is  

**Do not edit files inside spm_bids_pipeline**  
  
It is advisable to have a "scripts" directory somewhere where your wrapper functions are located  
In this wrapper function make sure to have an   
**addpath(genpath('your_path_2_spm_bids_pipeline'))**

inside derivatives you have spm_preprocessing where all the data goes using \sub-01 etc  
inside sub-01 you have \anat and \func   

Before you start you need to create a **get_study_specs.m** in your "scripts" folder. There is an example inside the repository "get_study_specs.txt"  
Copy it to "scripts" edit it and rename it to "get_study_specs.m"  

Then you should create a wrapper file such as "do_em_all.m" that will sequentially call the subfunctions in the repository  
There is a commented example in the repo (with ext .txt)

This now also includes quality control routines found in \qc by Alexandra Tinnermann
See full_diagn_qc.txt for how to call these functions. 
