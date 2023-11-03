SPM based preprocessing and analyses for data in BIDS format  
The data are organized like this:  

# Organisation  
c:\peep - the top data directory  
inside you need \derivatives and \rawdata (in the latter nothing happens)  
insde derivatives you have spm_preprocessing where all the data goes using \sub-01 etc  
inside sub-01 you have \anat and \func  

You also need to create a **get_study_specs.m**. There are a few examples to help you   

**The first are the path definitions and need to be adapted**  
%% path definitions  
path.baseDir     = 'd:\peep\';  
path.templateDir = 'd:\peep\spm_bids_pipeline\templates\'; % the templates come with the git repository, so should be here  
  
path.rawDir          = fullfile(path.baseDir, 'rawdata'); % that should be the folder with the not-to-touch data sets  
path.derivDir        = fullfile(path.baseDir, 'derivatives');  
path.preprocDir      = fullfile(path.baseDir, 'derivatives', 'spm_preprocessing');  
path.firstlevelDir   = fullfile(path.baseDir, 'derivatives', 'spm_firstlevel');  
path.secondlevelDir  = fullfile(path.baseDir, 'derivatives', 'spm_secondlevel');  
  
**Don't touch the following unless you know what you are doing**  
%% various predefined names (do not change)  
vars.skullStripID    = 'skull-strip-T1.nii';  
vars.T1maskID        = 'brain_mask.nii';  
vars.templateID      = 'cb_Template_%d_Dartel.nii';  
vars.templateT1ID    = 'cb_Template_T1.nii';  
%vars.groupMaskID     = 'brainmask.nii';  
vars.groupMaskID     = 'neuromorphometrics.nii';  

**the next section needs to be edited to match your data**  
vars.max_procs   = 12 %How many physical cores of your CPU you can use  
vars.task        = 'peep'; %Name of your task in all the filenames  
vars.nRuns       = 4; %Number of Runs  
vars.nSess       = 2; %Number of Sessions  
% get info for slice timing correction  
vars.sliceTiming.so    = [1722,1662,1602,1543,1485,1425,1365,1305,1248,1188,1127,1067,1010,950,890,832,772,712,653,595,535,475,415,358,298,237,177,120,60,0,...  
    1722,1662,1602,1543,1485,1425,1365,1305,1248,1188,1127,1067,1010,950,890,832,772,712,653,595,535,475,415,358,298,237,177,120,60,0]; %vector of acquistion times in ms  
vars.sliceTiming.tr       = 1.8;  %TR in s  
vars.sliceTiming.nslices  = 60; %Number of slices, should be size(vars.sliceTiming.so,2)  
vars.sliceTiming.refslice = 900; %Reference slice for slice timing  
  
  
# Examples  
## do_em_all_basic   
  
The workflow of a basic brain EPI+T1 preprocessing  
   
Final EPI images for analysis: __rasub*__  
  
## do_em_all_realign_two_step  
as above, but realignment of the brain is done in two steps:   
Initailly a single run is performed to create a mean  
After all the masks have been created the realignment (2-pass) is done with a brainmask eg. to discard eye movements  
Final EPI images for analysis: __rasub*__  
  
## do_em_all_bs  
realignment is done specificlly for the brainstem  
Final EPI images for analysis: __brasub*__  

## do_em_all_bs_physio  
realignment is done specificlly for the brainstem but before motion is estimated physio and CSF regressors are removed from the data  
Final EPI images for analysis: __pbrasub*__  
  
  
# TODO:  
Currently task TSV files are interpreted as times in TRs not seconds  
