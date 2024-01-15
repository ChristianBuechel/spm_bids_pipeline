SPM based preprocessing and analyses for data in BIDS format  
Includes a function to download the data from the DICOM data base to the project directory

The data are organized like this:  

# Organisation  
c:\peep - the top data directory  
inside you need \derivatives and \rawdata (in the latter nothing happens)
here you also clone this repository and should create a /scripts directory where your *.m functions go
**Do not edit the spm_bids_pipeline**

inside derivatives you have spm_preprocessing where all the data goes using \sub-01 etc  
inside sub-01 you have \anat and \func  

Before you start you need to create a **get_study_specs.m** in your /scripts folder. Here is an example to help you   


**The first section deals with data import from the DICOM data base**  

import.prisma        = {{22147, 22165},{22197, 22214}}; % translates to PRISMA_nnn What is in curly brackets belongs to 1 volunteer (ie 2 sessions) 
import.prisma_no     = [2             ,3             ]; % assign BIDS subject numbers
                                    
import.user          = 'buechel';
import.server        = 'aither.nin.uke.uni-hamburg.de';

import.data(1).dir        = 'func'; %needs to be BIDS conform
import.data(1).type       = 'epi';
import.data(1).seq        = 'ninEPI_bold_v12B, 1.5mm3, mb3 '; %protocol name (trailing space makes it unique) 
import.data(1).cond       = 'n > 600'; % heuristic to get only valid runs (e.g. more than 1000 volumes)

import.data(2).dir        = 'anat'; % valid BIDS dir name
import.data(2).type       = 'T1w'; % valid BIDS file name
import.data(2).seq        = 'mprage, HR64 ';
import.data(2).cond       = 'n == 240'; % heuristic to get only valid runs (e.g. exactly 240 slices)
 
% import.data(3).dir        = 'fmap';
% import.data(3).type       = 'phasediff';
% import.data(3).seq        = 'gre_field_map, 2mm ';
% import.data(3).cond       = 'n == 40';
% 
% import.data(4).dir        = 'fmap';
% import.data(4).type       = 'magnitude';
% import.data(4).seq        = 'gre_field_map, 2mm ';
% import.data(4).cond       = 'n == 80';

import.dummies            = 3; % how many dummies will be deleted in the 4D epi file


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
  
## do_em_all_realign_two_step  (RECOMMENDED)
as above, but realignment of the brain is done in two steps:   
Initailly a single run is performed to create a mean  
After all the masks have been created the realignment (2-pass) is done with a brainmask eg. to discard eye movements  
Final EPI images for analysis: __rasub*__  
  
## do_em_all_bs  
realignment is done specificlly for the brainstem  
Final EPI images for analysis: __brasub*__  

## do_em_all_bs_physio  (EXPERIMENTAL)
realignment is done specificlly for the brainstem but before motion is estimated physio and CSF regressors are removed from the data  
Final EPI images for analysis: __pbrasub*__  
  
  
