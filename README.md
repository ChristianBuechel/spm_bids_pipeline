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
**addpath('path_2_som_bids_pipeline')**

inside derivatives you have spm_preprocessing where all the data goes using \sub-01 etc  
inside sub-01 you have \anat and \func   

Before you start you need to create a **get_study_specs.m** in your "scripts" folder. There is an example inside the repository "get_study_specs.txt"  
Copy it to "scripts" edit it and rename it to "get_study_specs.m"  

  
Then you should create a wrapper file such as "do_em_all.m" that will sequentially call the subfunctions in the repository  
There atre some examples in the repo (all with ext .txt)
# Examples  
## do_em_all_basic   
  
The workflow of a basic brain EPI+T1 preprocessing  

## do_em_all_realign_two_step  (RECOMMENDED)  
realignment of the brain is done in two steps:  
Initailly a simple run is performed to create a mean image  
After all the masks have been created the realignment (2-pass) is done again, but with a brainmask eg. to discard eye movements  

% perform slice timing correction  
slice_timing(all_sub_ids);  
  
% perform initial realignment (one pass, write mean only)  
realign_1_2(all_sub_ids);  
  
% perform non-linear coregistartion to get from EPI space to T1 space  
nonlin_coreg(all_sub_ids);  
  
% perform skull stripping  
skullstrip(all_sub_ids);  
  
% perform spatial normalization of T1 to Template space using DARTEL  
create_dartel(all_sub_ids);  
  
% create various warp fields to map from EPI --> T1 --> Template space  
create_trans(all_sub_ids);  
  
% Create a mnask for the 1st level GLM  
create_mask(all_sub_ids);   
  
% now we repeat the realignemnt but with a brain mask and write out all images  
realign_2_2(all_sub_ids);  
  
% Warp a few images to Template space (skullstrip and mean EPI)  
warp_images(all_sub_ids);  
  
% Create a mean skullstrip and a mean of mean EPIs for all_sub_ids  
create_means(all_sub_ids);  
  
Final EPI images for brain analysis: __rasub*__  
  
  
  
## do_em_all_bs  
additional realignment is done specificlly for the brainstem  
do this AFTER you have done the above ie preprocessing for the brain  
  
% create subject specific brainstem masks  
create_bs_mask(all_sub_ids);  
% does realignment with a weight of 1 over the braistem and zeros elsewhere  
realign_bs(all_sub_ids);  
  
Final EPI images for brainstem analysis: __brasub*__  

