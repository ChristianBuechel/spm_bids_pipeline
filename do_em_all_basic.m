function do_em_all_basic
% test the new batch system 

% define subjects that should be included
all_sub_ids  = [6]; 

% perform slice timing correction
%slice_timing(all_sub_ids);

% perform full realignment (for slice time corrected data)
%realign(all_sub_ids);

% perform non-linear coregistartion to get from EPI space to T1 space
% nonlin_coreg(all_sub_ids);

% perform skull stripping
% skullstrip(all_sub_ids);

% perform spatial normalization of T1 to Template space using DARTEL
% create_dartel(all_sub_ids);

% create various warp fields to map from EPI --> T1 --> Template space
%create_trans(all_sub_ids);

% Create a mnask for the 1st level GLM
%create_mask(all_sub_ids);

% Warp a few images to Template space (skullstrip and mean EPI)
% warp_images(all_sub_ids);

% Create a mean skullstrip and a mean of mean EPIs
% create_means(all_sub_ids);

% Create 6 WM and 6 CSF noise regressors (and 6 regressors picking up noise
% from the posterior part of the lateral ventricle)
% create_noise(all_sub_ids);
  
% use coreg to display images for QC 
% check_images(all_sub_ids);

% Display realignment parameters (same as SPM) and check 1st image from all
% runs
check_realign(all_sub_ids);

% do analyses
analyses(all_sub_ids);

end


