function y = get_roi_ts(xY, V_epi, V_wmean_scan, N_back_scan,V_part)
% function y = get_roi_ts(xY, V_epi, V_wmean_scan, N_back_scan,V_part)
% universal function that returns voxel time series from the EPIs based on
% the xY ROI definition. IMportantly a deformation field is specified that
% should map the mm coordinates from the ROI structure back to the EPI
% space
% although V_part is part of xY (xY.spec), V_part has the eroded data in
% RAM which is crucial for the threshold masking

[~, XYZmm, ~] = spm_ROI_ind(xY, V_wmean_scan); % get the relevant ROI voxels in mm

if (isfield(xY,'thresh')) & strcmp(xY.def,'mask') % helpful if we want only voxels in mask above a certain threshold
    XYZvoxM       = inv(V_part.mat)*[XYZmm; ones(1,size(XYZmm,2))]; % mm--> vox
    XYZvoxM       = XYZvoxM(1:3,:); % remove ones 
    mask          = spm_get_data(V_part,XYZvoxM); % extract voxel values
    prune_ind     = find(mask>xY.thresh); % threshold
    XYZmm         = XYZmm(:,prune_ind); % kickout voxels below threshold
end

[oXYZvox, ~]  = transform_back(XYZmm, V_wmean_scan, V_epi(1), N_back_scan);

y         = spm_get_data(V_epi,oXYZvox); % get relevant EPI time-series


