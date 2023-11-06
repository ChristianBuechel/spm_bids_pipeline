function y = get_roi_ts(xY, V_epi, V_wmean_scan, N_back_scan,V_part)
[~, XYZmm, ~] = spm_ROI_ind(xY, V_wmean_scan); % get the voxels in the image

if (isfield(xY,'thresh')) & strcmp(xY.def,'mask')
    XYZvoxM       = inv(V_part.mat)*[XYZmm; ones(1,size(XYZmm,2))];
    XYZvoxM       = XYZvoxM(1:3,:); %map mm to vox
    mask          = spm_get_data(V_part,XYZvoxM); %extract voxel values
    prune_ind     = find(mask>xY.thresh); % threshold
    XYZmm         = XYZmm(:,prune_ind); % kickout voxels below threshold
end

[oXYZvox, ~]  = transform_back(XYZmm, V_wmean_scan, V_epi(1), N_back_scan);

y         = spm_get_data(V_epi,oXYZvox); % get relevant EPI time-series


