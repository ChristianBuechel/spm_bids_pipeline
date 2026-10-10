function extract_roi(epi_files, wmean_file, back_file, roi_def, ex_var, fname)
% Identifies voxels in rois (sphere) 
% roi is a n x 4 matrix with columns x,y,z radius
% epi_files is a char array
% EPI and segments are in different spaces (e.g. through nonlin coreg or nonlin coreg + Dartel)
% Therefore one has to supply not only the epi_files, but also an epi file normalzed to the space of the segments
% (--> cmean* for nonlin coreg, wmean* for nonlin coreg + Dartel)
% in addition one needs to provide the deformation (back_file) FROM epi to segments (--> y_mean* for nonlin coreg, y_epi_2_template* for nonlin coreg + Dartel)
% voxels in the rois are then read from the EPIs and PCs (ex_var) extracted
% if ex_var <0 the number of PCs that explain ex_var of the total variance
% are extracted
% if exvar>1 ex_var PCs are returned
% the time-series are saved to roi_noise.mat in the directrory of the
% EPIs so session specific, we might have to change this for BIDS
% ROI voxels outside the field of view of the EPI (slab acquisitions) are ignored; if fewer than
% max(2,ceil(ex_var)) voxels are left, a warning and a message are given, roi(i).data is empty (no
% regressors for this ROI) and the function goes on with the next ROI.
% roi(i).n_voxels = [voxels used, voxels in the ROI]. Whoever uses the file has to handle empty roi(i).data.

%__________________________________________________________________________
%
n_rois = size(roi_def,1);
[epi_dir,~,~]    = fileparts(epi_files(1,:));
V_epi            = spm_vol(epi_files);
V_wmean_scan     = spm_vol(wmean_file);
N_back_scan      = nifti(back_file);

for i=1:n_rois
    name       = num2str(i);    
    xY(i).name = name;
    xY(i).spec = roi_def(i,4);
    xY(i).def  = 'sphere';
    xY(i).xyz  = roi_def(i,1:3)';
    
    yy          = get_roi_ts(xY(i), V_epi, V_wmean_scan, N_back_scan);
    % voxels of the ROI that are inside the field of view of the EPI (outside: NaN) and not empty (0)
    % (a slab acquisition does not cover all ROIs)
    ok          = all(isfinite(yy),1) & any(yy,1);
    n_min       = max(2,ceil(ex_var));
    roi(i).name = name;
    roi(i).n_voxels = [nnz(ok) size(yy,2)];     % voxels used, voxels in the ROI
    if nnz(ok) < n_min
        warning('extract_roi:outsideFOV','ROI %s (%.1f %.1f %.1f, radius %g mm): only %d of %d voxels are inside the field of view of the EPI, %d needed: NO regressors for this ROI (roi(%d).data is empty).', ...
            name,xY(i).xyz(1),xY(i).xyz(2),xY(i).xyz(3),xY(i).spec,nnz(ok),size(yy,2),n_min,i);
        fprintf('extract_roi: ROI %s outside the field of view of %s (%d of %d voxels), roi(%d).data is empty\n',name,spm_file(V_epi(1).fname,'filename'),nnz(ok),size(yy,2),i);
        roi(i).data = [];
        continue
    end
    if nnz(ok) < size(yy,2)
        fprintf('extract_roi: ROI %s: %d of %d voxels are inside the field of view of the EPI, the others are ignored\n',name,nnz(ok),size(yy,2));
    end
    roi(i).data = get_pcs(yy(:,ok),ex_var); % get principal components
end


%-Save
%==========================================================================
roi_f = spm_file(V_epi(1).fname,'ext','.mat','prefix',fname);
save(roi_f,'roi');
