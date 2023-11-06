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
    roi(i).data = get_pcs(yy,ex_var); % get principal components
    roi(i).name = name;
end


%-Save
%==========================================================================
roi_f = spm_file(V_epi(1).fname,'ext','.mat','prefix',fname);
save(roi_f,'roi');
