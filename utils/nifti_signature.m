function sig = nifti_signature(files)

% function sig = nifti_signature(files)
% geometry of 3D NIfTI volumes for classify_geometry: one row per file
%   [n_slices, slice thickness (mm), in-plane voxel size (mm)]
% (3rd voxel dimension = slice thickness incl. gap, as written by the SPM DICOM
% conversion of mosaic EPIs)

files = cellstr(files);
sig   = nan(numel(files),3);
for ii = 1:numel(files)
    V = spm_vol(files{ii});
    if numel(V) > 1, V = V(1); end
    vox = sqrt(sum(V.mat(1:3,1:3).^2,1));
    sig(ii,:) = [V.dim(3) vox(3) vox(1)];
end
end
