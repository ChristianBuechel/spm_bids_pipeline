function epi_def_from_vaso(y_vaso,mean_epi,ref_bold,y_epi)

% function epi_def_from_vaso(y_vaso,mean_epi,ref_bold,y_epi)
% deformation field of the mean EPI of the main runs, from the deformation field of the vaso T1-like
% image (nonlinear coregistration, see sbp_coreg_nlin_vaso_T1).
%   y_vaso     deformation field (y_*.nii, written by SPM's normalisation: for each voxel of the T1
%              / template grid the coordinates [mm] in the vaso image)
%   mean_epi   mean image of the main EPI runs (mean<P><epi>)
%   ref_bold   mean bold image of the vasoref run (defines the space of the vaso images)
%   y_epi      output: the same kind of deformation field for the mean EPI (y_mean<P><epi>)
% The vasoref and the main runs have the same sequence, geometry and distortion, but the head
% may have moved between them: a rigid body registration (NMI) of the mean EPI onto the mean
% bold of the vasoref run gives the matrix M with  (position in the EPI) = M * (position in the vaso
% image), the coordinates in y_vaso are transformed with M. Headers are not changed.

VG = spm_vol(ref_bold);  VG = VG(1);
VF = spm_vol(mean_epi);  VF = VF(1);
flags = struct('cost_fun','nmi','sep',[4 2],'tol',[0.02 0.02 0.02 0.001 0.001 0.001], ...
               'fwhm',[7 7],'graphics',0);
x = spm_coreg(VG,VF,flags);
M = spm_matrix(x(:)');                            % as in spm_run_coreg: M\mat(EPI) is the aligned header
fprintf('epi_def_from_vaso: rigid transformation EPI -> vaso [%s]\n',sprintf('%.2f ',x));

copyfile(y_vaso,y_epi);
N  = nifti(y_epi);                                % x y z 1 3 (the copy is overwritten below)
Nv = nifti(y_vaso);
Y  = cell(1,3);
for c = 1:3, Y{c} = double(Nv.dat(:,:,:,1,c)); end
for c = 1:3
    N.dat(:,:,:,1,c) = single(M(c,1)*Y{1} + M(c,2)*Y{2} + M(c,3)*Y{3} + M(c,4));   % NaN stay NaN
end
end
