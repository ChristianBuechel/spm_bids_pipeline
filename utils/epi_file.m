function f = epi_file(epi,vars,kind)

% function f = epi_file(epi,vars,kind)
% name of an EPI related file, derived from the BIDS name of the raw EPI (epi, char or
% cellstr, e.g. .../sub-01_ses-01_task-x_run-01_bold.nii) and the preprocessing steps
% (see epi_prefix, vars.preproc.steps). P is the realign input prefix, e.g. 'an'.
%
%   kind             file                                   written by
%   'nordic_in'      <epi> with the prefix of the steps before nordic (input of sbp_nordic)
%   'slicetiming_in' <epi> with the prefix of the steps before slice timing (input of sbp_slice_timing)
%   'nordic'      n + <epi> (prefix of nordic input)     sbp_nordic
%   'slicetiming' a + <epi> (prefix of slice timing in)  sbp_slice_timing
%   'in'          P<epi>          images realigned etc.  (nordic / slice timing output)
%   'mat'         P<epi>.mat      sidecar with the realigned header
%   'rp'          rp_P<epi>.txt   realignment parameters
%   'res'         rP<epi>         resliced images (realign 2 / apply vdm)
%   'smooth_in'   (prefix of the steps before smooth, e.g. rn<epi>)    input of sbp_smooth
%   'smooth'      s + smooth_in   smoothed images (sbp_smooth), e.g. srn<epi>
%   'analysis'    prefix after ALL steps (rn<epi>, or srn<epi> if smoothed): images of the first level analysis
%   'mean'        meanP<epi>      mean image of realign 1 / 2
%   'umean'       umeanP<epi>     unwarped mean (fieldmap)
%   'cmean'       cmeanP<epi>     mean coregistered to the T1
%   'wmean'       wmeanP<epi>     mean in template space
%   'wcmean'      wcmeanP<epi>    (intermediate name in sbp_warp_vdm)
%   'ymean'       y_meanP<epi>    deformation field of the nonlinear coregistration
%   'vdm'         vdm5_P<epi>     voxel displacement map
%   'noise_wm_csf' noise_wm_csf_rP<epi>.mat   WM/CSF noise regressors   (sbp_create_noise)
%   'noise_roi'    noise_roi_rP<epi>.mat      ventricle ROI regressors  (sbp_create_noise)
%   'noise_mov'    noise_mov_rP<epi>.mat      movement regressors       (qc/check_movement)
% (the noise files are named after the resliced images they were computed from)
% The fixed part of each prefix (mean, c, w, ...) is the one of the step that creates it.

switch kind
    case 'nordic_in',      f = spm_file(epi,'prefix',epi_prefix(vars,'nordic'));       % input of sbp_nordic
    case 'slicetiming_in', f = spm_file(epi,'prefix',epi_prefix(vars,'slicetiming'));  % input of sbp_slice_timing
    case 'nordic',      f = spm_file(epi,'prefix',['n' epi_prefix(vars,'nordic')]);
    case 'slicetiming', f = spm_file(epi,'prefix',['a' epi_prefix(vars,'slicetiming')]);
    case 'in',          f = spm_file(epi,'prefix',epi_prefix(vars,'realign'));
    case 'mat',         f = spm_file(epi,'prefix',epi_prefix(vars,'realign'),'ext','.mat');
    case 'rp',          f = spm_file(epi,'prefix',['rp_' epi_prefix(vars,'realign')],'ext','.txt');
    case 'res',         f = spm_file(epi,'prefix',['r' epi_prefix(vars,'realign')]);
    case 'smooth_in',   f = spm_file(epi,'prefix',epi_prefix(vars,'smooth'));              % input of sbp_smooth
    case 'smooth',      f = spm_file(epi,'prefix',['s' epi_prefix(vars,'smooth')]);        % sbp_smooth
    case 'analysis',    f = spm_file(epi,'prefix',epi_prefix(vars,'analysis'));            % images of the first level analysis (after all steps)
    case 'mean',        f = spm_file(epi,'prefix',['mean' epi_prefix(vars,'realign')]);
    case 'umean',       f = spm_file(epi,'prefix',['umean' epi_prefix(vars,'realign')]);
    case 'cmean',       f = spm_file(epi,'prefix',['cmean' epi_prefix(vars,'realign')]);
    case 'wmean',       f = spm_file(epi,'prefix',['wmean' epi_prefix(vars,'realign')]);
    case 'wcmean',      f = spm_file(epi,'prefix',['wcmean' epi_prefix(vars,'realign')]);
    case 'ymean',       f = spm_file(epi,'prefix',['y_mean' epi_prefix(vars,'realign')]);
    case 'vdm',         f = spm_file(epi,'prefix',['vdm5_' epi_prefix(vars,'realign')]);
    case 'noise_wm_csf', f = spm_file(epi,'prefix',['noise_wm_csf_r' epi_prefix(vars,'realign')],'ext','.mat'); % sbp_create_noise
    case 'noise_roi',    f = spm_file(epi,'prefix',['noise_roi_r' epi_prefix(vars,'realign')],'ext','.mat');    % sbp_create_noise
    case 'noise_mov',    f = spm_file(epi,'prefix',['noise_mov_r' epi_prefix(vars,'realign')],'ext','.mat');    % qc/check_movement
    otherwise
        error('Unknown kind "%s" (see help epi_file)',kind);
end
end
