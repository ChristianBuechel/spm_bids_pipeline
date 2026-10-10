function f = get_noise_file(BIDS,sub,ses,vars)

% function f = get_noise_file(BIDS,sub,ses,vars)
% name of the noise level file (noise_<name>.json, written by NIFTI_NORDIC_fast when the noise
% volumes of a series are split off) of the reference run (vars.ref_task, default 'vasoref') of
% this session. NIFTI_NORDIC_fast uses it (ARG.noise_level_file) for series without own noise volumes.
% The file is only written when sbp_vaso_T1 (NORDIC of the reference runs) was run, so this returns
% the name, not a guarantee that the file exists. Returns '' if the session has no reference run
% or vars.nordic.ref_noise is 0.

f = '';
if isfield(vars,'nordic') && isfield(vars.nordic,'ref_noise') && ~vars.nordic.ref_noise
    return
end
if isfield(vars,'ref_task'), ref_task = vars.ref_task; else, ref_task = 'vasoref'; end
mag = get_func_files(BIDS,sub,ses,[],'task',ref_task,'type','bold','part','mag','expect',NaN);
if isempty(mag)
    mag = get_func_files(BIDS,sub,ses,[],'task',ref_task,'type','cbv','part','mag','expect',NaN);
end
if isempty(mag), return, end
f = spm_file(mag{1},'prefix','noise_','ext','json');   % in the folder of the raw image
end
