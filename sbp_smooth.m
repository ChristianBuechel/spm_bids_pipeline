function sbp_smooth(all_sub_ids)
% function sbp_smooth(all_sub_ids)
% smooths the 4D EPI images of all runs (of vars.task, brain), e.g. a few mm before a first level
% analysis in native space. The smoothed images get the fixed prefix 's': s<in>.nii, where <in> is the
% image of the step before smoothing (epi_file(epi,vars,'smooth_in'), normally the resliced images
% r<P><epi>.nii of sbp_realign_2_2 / sbp_apply_vdm, e.g. rn<epi>.nii --> srn<epi>.nii).
% The step has to be listed in vars.preproc.steps AFTER 'realign', e.g.
%   vars.preproc.steps = {'nordic','realign','smooth'};     (see epi_steps, epi_prefix)
% The analysis (sbp_analyses) then uses the smoothed images (epi_file kind 'analysis'). The noise regressors
% (sbp_create_noise, WM / CSF components) and the realignment parameters are NOT affected: they are
% computed from the unsmoothed resliced images r<P><epi>.nii.
%
% Settings (get_study_specs):
%   vars.smooth.fwhm   FWHM of the Gaussian kernel in mm, scalar or 1x3 (default 1.5)
%   all_sub_ids   subject numbers, e.g. 1:5 (run in parallel if vars.parallel = 1)

specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

if ~ismember('smooth',epi_steps(vars))
    error(['sbp_smooth is run, but vars.preproc.steps in get_study_specs does not list "smooth". Add it after "realign", ' ...
        'e.g. {''nordic'',''realign'',''smooth''}, otherwise the analysis will not use the smoothed images (s<name>.nii).']);
end

fwhm = 1.5;
if isfield(vars,'smooth') && isfield(vars.smooth,'fwhm'), fwhm = vars.smooth.fwhm; end
if isscalar(fwhm), fwhm = repmat(fwhm,1,3); end

BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id = all_sub_ids(sub);
    gi     = 1;
    for ses = 1:vars.nSess
        for run = get_run_numbers(BIDS,sub_id,ses)
            epi    = get_func_files(BIDS,sub_id,ses,run);
            frames = cellstr(epi_frames(epi,vars,'smooth_in')); % all frames of the input image (error if it does not exist)
            fprintf('sub-%02d ses-%02d run-%02d: %s --> %s (%.1f mm)\n',sub_id,ses,run,spm_file(frames{1},'filename'),['s' spm_file(frames{1},'filename')],fwhm(1));
            matlabbatch{gi,sub}.spm.spatial.smooth.data   = frames;
            matlabbatch{gi,sub}.spm.spatial.smooth.fwhm   = fwhm;
            matlabbatch{gi,sub}.spm.spatial.smooth.dtype  = 0; % as the input
            matlabbatch{gi,sub}.spm.spatial.smooth.im     = 0;
            matlabbatch{gi,sub}.spm.spatial.smooth.prefix = 's'; % fixed, see epi_prefix
            gi = gi + 1;
        end
    end
end

% run matlabbatch
n_procs = n_subs; % to not block too many cores on the server

if n_procs > vars.max_procs
    n_procs = vars.max_procs;
end

if vars.parallel == 1
    run_spm_parallel(matlabbatch, n_procs);
else
    run_spm_sequential(matlabbatch);
end
