function sbp_nordic(all_sub_ids)
% function sbp_nordic
% does NORDIC noise filtering on data in image space (magnitude + phase of the main
% task, brain). The denoised magnitude and phase images are written next to the input
% with the prefix 'n': n<name>.nii (the prefix is prepended to the name of the input,
% see epi_prefix: vars.preproc.steps says what was run before NORDIC, normally nothing).
% The vasoref runs are handled by sbp_vaso_T1 (run it BEFORE sbp_nordic: it provides the noise level
% from the noise volumes of the vasoref runs, which is used here if the runs have no noise volumes).

% add paths
specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

if ~ismember('nordic',epi_steps(vars))
    error(['sbp_nordic is run, but vars.preproc.steps in get_study_specs does not list "nordic". Add it (e.g. {''nordic''} or ' ...
        '{''nordic'',''slicetiming''}), otherwise realignment and the later steps will not use the denoised images (n<name>.nii).']);
end

types = {'bold','cbv'};

% noise level for NORDIC (see NIFTI_NORDIC_fast): runs with noise volumes at their end use them
% (vars.nordic.func_noisevols: number or 'auto' = detected, default; split off into noise_<name>.nii);
% all other runs use the noise level of the vasoref run of the same session, if sbp_vaso_T1 was run
% before, otherwise NORDIC assumes unit noise.
n_noise = 'auto';
if isfield(vars,'nordic') && isfield(vars.nordic,'func_noisevols'), n_noise = vars.nordic.func_noisevols; end

arg = [];
arg.temporal_phase     = 1;
arg.phase_filter_width = 10;
arg.make_complex_nii   = 0; % only the denoised magnitude is written (the phase is input, but not used for the stats: saves space)

for sub = 1:n_subs
    sub_id    = all_sub_ids(sub);

    gi = 1;
    for ses = 1:vars.nSess
        for run = get_run_numbers(BIDS,sub_id,ses)
            for itype = 1:numel(types)
                fprintf('ses %d run %d\n',ses,run);
                f_magn = get_func_files(BIDS,sub_id,ses,run,'part','mag','type',types{itype},'expect',NaN); % brain only; empty if this run is of the other type
                if ~isempty(f_magn)
                    f_phase = get_func_files(BIDS,sub_id,ses,run,'part','phase','type',types{itype},'expect',NaN);
                    if ~isempty(f_phase)
                        in_magn  = char(epi_file(f_magn, vars,'nordic_in')); % raw images, unless earlier steps were run
                        in_phase = char(epi_file(f_phase,vars,'nordic_in'));
                        if ~isfile(in_magn) || ~isfile(in_phase)
                            error('sub-%02d ses-%02d run-%02d: NORDIC input not found (%s). Check vars.preproc.steps.',sub_id,ses,run,in_magn);
                        end
                        matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.string    = in_magn;
                        matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.string    = in_phase;
                        matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{3}.string    =  'n'; % output prefix of NORDIC (fixed), see epi_prefix
                        arg_run = arg;
                        arg_run.noise_volume_split = n_noise;   % own noise volumes (if any), they win
                        arg_run.noise_level_file   = get_noise_file(BIDS,sub_id,ses,vars);  % else the noise level of the vasoref run ('' = none)
                        matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{4}.evaluated = arg_run;
                        matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
                        matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'NIFTI_NORDIC_fast';
                        gi    = gi + 1;
                    end
                end
            end
        end
    end
end

% run matlabbatch
n_procs = n_subs; % to not block to many cores on the server

if n_procs > vars.max_procs
    n_procs = vars.max_procs;
end

if vars.parallel == 1
    run_spm_parallel(matlabbatch, n_procs);
    %run_spm_multiple(matlabbatch, n_procs);
else
    run_spm_sequential(matlabbatch);
end
