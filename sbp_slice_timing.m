function sbp_slice_timing(all_sub_ids)
% function slice_timing(all_sub_ids)
% corrects slice timing in epi images (output prefix 'a', fixed: a<name>.nii)
% The input are the images of the steps before it (vars.preproc.steps, see epi_prefix),
% e.g. the NORDIC output n<name> if sbp_nordic was run.
% For 3D EPI (vars.sliceTiming.so is empty) slice timing does not make sense: nothing
% is done and nothing is written. Set vars.preproc.steps in get_study_specs so that the
% later steps (realignment, ...) work on the images of the steps that were run, e.g.
% {'nordic'} for 3D EPI with NORDIC, {} for 3D EPI without it.

specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

if isempty(vars.sliceTiming.so) % 3D epifile
    fprintf('3D EPI (vars.sliceTiming.so is empty): no slice timing correction, nothing to do.\n');
    if ismember('slicetiming',epi_steps(vars))
        warning(['vars.preproc.steps lists "slicetiming" (or is not set, which means {''slicetiming''}), but there is no slice timing ' ...
            'for 3D EPI: the later steps would look for images that do not exist. Set vars.preproc.steps in get_study_specs, ' ...
            'e.g. {''nordic''} or {}.']);
    end
    return
end

BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    cnt        = 1;
    all_niftis = {}; % reset: subjects can have different numbers of runs
    for ses = 1:vars.nSess
        for run = get_run_numbers(BIDS,sub_id,ses)
            epi                =  get_func_files(BIDS,sub_id,ses,run);
            all_niftis{cnt,1}  =  cellstr(epi_frames(epi,vars,'slicetiming_in')); %#ok<AGROW> % frames of the input of this step
            cnt                =  cnt + 1;
        end
    end

    gi  = 1;
    matlabbatch{gi,sub}.spm.temporal.st.scans    = all_niftis;
    matlabbatch{gi,sub}.spm.temporal.st.nslices  = vars.sliceTiming.nslices;
    matlabbatch{gi,sub}.spm.temporal.st.tr       = vars.sliceTiming.tr;
    matlabbatch{gi,sub}.spm.temporal.st.ta       = 0;
    matlabbatch{gi,sub}.spm.temporal.st.so       = vars.sliceTiming.so;
    matlabbatch{gi,sub}.spm.temporal.st.refslice = vars.sliceTiming.refslice;
    matlabbatch{gi,sub}.spm.temporal.st.prefix   = 'a'; % fixed output prefix of slice timing
end

% run matlabbatch
n_procs = n_subs; % to not block too many cores on the server

if n_procs > vars.max_procs
    n_procs = vars.max_procs;
end

if vars.parallel == 1
    run_spm_parallel(matlabbatch, n_procs);
    %run_spm_multiple(matlabbatch, n_procs);
else
    run_spm_sequential(matlabbatch);
end
