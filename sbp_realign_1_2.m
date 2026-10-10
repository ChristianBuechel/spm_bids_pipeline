function realign_1_2(all_sub_ids,runs)
% function realign_1_2(all_sub_ids)
% realign epi images (first pass)
% This is part 1 of 2 here we simply create a mean EPI for all the upcoming
% spatial things, but also update the .mat info
% Part 2 will then use a brain mask to get rid of the eyes
%
% The images realigned are the output of the steps before (vars.preproc.steps, see
% epi_prefix: P<name>.nii, e.g. an<name>.nii for NORDIC + slice timing). All runs and
% sessions of a subject are realigned together to the first image. The mean image
% (mean<P><name>.nii, of ses-01 run-01) and the realignment parameters (rp_<P><name>.txt)
% get this prefix too, all later steps find them with epi_file.

specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;
if nargin < 2
    runs = []; % all existing runs
end
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    cnt        = 1;
    all_niftis = {}; % reset: subjects can have different numbers of runs
    for ses = 1:vars.nSess
        if isempty(runs)
            ses_runs = get_run_numbers(BIDS,sub_id,ses);
        else
            ses_runs = runs;
        end
        for run = ses_runs
            epi               =  get_func_files(BIDS,sub_id,ses,run);
            if cnt == 1 && ~strcmp(char(epi),char(get_ref_epi(BIDS,sub_id)))
                % SPM names the mean after the FIRST image, the later steps expect it for ses-01 run-01
                error('sub-%02d: the first run to realign (ses-%02d run-%02d) is not ses-01 run-01, but the later steps look for the mean image of ses-01 run-01.',sub_id,ses,run);
            end
            all_niftis{cnt,1}  =  cellstr(epi_frames(epi,vars,'in')); %#ok<AGROW> % frames of P<name>
            cnt                =  cnt + 1;
        end
    end
    % fill with default options
    gi  = 1;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.data             = all_niftis;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.quality = 0.95; % was 0.9 in SPM12
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.sep     = 1.5; % was 4 in SPM12
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.fwhm    = 1; % was 5 in SPM12;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.rtm     = 0; % one pass is ok
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.interp  = 2;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.weight  = '';
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.which   = [0 1]; %only mean
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.interp  = 4;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.mask    = 1; % if set to 1 takes out voxels from people who move out of FOV
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.prefix  = 'r';    
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
