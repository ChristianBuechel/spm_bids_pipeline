function sbp_realign_2_2(all_sub_ids)
% function sbp_realign_2_2(all_sub_ids)
% This is part 2 of 2 of realignement, where we take a brain mask and actually reslice all images
% The images realigned are the same as in part 1 (P<name>.nii, P: see epi_prefix, e.g. 'an' for
% NORDIC + slice timing). The resliced images get the fixed prefix 'r' (r<P><name>.nii), the mean
% image (mean<P><name>.nii) and the realignment parameters (rp_<P><name>.txt) of part 1 are
% replaced by the ones of this run. The brain mask for the weighting is the one of sbp_create_mask
% (bins3r<vars.T1maskID> in the func folder).

specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    epi        = get_ref_epi(BIDS,sub_id);
    func_dir   = spm_file(epi,'path');

    cnt        = 1;
    all_niftis = {}; % reset: subjects can have different numbers of runs
    for ses = 1:vars.nSess
        for run = get_run_numbers(BIDS,sub_id,ses)
            epi                =  get_func_files(BIDS,sub_id,ses,run);
            if cnt == 1 && ~strcmp(char(epi),char(get_ref_epi(BIDS,sub_id)))
                % SPM names the mean after the FIRST image, the later steps expect it for ses-01 run-01
                error('sub-%02d: the first run to realign (ses-%02d run-%02d) is not ses-01 run-01, but the later steps look for the mean image of ses-01 run-01.',sub_id,ses,run);
            end
            all_niftis{cnt,1}  =  cellstr(epi_frames(epi,vars,'in')); %#ok<AGROW> % frames of P<name>
            cnt                =  cnt + 1;
        end
    end
    weight_mask = fullfile(func_dir,spm_file(['r' vars.T1maskID],'prefix','bins3')); % made by sbp_create_mask
    if ~isfile(char(weight_mask))
        error('sub-%02d: brain mask %s not found, was sbp_create_mask run?',sub_id,char(weight_mask));
    end
    % fill with default options
    gi  = 1;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.data             = all_niftis;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.quality = 0.95;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.sep     = 1.5;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.fwhm    = 1;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.rtm     = 1; % two pass now
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.interp  = 2;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.weight  = cellstr(weight_mask);
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.which   = [2 1]; % all images and mean
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.interp  = 4;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.mask    = 1; % if set to 1 takes out people who move out of FOV
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
