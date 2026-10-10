function realign_2_2_vdm(all_sub_ids)
% function realign_2_2_vdm(all_sub_ids)
% realign epi images using a brainmask
% This is part 2 of 2 here we take a brain mask but do not reslice all
% images (this is done when applying the vdms in the next step)
% The images are the P<epi> images (P: see epi_prefix); the realignment parameters
% rp_<P><epi>.txt and the header sidecars P<epi>.mat are updated. The brain mask is the
% one of sbp_create_mask (bins3r<vars.T1maskID> in the func folder).

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
            all_niftis{cnt,1}  =  cellstr(epi_frames(epi,vars,'in')); %#ok<AGROW> % frames of P<epi>
            cnt                =  cnt + 1;
        end
    end
    weight_mask = fullfile(func_dir,spm_file(['r' vars.T1maskID],'prefix','bins3')); % made by sbp_create_mask
    if ~isfile(char(weight_mask))
        error('sub-%02d: brain mask %s not found, was sbp_create_mask run?',sub_id,char(weight_mask));
    end
    % fill with default options
    gi  = 1;    
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.data = all_niftis;
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.quality = 0.95;
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.sep = 1.5;
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.fwhm = 1;
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.rtm = 1;
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.interp = 2;
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.weight = cellstr(weight_mask);
end

%save('realign_2_2_vdm','matlabbatch');
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
