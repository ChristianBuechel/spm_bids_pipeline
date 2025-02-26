function realign_2_2_vdm(all_sub_ids)
% function realign_2_2_vdm(all_sub_ids)
% realign epi images using a brainmask
% This is part 2 of 2 here we take a brain mask but do not reslice all
% images (this is done when applying the vdms in the next step)

[path,vars]  = get_study_specs; 
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    epi        = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    epi        =  epi(1); % if there are brain and spinal, just take brain
    func_dir   = spm_file(epi,'path');

    cnt        = 1;
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            epi                =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            epi                =  epi(1); % if there are brain and spinal, just take brain
            run_niftis         =  spm_select('ExtFPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^a'),Inf);
            all_niftis{cnt,1}  =  cellstr(run_niftis);
            cnt                =  cnt + 1;
        end
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
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.weight = fullfile(func_dir,'bins3rbrain_mask.nii');
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
