function apply_vdm(all_sub_ids)
% function realign(all_sub_ids)
% realign and reslice epi images
% This is part 2 of 2 here we take a brain mask and actually reslice all images

[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);

    gi  = 1;
    cnt = 1;
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            epi                =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            epi                =  epi(1); % if there are brain and spinal, just take brain
            run_niftis         =  spm_select('ExtFPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^a'),Inf);
            matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.data(cnt).scans   = cellstr(run_niftis);
            vdm_file = spm_file(epi,'prefix','vdm5_a');
            matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.data(cnt).vdmfile = cellstr(vdm_file);
            cnt                =  cnt + 1;
        end
    end
    
    matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.roptions.pedir = 2; % Y direction
    matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.roptions.which = [2 1];
    matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.roptions.rinterp = 4;
    matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.roptions.wrap = [0 0 0];
    matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.roptions.mask = 1;
    matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.roptions.prefix = 'r';
end

% save('apply_vdm','matlabbatch');
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
