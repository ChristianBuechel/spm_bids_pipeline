function check_images(all_sub_ids)

% function check_images(all_sub_ids);
% displays resulting images in coreg to see whether everything is OK


% resolve paths
[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs        = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    
    
    struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'type','T1w');
    struc_file = struc_file(1); % if there are more than 1
    anat_dir   = spm_file(struc_file,'path');

    epi        = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    func_dir   = spm_file(epi,'path');
    
    skull_file      = fullfile(anat_dir, vars.skullStripID);
    bm_file         = fullfile(anat_dir, vars.T1maskID); % brainmask in T1 space
    mean_file       = spm_file(epi,'prefix','meana');
    rbm_file        = fullfile(func_dir, spm_file(vars.T1maskID,'prefix','r')); %brainmask in EPI space
    wm_file         = spm_file(epi,'prefix','wmeana');

    template_file   = fullfile(path.templateDir,vars.templateT1ID);
    
    matlabbatch{1}.spm.util.checkreg.data = [skull_file;bm_file;mean_file;rbm_file;wm_file;template_file];
    % Images in the same row are in the same space !
    
    spm_jobman('run',matlabbatch);
    fprintf('Checking Subject %d \nPress enter in command window to continue\n',sub_id);
    input('');
end

end






