function warp_images(all_sub_ids)
% function warp_images
% warps the skullstrip to the template space
% also warps the mean EPI to the template space


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
    epi        = epi(1); % if there are brain and spinal, just take brain
    func_dir   = spm_file(epi,'path');

    strip_file  = fullfile(anat_dir, vars.skullStripID);
    
    c1_file     = spm_file(struc_file,'prefix','c1');    
    c2_file     = spm_file(struc_file,'prefix','c2');
    u_rc1_file  = spm_file(struc_file,'prefix','u_rc1');

    mean_file   = spm_file(epi,'prefix','meana');
    
    flowfield_epi =  fullfile(func_dir,'y_epi_2_template.nii');
    
    gi    = 1;
    matlabbatch{gi,sub}.spm.tools.dartel.crt_warped.flowfields = u_rc1_file;
    matlabbatch{gi,sub}.spm.tools.dartel.crt_warped.images = {strip_file,c1_file,c2_file};
    matlabbatch{gi,sub}.spm.tools.dartel.crt_warped.jactransf = 0;
    matlabbatch{gi,sub}.spm.tools.dartel.crt_warped.K = 6;
    matlabbatch{gi,sub}.spm.tools.dartel.crt_warped.interp = 1;
    gi = gi + 1;
    
    matlabbatch{gi,sub}.spm.util.defs.comp{1}.def = flowfield_epi;
    matlabbatch{gi,sub}.spm.util.defs.out{1}.pull.fnames = mean_file;
    matlabbatch{gi,sub}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
    matlabbatch{gi,sub}.spm.util.defs.out{1}.pull.interp = 4;
    matlabbatch{gi,sub}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{gi,sub}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{gi,sub}.spm.util.defs.out{1}.pull.prefix = 'w';
    gi = gi + 1;
    
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
