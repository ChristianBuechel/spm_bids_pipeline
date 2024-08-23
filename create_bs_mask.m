function create_bs_mask(all_sub_ids)
% function create_bs_mask
% warps a brainstem mask defined in template space to the EPI space
% smooths it with 3mm
% and make it a binary mask using imcalc

% resolve paths
[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);

n_subs        = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    
    epi        = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    epi        = epi(1); % if there are brain and spinal, just take brain
    func_dir   = spm_file(epi,'path');
    
    struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'type','T1w');
    struc_file = struc_file(1); % if there are more than 1
    anat_dir   = spm_file(struc_file,'path');
    c3_file    = spm_file(struc_file,'prefix','c3');
    u_rc1_file = spm_file(struc_file,'prefix','u_rc1');

    
    bs_file    = fullfile(path.templateDir,vars.brainstemID);
    
    flowfield_epi =  fullfile(func_dir,'y_epi_2_template.nii');
    flowfield_T1  =  fullfile(func_dir,'y_epi_2_T1.nii');
    
    % warp brainstem mask (template space) to EPI space
    gi    = 1;
    matlabbatch{gi,sub}.spm.util.defs.comp{1}.def = flowfield_epi;
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fnames = {bs_file};
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.weight = {''};
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.savedir.saveusr = func_dir;
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fov.bbvox.bb = [NaN NaN NaN;...
        NaN NaN NaN];
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fov.bbvox.vox = [NaN NaN NaN];
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.preserve = 0;
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.prefix = 'c';
    gi = gi + 1;
    
    % smooth it with 3mm
    matlabbatch{gi,sub}.spm.spatial.smooth.data(1,1) = {spm_file(bs_file,'prefix','c','path',func_dir)};
    matlabbatch{gi,sub}.spm.spatial.smooth.fwhm   = [3 3 3];
    matlabbatch{gi,sub}.spm.spatial.smooth.dtype  = 0;
    matlabbatch{gi,sub}.spm.spatial.smooth.im     = 0;
    matlabbatch{gi,sub}.spm.spatial.smooth.prefix = 's3';
    gi = gi + 1;
    
    % and binarize it
    matlabbatch{gi,sub}.spm.util.imcalc.input            = {spm_file(bs_file,'prefix','s3c','path',func_dir)};
    matlabbatch{gi,sub}.spm.util.imcalc.output           = spm_file(spm_file(bs_file,'prefix','bins3c'),'filename');
    matlabbatch{gi,sub}.spm.util.imcalc.outdir           = func_dir;
    matlabbatch{gi,sub}.spm.util.imcalc.expression       = 'i1>0';
    matlabbatch{gi,sub}.spm.util.imcalc.options.dmtx     = 0;
    matlabbatch{gi,sub}.spm.util.imcalc.options.mask     = 0;
    matlabbatch{gi,sub}.spm.util.imcalc.options.interp   = 1;
    matlabbatch{gi,sub}.spm.util.imcalc.options.dtype    = 2; %uint 8
    gi = gi + 1;
    
    % now warp brainstem mask from template to T1 space
    matlabbatch{gi,sub}.spm.util.defs.comp{1}.dartel.flowfield = u_rc1_file;
    matlabbatch{gi,sub}.spm.util.defs.comp{1}.dartel.times = [1 0];
    matlabbatch{gi,sub}.spm.util.defs.comp{1}.dartel.K = 6;
    matlabbatch{gi,sub}.spm.util.defs.comp{1}.dartel.template = {''};
    
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fnames = {bs_file};
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.weight = {''};
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.savedir.saveusr = anat_dir;
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fov.bbvox.bb = [NaN NaN NaN;...
        NaN NaN NaN];
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fov.bbvox.vox = [NaN NaN NaN];
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.preserve = 2;
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.prefix = 'c';
    gi = gi + 1;
    
    % and mask the c3 (CSF) image with it
    matlabbatch{gi,sub}.spm.util.imcalc.input            = [c3_file; spm_file(bs_file,'prefix','c','path',anat_dir)];
    matlabbatch{gi,sub}.spm.util.imcalc.output           = spm_file(spm_file(bs_file,'prefix','c3'),'filename');
    matlabbatch{gi,sub}.spm.util.imcalc.outdir           = anat_dir;
    matlabbatch{gi,sub}.spm.util.imcalc.expression       = 'i1.*i2';
    matlabbatch{gi,sub}.spm.util.imcalc.options.dmtx     = 0;
    matlabbatch{gi,sub}.spm.util.imcalc.options.mask     = 0;
    matlabbatch{gi,sub}.spm.util.imcalc.options.interp   = 1;
    matlabbatch{gi,sub}.spm.util.imcalc.options.dtype    = 2; %uint 8
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






