function sbp_create_mask(all_sub_ids)
% function sbp_create_mask(all_sub_ids)
% creates a brain mask based on GM and WM of the T1w (c1 + c2, vars.T1maskID, e.g. brain_mask.nii,
% in the anat folder). This mask is then pushed into EPI space with the deformation of the
% nonlinear coregistration (y_epi_2_T1.nii, sbp_create_trans) and written to the func folder
% as r<mask> (rbrain_mask.nii), smoothed (s3r<mask>) and binarised (bins3r<mask>), to be used
% as a mask for the 1st level GLMs. The skull-stripped T1 (vars.skullStripID) and the plain
% T1w are also pushed into EPI space (r<name>, anat folder), for checking.
% Needs c1<T1w>, c2<T1w> (T1 segmentation), the skull-stripped T1 (sbp_skullstrip) and y_epi_2_T1.nii.
%   all_sub_ids   subject numbers, e.g. 1:5 (run in parallel if vars.parallel = 1)

% add paths
specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

mask_name   = vars.T1maskID;
r_pref      = 'r';
r_mask_name = [r_pref mask_name];
skern       = 3; % smoothing the final mask in epi space


for sub = 1:n_subs
    sub_id    = all_sub_ids(sub);
    
    struc_file = get_anat(BIDS,sub_id);
    anat_dir   = spm_file(struc_file,'path');

    epi        = get_ref_epi(BIDS,sub_id);
    func_dir   = spm_file(epi,'path');
  
    c1_file     = spm_file(struc_file,'prefix','c1');    
    c2_file     = spm_file(struc_file,'prefix','c2');
    
    gi    = 1;
    matlabbatch{gi,sub}.spm.util.imcalc.input          = [c1_file; c2_file]; 
    matlabbatch{gi,sub}.spm.util.imcalc.output         = mask_name;
    matlabbatch{gi,sub}.spm.util.imcalc.outdir         = anat_dir;
    matlabbatch{gi,sub}.spm.util.imcalc.expression     = 'i1 + i2'; % simply add WM and GM
    matlabbatch{gi,sub}.spm.util.imcalc.options.dmtx   = 0;
    matlabbatch{gi,sub}.spm.util.imcalc.options.mask   = 0;
    matlabbatch{gi,sub}.spm.util.imcalc.options.interp = 1;
    matlabbatch{gi,sub}.spm.util.imcalc.options.dtype  = 2; % uint8 is sufficient for a binary mask
    gi    = gi + 1;

    matlabbatch{gi,sub}.spm.util.defs.comp{1}.def = fullfile(func_dir, 'y_epi_2_T1.nii');
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fnames = fullfile(anat_dir, mask_name);
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.weight = {''};
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.savedir.saveusr = func_dir;
    %matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fov.bbvox.bb = [-100 -150 -100;100 150 100]; % same as nlin_coreg
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fov.bbvox.bb = [NaN NaN NaN;NaN NaN NaN]; % same as nlin_coreg
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fov.bbvox.vox = [NaN NaN NaN];
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.preserve = 2; %option for labels, does a perfect job for binary image
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.prefix = r_pref;
    gi = gi + 1;
    
    matlabbatch{gi,sub}.spm.util.defs.comp{1}.def = fullfile(func_dir, 'y_epi_2_T1.nii');
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fnames = {char(fullfile(anat_dir, vars.skullStripID));char(struc_file)}; %also create a skullstrip and plain T1w in EPI space
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.weight = {''};
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.savedir.saveusr = anat_dir;
    %matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fov.bbvox.bb = [-100 -150 -100;100 150 100]; % same as nlin coreg
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fov.bbvox.bb = [NaN NaN NaN;NaN NaN NaN]; % same as nlin coreg
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fov.bbvox.vox = [1 1 1];
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.preserve = 0; % simple
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    matlabbatch{gi,sub}.spm.util.defs.out{1}.push.prefix = r_pref;
    gi = gi + 1;

    matlabbatch{gi,sub}.spm.spatial.smooth.data   = fullfile(func_dir, r_mask_name);
    matlabbatch{gi,sub}.spm.spatial.smooth.fwhm   = repmat(skern,1,3);
    matlabbatch{gi,sub}.spm.spatial.smooth.prefix = ['s' num2str(skern)];
    gi = gi + 1;
    
    matlabbatch{gi,sub}.spm.util.imcalc.input            = {spm_file(r_mask_name,'prefix','s3','path',func_dir)};
    matlabbatch{gi,sub}.spm.util.imcalc.output           = spm_file(spm_file(r_mask_name,'prefix','bins3'),'filename');
    matlabbatch{gi,sub}.spm.util.imcalc.outdir           = func_dir;
    matlabbatch{gi,sub}.spm.util.imcalc.expression       = 'i1>0';
    matlabbatch{gi,sub}.spm.util.imcalc.options.dmtx     = 0;
    matlabbatch{gi,sub}.spm.util.imcalc.options.mask     = 0;
    matlabbatch{gi,sub}.spm.util.imcalc.options.interp   = 1;
    matlabbatch{gi,sub}.spm.util.imcalc.options.dtype    = 2; %uint 8

     
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