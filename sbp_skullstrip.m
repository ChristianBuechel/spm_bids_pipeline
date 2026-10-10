function sbp_skullstrip(all_sub_ids)
% function sbp_skullstrip(all_sub_ids)
% strips the skull of the T1w image (get_anat): the T1w is multiplied with a brain mask,
% (c1 + c2) > 0.2, i.e. GM + WM of the segmentation of the T1w. The result is written
% to the anat folder as vars.skullStripID (default skull-strip-T1.nii).
%
% Needs c1<T1w> and c2<T1w> in the anat folder, i.e. run after the segmentation of the T1w
% (sbp_coreg_nlin_vaso_T1 or sbp_coreg_nlin, sbp_seg_T1 in the VDM workflow).
%   all_sub_ids   subject numbers, e.g. 1:5 (run in parallel if vars.parallel = 1)

% resolve paths
specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

output_name   = vars.skullStripID;

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    
    struc_file = get_anat(BIDS,sub_id);
    anat_dir   = spm_file(struc_file,'path');
    
    
    c1_file     = spm_file(struc_file,'prefix','c1');
    c2_file     = spm_file(struc_file,'prefix','c2');
    
    gi    = 1;
    
    Vfnames      = [struc_file;c1_file;c2_file];
    matlabbatch{gi,sub}.spm.util.imcalc.input            = Vfnames;
    matlabbatch{gi,sub}.spm.util.imcalc.output           = output_name;
    matlabbatch{gi,sub}.spm.util.imcalc.outdir           = anat_dir;
    matlabbatch{gi,sub}.spm.util.imcalc.expression       = 'i1.*((i2+i3)>0.2)';
    matlabbatch{gi,sub}.spm.util.imcalc.options.dmtx     = 0;
    matlabbatch{gi,sub}.spm.util.imcalc.options.mask     = 0;
    matlabbatch{gi,sub}.spm.util.imcalc.options.interp   = 1;
    matlabbatch{gi,sub}.spm.util.imcalc.options.dtype    = 4;
    
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

