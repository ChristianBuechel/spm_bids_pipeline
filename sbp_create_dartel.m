function sbp_create_dartel(all_sub_ids)
% function sbp_create_dartel(all_sub_ids)
% estimates the Dartel flow field from the T1w (GM and WM, rc1 and rc2) to the template
% (vars.templateID in path.templateDir, 6 templates). Writes u_rc1<T1w>.nii in the anat folder.
%
% Needs rc1<T1w> and rc2<T1w> (Dartel imported segments of the T1w), i.e. run after the
% segmentation (sbp_coreg_nlin_vaso_T1 or sbp_coreg_nlin, sbp_seg_T1 in the VDM workflow).
% The flow field is used by sbp_create_trans.
%   all_sub_ids   subject numbers, e.g. 1:5 (run in parallel if vars.parallel = 1)


% resolve paths
specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

BIDS         = spm_BIDS(path.preprocDir);

template_t   = vars.templateID;

n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    
    struc_file = get_anat(BIDS,sub_id);
    anat_dir   = spm_file(struc_file,'path');

    rc1_file     = spm_file(struc_file,'prefix','rc1');    
    rc2_file     = spm_file(struc_file,'prefix','rc2');    

    gi    = 1;
    %matlabbatch{gi,sub}.spm.tools.dartel.warp1.images = {cellstr(rc1_file{1}),cellstr(rc2_file{1})};
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.images = {rc1_file rc2_file}';
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.rform = 0;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(1).its = 3;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(1).K = 0;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(1).template = {fullfile(path.templateDir,sprintf(template_t,1))};
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(2).its = 3;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(2).K = 0;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(2).template = {fullfile(path.templateDir,sprintf(template_t,2))};
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(3).its = 3;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(3).K = 1;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(3).template = {fullfile(path.templateDir,sprintf(template_t,3))};
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(4).its = 3;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(4).K = 2;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(4).template = {fullfile(path.templateDir,sprintf(template_t,4))};
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(5).its = 3;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(5).K = 4;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(5).template = {fullfile(path.templateDir,sprintf(template_t,5))};
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(6).its = 3;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(6).K = 6;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(6).template = {fullfile(path.templateDir,sprintf(template_t,6))};
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.optim.its = 3;
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
