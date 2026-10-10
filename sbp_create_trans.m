function sbp_create_trans(all_sub_ids)
% function sbp_create_trans(all_sub_ids)
% creates the deformation fields of the NLIN workflow (nonlinear coregistration + Dartel), written
% to the func folder of the reference EPI (get_ref_epi):
%   y_epi_2_T1.nii         from the nonlinear coregistration (y_mean<P><epi>), defined on the grid of the T1w
%                          (used to get from the T1 to the EPI)
%   y_epi_2_template.nii   the same combined with the Dartel flow field u_rc1<T1w> (EPI --> T1 --> template)
% Needs y_mean<P><epi> (sbp_coreg_nlin_vaso_T1 or sbp_coreg_nlin) and u_rc1<T1w> (sbp_create_dartel).
%   all_sub_ids   subject numbers, e.g. 1:5 (run in parallel if vars.parallel = 1)

% resolve paths
specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

BIDS         = spm_BIDS(path.preprocDir);
n_subs        = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);

    struc_file = get_anat(BIDS,sub_id);
    anat_dir   = spm_file(struc_file,'path');

    epi        = get_ref_epi(BIDS,sub_id);
    func_dir   = spm_file(epi,'path');

    u_rc1_file      = spm_file(struc_file,'prefix','u_rc1');    
    nlin_coreg_file = epi_file(epi,vars,'ymean'); % def field from nlin coreg (y_mean<P><epi>, see epi_prefix)
    %mean_nlin_reg   = epi_file(epi,vars,'cmean');
            
    gi    = 1;
    
    % get deformation from EPI to T1 space
    matlabbatch{gi,sub}.spm.util.defs.comp{1}.def = nlin_coreg_file;
    matlabbatch{gi,sub}.spm.util.defs.comp{2}.id.space = struc_file;  % as this will be used to get from T1 to EPI, we should use the T1 to define the dimensions
    matlabbatch{gi,sub}.spm.util.defs.out{1}.savedef.ofname = 'epi_2_T1';
    matlabbatch{gi,sub}.spm.util.defs.out{1}.savedef.savedir.saveusr = func_dir;
    gi = gi + 1;
    
    % get deformation from EPI to T1 space
    % ie combine nonlin coreg with Dartel (EPI --> T1 --> Template)
    matlabbatch{gi,sub}.spm.util.defs.comp{1}.def = nlin_coreg_file;
    matlabbatch{gi,sub}.spm.util.defs.comp{2}.dartel.flowfield = u_rc1_file;
    matlabbatch{gi,sub}.spm.util.defs.comp{2}.dartel.times = [1 0];
    matlabbatch{gi,sub}.spm.util.defs.comp{2}.dartel.K = 6;
    matlabbatch{gi,sub}.spm.util.defs.comp{2}.dartel.template = {''};
    matlabbatch{gi,sub}.spm.util.defs.out{1}.savedef.ofname = 'epi_2_template';
    matlabbatch{gi,sub}.spm.util.defs.out{1}.savedef.savedir.saveusr = func_dir;
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
