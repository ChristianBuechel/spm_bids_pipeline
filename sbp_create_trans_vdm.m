function create_trans_vdm(all_sub_ids)
% function create_trans
% creates all required spatial transformation matrices based on nonlinear coreg and dartel
% when we use vdms (fieldmaps) some deformation fields are simply 1:1
% mappings, but we keep to concept for simplicity

% resolve paths
specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);

    struc_file = get_anat(BIDS,sub_id);

    epi        = get_ref_epi(BIDS,sub_id);
    func_dir   = spm_file(epi,'path');

    u_rc1_file      = spm_file(struc_file,'prefix','u_rc1');    
    %mean_reg        = epi_file(epi,vars,'mean');
    mean_reg        = epi_file(epi,vars,'cmean'); % what does this do? (not used, see comments below)
     
    %needs to be bigger?
    
    gi    = 1;
    % get deformation from EPI to T1 space (1:1)
    %matlabbatch{gi,sub}.spm.util.defs.comp{1}.id.space = mean_reg; %
    matlabbatch{gi,sub}.spm.util.defs.comp{1}.id.space = struc_file; % T1w as we use it to get from T1 space to EPI space
    matlabbatch{gi,sub}.spm.util.defs.out{1}.savedef.ofname = 'epi_2_T1';
    matlabbatch{gi,sub}.spm.util.defs.out{1}.savedef.savedir.saveusr = func_dir;
    gi = gi + 1;
    
    % get deformation from EPI to T1 space
    % ie combine nonlin coreg with Dartel (EPI --> T1 --> Template)
    matlabbatch{gi,sub}.spm.util.defs.comp{1}.dartel.flowfield = u_rc1_file;
    matlabbatch{gi,sub}.spm.util.defs.comp{1}.dartel.times = [1 0];
    matlabbatch{gi,sub}.spm.util.defs.comp{1}.dartel.K = 6;
    matlabbatch{gi,sub}.spm.util.defs.comp{1}.dartel.template = {''};
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
