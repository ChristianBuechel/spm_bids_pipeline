function sbp_warp(all_sub_ids)
% function sbp_warp(all_sub_ids)
% warps the skull-stripped T1 and its GM / WM segments to the template space with the Dartel
% flow field u_rc1<T1w> (sbp_create_dartel): w<vars.skullStripID>, wc1<T1w>, wc2<T1w> in the anat folder.
% Also warps the mean EPI of realign 2 (mean<P><epi>, P: see epi_prefix) to the template space with
% y_epi_2_template.nii (sbp_create_trans): wmean<P><epi> in the func folder.
% Needs sbp_skullstrip, sbp_create_dartel, sbp_create_trans and sbp_realign_2_2.
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

    strip_file  = fullfile(anat_dir, vars.skullStripID);
    
    c1_file     = spm_file(struc_file,'prefix','c1');    
    c2_file     = spm_file(struc_file,'prefix','c2');
    u_rc1_file  = spm_file(struc_file,'prefix','u_rc1');

    mean_file   = epi_file(epi,vars,'mean'); % mean<P><epi> (of realign 2), is written as wmean<P><epi> (prefix w, see epi_prefix)

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
