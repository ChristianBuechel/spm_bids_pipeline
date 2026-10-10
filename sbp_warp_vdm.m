function warp_images_vdm(all_sub_ids)
% function warp_images
% warps the skullstrip to the template space
% also warps the mean EPI (cmean<P><epi>, P: see epi_prefix) to the template space. The result
% (wcmean<P><epi>) is renamed to wmean<P><epi>, the name the nlin workflow produces too


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

    mean_file   = epi_file(epi,vars,'cmean'); % cmean<P><epi>: for vdm pipeline we need the vdm corrected and lin coreg mean as the y_epi_2_T1 is a 1:1 mapping
    mean_dir    = spm_file(mean_file,'path');
    
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
    
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = epi_file(epi,vars,'wcmean'); % it is now called wcmean<P><epi> (w + cmean<P><epi>) ...
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = mean_dir;
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(1).pattern = '^.*$'; % ... the whole name is replaced by
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(1).repl = char(spm_file(epi_file(epi,vars,'wmean'),'basename')); % wmean<P><epi>, the name the later steps (and the nlin workflow) use
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
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
