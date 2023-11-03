function skullstrip(all_sub_ids)
% function skullstrip
% strips the skull

% resolve paths
[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);
run_parallel = 1;

output_name   = vars.skullStripID; 

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    
    struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'type','T1w');
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

if run_parallel == 1
    run_spm_parallel(matlabbatch, n_procs);
else
    spm_jobman('run',matlabbatch);
end
end






