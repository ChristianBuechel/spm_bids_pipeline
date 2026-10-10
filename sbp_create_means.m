function sbp_create_means(all_sub_ids)

% function sbp_create_means(all_sub_ids)
% creates the mean of the skull-stripped T1 (wskull-strip-T1.nii), of the GM segment (wc1) and of the mean EPI
% (wmean<P><epi>, P: see epi_prefix) in template space, and their variance images, over the subjects (path.preprocDir)
% Needs sbp_warp for all subjects; at least 2 subjects


% resolve paths
specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

BIDS         = spm_BIDS(path.preprocDir);
n_subs        = length(all_sub_ids);
if n_subs < 2
    error('need at least 2 images to create means')
end
all_wskull_files = [];
all_wmean_files  = [];
all_wc1_files    = [];

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    
    struc_file = get_anat(BIDS,sub_id);
    anat_dir   = spm_file(struc_file,'path');

    epi        = get_ref_epi(BIDS,sub_id);
    func_dir   = spm_file(epi,'path');
    
    wskull_file     = fullfile(anat_dir, 'wskull-strip-T1.nii');
    wmean_file      = epi_file(epi,vars,'wmean'); % wmean<P><epi> (sbp_warp)
    wc1_file        = spm_file(struc_file,'prefix','wc1');
    
    all_wskull_files  = [all_wskull_files;wskull_file];
    all_wmean_files   = [all_wmean_files;wmean_file];
    all_wc1_files     = [all_wc1_files;wc1_file];
    
end
mbr = 1;
matlabbatch{mbr}.spm.util.imcalc.input = all_wskull_files;
matlabbatch{mbr}.spm.util.imcalc.output = 'mean_wskull';
matlabbatch{mbr}.spm.util.imcalc.outdir = cellstr(path.preprocDir);
matlabbatch{mbr}.spm.util.imcalc.expression = 'mean(X)';
matlabbatch{mbr}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{mbr}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{mbr}.spm.util.imcalc.options.mask = 0;
matlabbatch{mbr}.spm.util.imcalc.options.interp = 1;
matlabbatch{mbr}.spm.util.imcalc.options.dtype = 4;
mbr = mbr+1;

matlabbatch{mbr} = matlabbatch{1};
matlabbatch{mbr}.spm.util.imcalc.input = all_wmean_files;
matlabbatch{mbr}.spm.util.imcalc.output = 'mean_wmean';
matlabbatch{mbr}.spm.util.imcalc.expression = 'mean(X)';
mbr = mbr+1;

matlabbatch{mbr} = matlabbatch{1};
matlabbatch{mbr}.spm.util.imcalc.input = all_wc1_files;
matlabbatch{mbr}.spm.util.imcalc.output = 'mean_wc1';
matlabbatch{mbr}.spm.util.imcalc.expression = 'mean(X)';
mbr = mbr+1;

matlabbatch{mbr} = matlabbatch{1};
matlabbatch{mbr}.spm.util.imcalc.input = all_wskull_files;
matlabbatch{mbr}.spm.util.imcalc.output = 'var_wskull';
matlabbatch{mbr}.spm.util.imcalc.expression = 'var(X)';
mbr = mbr+1;

matlabbatch{mbr} = matlabbatch{1};
matlabbatch{mbr}.spm.util.imcalc.input = all_wmean_files;
matlabbatch{mbr}.spm.util.imcalc.output = 'var_wmean';
matlabbatch{mbr}.spm.util.imcalc.expression = 'var(X)';
mbr = mbr+1;

matlabbatch{mbr} = matlabbatch{1};
matlabbatch{mbr}.spm.util.imcalc.input = all_wc1_files;
matlabbatch{mbr}.spm.util.imcalc.output = 'var_wc1';
matlabbatch{mbr}.spm.util.imcalc.expression = 'var(X)';
mbr = mbr+1;

spm_jobman('run',matlabbatch);

end






