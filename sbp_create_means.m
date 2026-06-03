function create_means(all_sub_ids)

% function create_means(all_sub_ids)
% creates mean skullstrip and mean epi


% resolve paths
[path,vars]  = get_study_specs;
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
    
    struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'type','T1w');
    struc_file = struc_file(1); % if there are more than 1
    anat_dir   = spm_file(struc_file,'path');

    epi        = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    epi        = epi(1); % if there are brain and spinal, just take brain
    func_dir   = spm_file(epi,'path');
    
    wskull_file     = fullfile(anat_dir, 'wskull-strip-T1.nii');
    wmean_file      = spm_file(epi,'prefix','wmeana');
    wc1_file        = fullfile(anat_dir,sprintf('wc1sub-%02d_T1w.nii',sub_id));
    
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






