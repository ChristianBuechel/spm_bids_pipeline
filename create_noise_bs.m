function create_noise_bs(all_sub_ids)
% function create_noise
% creates noise regressors for 1st ´level analysis
% 1) CSF for brainstem

% add paths
[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);
run_parallel = 1;

for sub = 1:n_subs
    sub_id    = all_sub_ids(sub);
    
    struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'type','T1w');
    struc_file = struc_file(1); % if there are more than 1
    anat_dir   = spm_file(struc_file,'path');
    
    epi        = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    func_dir   = spm_file(epi,'path');
    
    bs_file    = fullfile(path.templateDir,'brainstem_mask.nii');
    c_files{1} = spm_file(bs_file,'prefix','c3','path',anat_dir);
   
    cm_file        = spm_file(epi,'prefix','cmeana');
    epi_2_T1       = fullfile(func_dir, 'y_epi_2_T1.nii');
    
    
    
    
    
    thresh = [0.1]; ex_var = 6; %low threshold as erosion takes away much of the bs
    gi = 1;
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            fprintf('ses %d run %d\n',ses,run);
            epi       =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            epi_files =  spm_select('ExtFPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^bra'),Inf);
            %extract_segment(epi_files,cm_file,epi_2_T1,c_files,[],thresh,ex_var);
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.evaluated = epi_files;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.string    = char(cm_file);
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{3}.string    = char(epi_2_T1);
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{4}.evaluated = c_files;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{5}.evaluated = [];
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{6}.evaluated = thresh;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{7}.evaluated = ex_var;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{8}.string    = 'noise_csf_';
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'extract_segment';
            gi    = gi + 1;
        end
    end
    
end

% run matlabbatch
n_procs = n_subs; % to not block to many cores on the server

if n_procs > vars.max_procs
    n_procs = vars.max_procs;
end

if run_parallel == 1
    run_spm_parallel(matlabbatch, n_procs);
    %run_spm_multiple(matlabbatch,n_procs)
else
    spm_jobman('run',matlabbatch);
end

end