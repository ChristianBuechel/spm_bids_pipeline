function create_noise(all_sub_ids)
% function create_noise
% creates noise regressors for 1st ´level analysis
% 1) WM and CSF
% 2) ROI at the posterior tip of the lateral ventricles (as in Horing et al 202X)

% add paths
[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id    = all_sub_ids(sub);
    
    struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'type','T1w');
    struc_file = struc_file(1); % if there are more than 1
    anat_dir   = spm_file(struc_file,'path');
    
    epi        = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    epi        = epi(1); % if there are brain and spinal, just take brain

    func_dir   = spm_file(epi,'path');
    
    c_files{1}  = spm_file(struc_file,'prefix','c2');
    c_files{2}  = spm_file(struc_file,'prefix','c3');
    
    wm_file        = spm_file(epi,'prefix','wmeana');
    cm_file        = spm_file(epi,'prefix','cmeana');
    epi_2_T1       = fullfile(func_dir, 'y_epi_2_T1.nii');
    epi_2_templ    = fullfile(func_dir, 'y_epi_2_template.nii');
    
    
    thresh = [0.9 0.9]; ex_var = 6;
    gi = 1;
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            fprintf('ses %d run %d\n',ses,run);
            epi       =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            epi       =  epi(1); % if there are brain and spinal, just take brain
            epi_files =  spm_select('ExtFPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^ra'),Inf);
            %extract_segment(epi_files,cm_file,epi_2_T1,c_files,[],thresh,ex_var);
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.evaluated = epi_files;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.string    = char(cm_file);
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{3}.string    = char(epi_2_T1);
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{4}.evaluated = c_files;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{5}.evaluated = [];
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{6}.evaluated = thresh;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{7}.evaluated = ex_var;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{8}.string    = 'noise_wm_csf_';
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'extract_segment';
            gi    = gi + 1;
            %extract_roi(epi_files,wm_file,epi_2_templ,[27 -43 16.5 4;-25.5 -46 16 4],ex_var); %coord from Bjoern
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.evaluated = epi_files;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.string    = char(wm_file);
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{3}.string    = char(epi_2_templ);
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{4}.evaluated = [27 -43 16.5 4;-25.5 -46 16 4];
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{5}.evaluated = ex_var;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{6}.string = 'noise_roi_';
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'extract_roi';
            gi    = gi + 1;
        end
    end
    
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
