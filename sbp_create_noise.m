function sbp_create_noise(all_sub_ids)
% function sbp_create_noise(all_sub_ids)
% creates noise regressors for 1st level analysis
% 1) WM and CSF (principal components of the time series in the eroded c2 / c3 cores of the T1w)
% 2) ROI at the posterior tip of the lateral ventricles (as in Horing et al 202X); ROIs that are outside
%    the field of view of the EPI (slab) give empty regressors and a message (see extract_roi)
% Works on the resliced images (r<P><epi>.nii from sbp_realign_2_2 / sbp_apply_vdm, P: see epi_prefix)
% and writes noise_wm_csf_r<P><epi>.mat and noise_roi_r<P><epi>.mat next to them
% (names: epi_file(epi,vars,'noise_wm_csf') and 'noise_roi').
% Needs c2<T1w>, c3<T1w>, cmean<P><epi> and y_epi_2_T1.nii (and wmean<P><epi>, y_epi_2_template.nii for the ROIs).
% Optional settings (get_study_specs) for the WM / CSF components (see extract_segment):
%   vars.noise.smooth_fwhm  smoothing [mm] of the time series before the components are extracted (0)
%   vars.noise.n_erode      number of erosions of the c2 / c3 cores (1)
%   all_sub_ids   subject numbers, e.g. 1:5 (run in parallel if vars.parallel = 1)

% add paths
specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

smooth_fwhm = 0;  n_erode = 1;           % defaults of extract_segment (as before)
if isfield(vars,'noise')
    if isfield(vars.noise,'smooth_fwhm'), smooth_fwhm = vars.noise.smooth_fwhm; end
    if isfield(vars.noise,'n_erode'),     n_erode     = vars.noise.n_erode;     end
end

for sub = 1:n_subs
    sub_id    = all_sub_ids(sub);
    
    struc_file = get_anat(BIDS,sub_id);
    anat_dir   = spm_file(struc_file,'path');
    
    epi        = get_ref_epi(BIDS,sub_id);

    func_dir   = spm_file(epi,'path');
    
    c_files{1}  = spm_file(struc_file,'prefix','c2');
    c_files{2}  = spm_file(struc_file,'prefix','c3');
    
    wm_file        = epi_file(epi,vars,'wmean'); % mean EPI in template space (sbp_warp / sbp_warp_vdm), P: see epi_prefix
    cm_file        = epi_file(epi,vars,'cmean'); % mean EPI in T1 space (sbp_coreg_nlin / sbp_coreg_lin)
    epi_2_T1       = fullfile(func_dir, 'y_epi_2_T1.nii');
    epi_2_templ    = fullfile(func_dir, 'y_epi_2_template.nii');
    
    
    thresh = [0.9 0.9]; ex_var = 6;
    gi = 1;
    for ses = 1:vars.nSess
        for run = get_run_numbers(BIDS,sub_id,ses)
            fprintf('ses %d run %d\n',ses,run);
            epi       =  get_func_files(BIDS,sub_id,ses,run);
            epi_files =  epi_frames(epi,vars,'res'); % resliced frames r<P><epi>; the noise files are named after them (epi_file kinds noise_*)
            %extract_segment(epi_files,cm_file,epi_2_T1,c_files,[],thresh,ex_var);
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.evaluated = epi_files;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.string    = char(cm_file);
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{3}.string    = char(epi_2_T1);
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{4}.evaluated = c_files;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{5}.evaluated = [];
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{6}.evaluated = thresh;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{7}.evaluated = ex_var;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{8}.string    = 'noise_wm_csf_';
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{9}.evaluated  = smooth_fwhm;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{10}.evaluated = n_erode;
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
