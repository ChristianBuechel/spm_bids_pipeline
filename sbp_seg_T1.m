function seg_T1(all_sub_ids)
% function seg_T1
% performs a simple segmentatoon of the T1 image and does the Dartel export

% resolve paths
[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
      
    struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'type','T1w');
    struc_file = struc_file(1); %in case there are more sessions with T1w
    
    gi    = 1;
    matlabbatch{gi,sub}.spm.spatial.preproc.channel.vols = cellstr(struc_file); % T1 image gets segmented 
    matlabbatch{gi,sub}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{gi,sub}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{gi,sub}.spm.spatial.preproc.channel.write = [0 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(1).tpm = {fullfile(path.templateDir,'enhanced_TPM.nii,1')}; 
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(1).native = [1 1]; % also dartel prepared
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(2).tpm = {fullfile(path.templateDir,'enhanced_TPM.nii,2')};
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(2).native = [1 1]; % also dartel prepared
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(3).tpm = {fullfile(path.templateDir,'enhanced_TPM.nii,3')};
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(4).tpm = {fullfile(path.templateDir,'enhanced_TPM.nii,4')};
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(5).tpm = {fullfile(path.templateDir,'enhanced_TPM.nii,5')};
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(5).native = [0 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(6).tpm = {fullfile(path.templateDir,'enhanced_TPM.nii,6')};
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{gi,sub}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{gi,sub}.spm.spatial.preproc.warp.reg    = [0 0.001 0.5 0.05 0.2];
    matlabbatch{gi,sub}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{gi,sub}.spm.spatial.preproc.warp.fwhm   = 0;
    matlabbatch{gi,sub}.spm.spatial.preproc.warp.samp   = 3;
    matlabbatch{gi,sub}.spm.spatial.preproc.warp.write  = [0 0]; %we do not need deformation fields
    gi = gi + 1;    
end

% run matlabbatch
n_procs = n_subs; % to not block too many cores on the server

if n_procs > vars.max_procs
    n_procs = vars.max_procs;
end

if vars.parallel == 1
    run_spm_parallel(matlabbatch, n_procs);
    %run_spm_multiple(matlabbatch, n_procs);
else
    run_spm_sequential(matlabbatch);
end
