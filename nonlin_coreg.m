function nonlin_coreg(all_sub_ids)
% function nonlin_coreg
% performs a nonlinear coregistration between mean EPI and T1
% writes a cmean* file that is the m,ean EPI in the space of the T1
% also does a segmentation of the T1, and the Dartel export, so we don't
% have to do this again

tmpTPM = 'TPM_tmp.nii';

% resolve paths
[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);
run_parallel = 1;

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
      
    struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'type','T1w');
    anat_dir   = spm_file(struc_file,'path');
    
    c1_file    = spm_file(struc_file, 'prefix','c1');
    c2_file    = spm_file(struc_file, 'prefix','c2');
    c3_file    = spm_file(struc_file, 'prefix','c3');
    c4_file    = spm_file(struc_file, 'prefix','c4');
    c5_file    = spm_file(struc_file, 'prefix','c5');
    c6_file    = spm_file(struc_file, 'prefix','c6');
    
    sc1_file    = spm_file(struc_file, 'prefix','sc1');
    sc2_file    = spm_file(struc_file, 'prefix','sc2');
    sc3_file    = spm_file(struc_file, 'prefix','sc3');
    sc4_file    = spm_file(struc_file, 'prefix','sc4');
    sc5_file    = spm_file(struc_file, 'prefix','sc5');
    sc6_file    = spm_file(struc_file, 'prefix','sc6');
    
    epi         = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    mean_file   = spm_file(epi,'prefix','meana');

    
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
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(6).tpm = {fullfile(path.templateDir,'enhanced_TPM.nii,6')};
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(6).native = [1 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{gi,sub}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{gi,sub}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{gi,sub}.spm.spatial.preproc.warp.reg    = [0 0.001 0.5 0.05 0.2];
    matlabbatch{gi,sub}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{gi,sub}.spm.spatial.preproc.warp.fwhm   = 0;
    matlabbatch{gi,sub}.spm.spatial.preproc.warp.samp   = 3;
    matlabbatch{gi,sub}.spm.spatial.preproc.warp.write  = [0 0]; %we do not need deformation fields
    gi = gi + 1;
    
    matlabbatch{gi,sub}.spm.spatial.smooth.data(1,1) = c1_file; % Now smooth segments
    matlabbatch{gi,sub}.spm.spatial.smooth.data(2,1) = c2_file;
    matlabbatch{gi,sub}.spm.spatial.smooth.data(3,1) = c3_file;
    matlabbatch{gi,sub}.spm.spatial.smooth.data(4,1) = c4_file;
    matlabbatch{gi,sub}.spm.spatial.smooth.data(5,1) = c5_file;
    matlabbatch{gi,sub}.spm.spatial.smooth.data(6,1) = c6_file;
    matlabbatch{gi,sub}.spm.spatial.smooth.fwhm   = [3 3 3];
    matlabbatch{gi,sub}.spm.spatial.smooth.dtype  = 0;
    matlabbatch{gi,sub}.spm.spatial.smooth.im     = 0;
    matlabbatch{gi,sub}.spm.spatial.smooth.prefix = 's';
    gi = gi + 1;
    
    matlabbatch{gi,sub}.spm.util.cat.vols(1,1) = sc1_file; % assemble all segments into a 4D file which is a TPM for later
    matlabbatch{gi,sub}.spm.util.cat.vols(2,1) = sc2_file;
    matlabbatch{gi,sub}.spm.util.cat.vols(3,1) = sc3_file;
    matlabbatch{gi,sub}.spm.util.cat.vols(4,1) = sc4_file;
    matlabbatch{gi,sub}.spm.util.cat.vols(5,1) = sc5_file;
    matlabbatch{gi,sub}.spm.util.cat.vols(6,1) = sc6_file;
    all_sc = matlabbatch{gi,sub}.spm.util.cat.vols; % handy for deletion
    matlabbatch{gi,sub}.spm.util.cat.name  = tmpTPM;
    matlabbatch{gi,sub}.spm.util.cat.dtype = 16;
    gi = gi + 1;
    
    reg = 1;    
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.subj.vol = cellstr(mean_file); % now normalize mean EPI using the generated TPM 
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.subj.resample = cellstr(mean_file);
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.eoptions.tpm(1) = fullfile(anat_dir,tmpTPM);
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.eoptions.affreg = 'subj';
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.eoptions.reg = reg*[0 1e-05 0.005 0.0005 0.002];%if reg is bigger nonlinear will be less aggressive
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.eoptions.fwhm = 3;
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.woptions.bb = [NaN NaN NaN; NaN NaN NaN];
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.woptions.vox = [1.5 1.5 1.5];
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.woptions.prefix = 'c';
    gi = gi + 1;
    
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = [c5_file; c6_file; all_sc; fullfile(anat_dir,tmpTPM)];
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false; % delete all the stuff we do not need anymore
    gi = gi + 1;
end

% run matlabbatch
n_procs = n_subs; % to not block too many cores on the server

if n_procs > vars.max_procs
    n_procs = vars.max_procs;
end

if run_parallel == 1
    run_spm_parallel(matlabbatch, n_procs);
else
    spm_jobman('run',matlabbatch);
end
end






