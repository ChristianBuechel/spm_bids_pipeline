function sbp_coreg_nlin_vaso_T1(all_sub_ids)
% function sbp_coreg_nlin_vaso_T1(all_sub_ids)
% Nonlinear coregistration of the EPI to the T1w image, using the signed VASO T1-like image
% (sbp_vaso_T1) as the source of the normalisation instead of the mean EPI. That image has the
% contrast of the T1w (CSF < GM < WM) and the distortions of the EPI. Replaces sbp_coreg_nlin
% (same outputs, so the following steps are unchanged).
%
% Per subject (session 1):
%  1) segmentation of the T1w (as in sbp_coreg_nlin: 6 classes, c1..c4, rc1/rc2 for Dartel),
%     the segments are smoothed (vars.coreg_nlin_vaso.tpm_smooth, 1 mm) and merged into a TPM
%  2) the signed vaso image is clipped to the percentile range of its non-zero voxels
%     (clip_signed_image, clip_<name>.nii). Later this will be the last step of sbp_vaso_T1
%  3) SPM's normalisation (= spm_preproc8 with intensity histograms instead of Gaussians per class)
%     of the clipped image to the TPM: deformation field y_clip_<name>.nii (T1 grid -> vaso image)
%     the warped image c<clip_name>.nii (in T1 space, for checking) is kept
%  4) the deformation of the mean EPI of the main runs: rigid registration of mean<P><epi> onto the
%     mean bold of the vasoref run, applied to the deformation of 3) (epi_def_from_vaso):
%     y_mean<P><epi>.nii. The mean EPI is written in T1 space: cmean<P><epi>.nii (P: see epi_prefix)
%
% Settings (get_study_specs, all optional, vars.coreg_nlin_vaso.<field>), tested on sub-01 / sub-03:
%   clip_pct     [0.5 99.5]  percentile range of the clipping
%   reg_scale    10          scaling of the regularisation of the warp [0 1e-5 0.005 0.0005 0.002]; larger = stiffer
%   samp         1.5         sampling distance [mm] of the warp (SPM default 3)
%   fwhm         1.5         smoothing [mm] of the image for the registration (default 3)
%   tpm_smooth   1           smoothing [mm] of the T1 segments that make the TPM (sbp_coreg_nlin: 3)
%   biasfwhm     Inf         no bias correction (the cbv/bold ratio image has no receive bias)
%   biasreg      0.0001
%   affreg       ''          no affine step (the images are aligned by the headers; SPM: '' = none, 'none' = unregularised)
%   vaso_image   'signed'    which vaso T1-like image (acq-vaso<..>): 'signed', 'signedpos', ...
%   ref_task     vars.ref_task ('vasoref'), task label of the reference runs

tmpTPM = 'TPM_tmp.nii';

% resolve paths
specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

% parameters (defaults, see above)
p = struct('clip_pct',[0.5 99.5],'reg_scale',10,'samp',1.5,'fwhm',1.5,'tpm_smooth',1,'biasfwhm',Inf,'biasreg',0.0001,'affreg','','vaso_image','signed');
if isfield(vars,'coreg_nlin_vaso')
    for f = fieldnames(vars.coreg_nlin_vaso)', p.(f{1}) = vars.coreg_nlin_vaso.(f{1}); end
end
if isfield(vars,'ref_task'), ref_task = vars.ref_task; else, ref_task = 'vasoref'; end
if isfield(p,'ref_task'), ref_task = p.ref_task; end

BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);

    struc_file = get_anat(BIDS,sub_id);
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

    epi         = get_ref_epi(BIDS,sub_id);
    mean_file   = epi_file(epi,vars,'mean');   % mean<P><epi> of the main runs (realign 1)
    y_file      = epi_file(epi,vars,'ymean');  % output: y_mean<P><epi>
    % (cmean<P><epi> is written by SPM as c + the name of mean_file)

    % the vaso image and the mean bold of the vasoref run (sbp_vaso_T1)
    vaso_file   = fullfile(char(anat_dir),bids_name(sub_id,1,'T1w','acq',['vaso' p.vaso_image]));
    clip_file   = spm_file(vaso_file,'prefix','clip_');
    y_vaso      = spm_file(clip_file,'prefix','y_');     % written by the normalisation
    ref_raw     = get_func_files(BIDS,sub_id,1,[],'task',ref_task,'type','bold','part','mag','expect',1);
    ref_bold    = spm_file(ref_raw,'prefix','meann');    % mean bold of the vasoref run (sbp_vaso_T1: NORDIC n, realign mean)
    if ~isfile(vaso_file), error('sub-%02d: %s not found, run sbp_vaso_T1 first',sub_id,vaso_file); end

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

    matlabbatch{gi,sub}.spm.spatial.smooth.data(1,1) = c1_file; % Now smooth segments (less than in sbp_coreg_nlin: a sharper template)
    matlabbatch{gi,sub}.spm.spatial.smooth.data(2,1) = c2_file;
    matlabbatch{gi,sub}.spm.spatial.smooth.data(3,1) = c3_file;
    matlabbatch{gi,sub}.spm.spatial.smooth.data(4,1) = c4_file;
    matlabbatch{gi,sub}.spm.spatial.smooth.data(5,1) = c5_file;
    matlabbatch{gi,sub}.spm.spatial.smooth.data(6,1) = c6_file;
    matlabbatch{gi,sub}.spm.spatial.smooth.fwhm   = p.tpm_smooth*[1 1 1];
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

    % clip the signed vaso image
    matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.string    = vaso_file;
    matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.string    = clip_file;
    matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{3}.evaluated = p.clip_pct;
    matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
    matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'clip_signed_image';
    gi = gi + 1;

    % normalise the clipped vaso image to the TPM (deformation y_clip_<name>, warped image c<clip_name>)
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.subj.vol = cellstr(clip_file);
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.subj.resample = cellstr(clip_file);
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.eoptions.biasreg = p.biasreg;
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = p.biasfwhm;
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.eoptions.tpm(1) = fullfile(anat_dir,tmpTPM);
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.eoptions.affreg = p.affreg;
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.eoptions.reg = p.reg_scale*[0 1e-05 0.005 0.0005 0.002]; % larger: nonlinear less aggressive
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.eoptions.fwhm = p.fwhm;
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.eoptions.samp = p.samp;
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.woptions.bb = [NaN NaN NaN; NaN NaN NaN];
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.woptions.vox = [NaN NaN NaN];
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    matlabbatch{gi,sub}.spm.spatial.normalise.estwrite.woptions.prefix = 'c';
    gi = gi + 1;

    % deformation of the mean EPI of the main runs: rigid registration to the vasoref mean bold, applied to the vaso deformation
    matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = char(y_vaso);
    matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.string = char(mean_file);
    matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{3}.string = char(ref_bold);
    matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{4}.string = char(y_file);
    matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
    matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'epi_def_from_vaso';
    gi = gi + 1;

    % mean EPI in T1 space (cmean<P><epi>)
    matlabbatch{gi,sub}.spm.spatial.normalise.write.subj.def = cellstr(y_file);
    matlabbatch{gi,sub}.spm.spatial.normalise.write.subj.resample = cellstr(mean_file);
    matlabbatch{gi,sub}.spm.spatial.normalise.write.woptions.bb = [NaN NaN NaN; NaN NaN NaN];
    matlabbatch{gi,sub}.spm.spatial.normalise.write.woptions.vox = [NaN NaN NaN];
    matlabbatch{gi,sub}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{gi,sub}.spm.spatial.normalise.write.woptions.prefix = 'c';
    gi = gi + 1;

    % delete all the stuff we do not need anymore (the large vaso deformation field is not needed after step 4)
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = [c5_file; c6_file; all_sc; cellstr(fullfile(anat_dir,tmpTPM)); cellstr(y_vaso)];
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
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
