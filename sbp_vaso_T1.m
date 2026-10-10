function sbp_vaso_T1(all_sub_ids)

% function sbp_vaso_T1(all_sub_ids)
% Preprocesses the short VASO reference runs (task vars.ref_task, default 'vasoref': one cbv
% (blood nulled) and one bold series per session, each with magnitude and phase) and builds
% T1-like images in the anat folder. They have the contrast of the T1w image and the
% distortions of the EPI, to improve the registration between the EPI and the T1w.
%
% For each session (jobs of the SPM batch system, subjects in parallel):
%  1) NORDIC on both series (magnitude + phase, output prefix 'n': n<name>.nii). The images
%     are 3D EPI, there is no slice timing. Independent of vars.preproc.steps.
%  2) realignment of the magnitude images of each series separately (estimate, only the mean
%     image is written: meann<name>.nii; phase images do not realign well because of wraps)
%  3) coregistration of the mean cbv to the mean bold, applied to all cbv magnitude images
%  4) the realignment/coregistration headers (sidecars n<name>.mat) are copied to the
%     phase images (acquired together with the magnitudes, stored as separate series)
%  5) vaso_t1_images: images in the anat folder, BIDS names with an acq label:
%       sub-XX_ses-YY_acq-vasomedian_T1w.nii   median(cbv)/median(bold)
%       sub-XX_ses-YY_acq-vasomeanstd_T1w.nii  mean/std over all volumes
%       sub-XX_ses-YY_acq-vasosigned_T1w.nii   sign restored with the phase (CSF negative)
%       sub-XX_ses-YY_acq-vasosignedpos_T1w.nii the same, positive values (0..1000)
%     They are marked as derived (descrip), sbp_align_images does not move them again.
%
% Settings in get_study_specs (all optional):
%   vars.nordic.vasoref_noisevols  number of noise volumes (no excitation) at the END of each series, or 'auto'
%                          (default: detected, see detect_noise_volumes).
%                          NORDIC splits them off: noise_<name>.nii/.json (noise level, used by sbp_nordic
%                          for runs without noise volumes), n<name> contains the signal volumes only
%   vars.nordic.func_noisevols  same for the main runs (sbp_nordic), default 'auto'
%   vars.nordic.ref_noise  0: sbp_nordic does not use the noise level of the vasoref run (1)
%   vars.ref_task          task label of the reference runs (default 'vasoref')
%   vars.vaso.thresh       relative signal threshold for the ratios (0.05)
%   vars.vaso.phase_fwhm   smoothing [mm] of the complex phase reference (0)
%   vars.vaso.phase_scale  phase value to radians (pi/4096, Siemens)

specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

if isfield(vars,'ref_task'), ref_task = vars.ref_task; else, ref_task = 'vasoref'; end
vaso = struct();
if isfield(vars,'vaso'), vaso = vars.vaso; end
n_noise = 'auto';   % number of noise volumes at the end of each series, or 'auto' (detect_noise_volumes)
if isfield(vars,'nordic') && isfield(vars.nordic,'vasoref_noisevols'), n_noise = vars.nordic.vasoref_noisevols; end

BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

arg = [];
arg.temporal_phase     = 1;
arg.phase_filter_width = 10;
arg.make_complex_nii   = 1; % also write the denoised phase
% NORDIC splits the noise volumes at the end off (noise_<name>.nii/.json next to the raw images) and
% measures the noise level on them. Without, it assumes unit noise, which was ~11% too low for this
% sequence: tSNR 2.3 instead of 4.2 (sub-03 test). n<name> has no noise volumes.
arg.noise_volume_split = n_noise;

for sub = 1:n_subs
    sub_id = all_sub_ids(sub);
    gi     = 1;
    for ses = 1:vars.nSess
        % the four images of this session: cbv/bold, magnitude/phase
        q  = {'task',ref_task,'expect',1};
        mag_c = char(get_func_files(BIDS,sub_id,ses,[],q{:},'type','cbv', 'part','mag'));
        ph_c  = char(get_func_files(BIDS,sub_id,ses,[],q{:},'type','cbv', 'part','phase'));
        mag_b = char(get_func_files(BIDS,sub_id,ses,[],q{:},'type','bold','part','mag'));
        ph_b  = char(get_func_files(BIDS,sub_id,ses,[],q{:},'type','bold','part','phase'));

        n_mag_c = spm_file(mag_c,'prefix','n');  n_ph_c = spm_file(ph_c,'prefix','n');
        n_mag_b = spm_file(mag_b,'prefix','n');  n_ph_b = spm_file(ph_b,'prefix','n');
        mean_c  = spm_file(n_mag_c,'prefix','mean');
        mean_b  = spm_file(n_mag_b,'prefix','mean');

        % 1) NORDIC (splits off the noise volumes: n<name> has signal volumes only, so the steps below use the whole series)
        for ser = {{mag_c,ph_c},{mag_b,ph_b}}
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.string    = ser{1}{1};
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.string    = ser{1}{2};
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{3}.string    = 'n';
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{4}.evaluated = arg;
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'NIFTI_NORDIC_fast';
            gi = gi + 1;
        end

        % 2) realignment of the magnitudes, each series for itself
        for ser = {n_mag_c,n_mag_b}
            matlabbatch{gi,sub}.spm.spatial.realign.estwrite.data             = {ser(1)}; % whole 4D series (all frames)
            matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.quality = 0.95;
            matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.sep     = 1.5;
            matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.fwhm    = 1;
            matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.rtm     = 1; % two pass, register to the mean
            matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.interp  = 2;
            matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
            matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.weight  = '';
            matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.which   = [0 1]; % only the mean
            matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.interp  = 4;
            matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
            matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.mask    = 1;
            matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.prefix  = 'r';
            gi = gi + 1;
        end

        % 3) coregister the mean cbv to the mean bold, apply to all cbv magnitude images
        matlabbatch{gi,sub}.spm.spatial.coreg.estimate.ref    = {mean_b};
        matlabbatch{gi,sub}.spm.spatial.coreg.estimate.source = {mean_c};
        matlabbatch{gi,sub}.spm.spatial.coreg.estimate.other  = {n_mag_c}; % all frames
        matlabbatch{gi,sub}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{gi,sub}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
        matlabbatch{gi,sub}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{gi,sub}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
        gi = gi + 1;

        % 4) headers of the magnitudes to the phase images
        matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.evaluated = {n_mag_c;n_mag_b};
        matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.evaluated = {n_ph_c;n_ph_b};
        matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
        matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'copy_sidecar';
        gi = gi + 1;

        % 5) the T1-like images
        anat_dir = char(spm_file(get_anat(BIDS,sub_id,ses),'path'));
        w = struct('mag_c',n_mag_c,'ph_c',n_ph_c,'mag_b',n_mag_b,'ph_b',n_ph_b,'ref',mean_b); % (no noise volumes in n<name>)
        for f = fieldnames(vaso)', w.(f{1}) = vaso.(f{1}); end
        for lab = {'median','meanstd','signed','signedpos'}
            w.out.(lab{1}) = fullfile(anat_dir,bids_name(sub_id,ses,'T1w','acq',['vaso' lab{1}]));
        end
        matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.evaluated = w;
        matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
        matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'vaso_t1_images';
        gi = gi + 1;
    end
end

% run matlabbatch
n_procs = n_subs; % to not block to many cores on the server

if n_procs > vars.max_procs
    n_procs = vars.max_procs;
end

if vars.parallel == 1
    run_spm_parallel(matlabbatch, n_procs);
else
    run_spm_sequential(matlabbatch);
end
