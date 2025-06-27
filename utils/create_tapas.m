% Makes Physio Regressors
function create_tapas(all_sub_ids)

[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    gi = 1;

    for ses = 1:vars.nSess
        for run = 1:4
            epi = spm_BIDS(BIDS,'data','sub',sprintf('%02d',all_sub_ids(sub)),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            nScans = numel(spm_vol(epi{:,1}));
            % make batch input
            matlabbatch{gi,sub}.spm.tools.physio.save_dir = {''};
            matlabbatch{gi,sub}.spm.tools.physio.log_files.vendor = 'BIDS';
            matlabbatch{gi,sub}.spm.tools.physio.log_files.cardiac = {fullfile(path.preprocDir, sprintf('sub-%02d/ses-%02d/func/sub-%02d_ses-%02d_task-%s_run-%02d_physio.tsv',all_sub_ids(sub),ses,all_sub_ids(sub),ses,vars.task,run))};
            matlabbatch{gi,sub}.spm.tools.physio.log_files.respiration = {fullfile(path.preprocDir, sprintf('sub-%02d/ses-%02d/func/sub-%02d_ses-%02d_task-%s_run-%02d_physio.tsv',all_sub_ids(sub),ses,all_sub_ids(sub),ses,vars.task,run))};
            matlabbatch{gi,sub}.spm.tools.physio.log_files.scan_timing = {''};
            matlabbatch{gi,sub}.spm.tools.physio.log_files.sampling_interval = [];
            matlabbatch{gi,sub}.spm.tools.physio.log_files.relative_start_acquisition = 0;
            matlabbatch{gi,sub}.spm.tools.physio.log_files.align_scan = 'first';
            matlabbatch{gi,sub}.spm.tools.physio.scan_timing.sqpar.Nslices = 60;
            matlabbatch{gi,sub}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
            matlabbatch{gi,sub}.spm.tools.physio.scan_timing.sqpar.TR = 1.8;
            matlabbatch{gi,sub}.spm.tools.physio.scan_timing.sqpar.Ndummies = 0;
            matlabbatch{gi,sub}.spm.tools.physio.scan_timing.sqpar.Nscans = nScans;
            matlabbatch{gi,sub}.spm.tools.physio.scan_timing.sqpar.onset_slice = 1;
            matlabbatch{gi,sub}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
            matlabbatch{gi,sub}.spm.tools.physio.scan_timing.sqpar.Nprep = 0;
            matlabbatch{gi,sub}.spm.tools.physio.scan_timing.sync.nominal = struct([]);
            matlabbatch{gi,sub}.spm.tools.physio.preproc.cardiac.modality = 'PPU';
            matlabbatch{gi,sub}.spm.tools.physio.preproc.cardiac.filter.no = struct([]);
            matlabbatch{gi,sub}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
            matlabbatch{gi,sub}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
            matlabbatch{gi,sub}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.max_heart_rate_bpm = 90;
            matlabbatch{gi,sub}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
            matlabbatch{gi,sub}.spm.tools.physio.preproc.respiratory.filter.passband = [0.01 2];
            matlabbatch{gi,sub}.spm.tools.physio.preproc.respiratory.despike = true;
            matlabbatch{gi,sub}.spm.tools.physio.model.output_multiple_regressors = [path.preprocDir,sprintf('/sub-%02d/',all_sub_ids(sub)),sprintf('ses-%02d/',ses),'func/',sprintf('physio_sub-%02d_ses-%02d_task-%s_run-%02d_bold.txt',all_sub_ids(sub),ses,vars.task,run)];
            matlabbatch{gi,sub}.spm.tools.physio.model.output_physio = [path.preprocDir,sprintf('/sub-%02d/',all_sub_ids(sub)),sprintf('ses-%02d/',ses),'func/',sprintf('info_physio_sub-%02d_ses-%02d_task-%s_run-%02d.mat',all_sub_ids(sub),ses,vars.task,run)];
            %matlabbatch{gi,sub}.spm.tools.physio.model.output_physio = 0;
            matlabbatch{gi,sub}.spm.tools.physio.model.orthogonalise = 'none';
            matlabbatch{gi,sub}.spm.tools.physio.model.censor_unreliable_recording_intervals = false;
            matlabbatch{gi,sub}.spm.tools.physio.model.retroicor.yes.order.c = 3;
            matlabbatch{gi,sub}.spm.tools.physio.model.retroicor.yes.order.r = 4;
            matlabbatch{gi,sub}.spm.tools.physio.model.retroicor.yes.order.cr = 1;
            matlabbatch{gi,sub}.spm.tools.physio.model.rvt.no = struct([]);
            matlabbatch{gi,sub}.spm.tools.physio.model.hrv.no = struct([]);
            matlabbatch{gi,sub}.spm.tools.physio.model.noise_rois.no = struct([]);
            matlabbatch{gi,sub}.spm.tools.physio.model.movement.no = struct([]);
            matlabbatch{gi,sub}.spm.tools.physio.model.other.no = struct([]);
            matlabbatch{gi,sub}.spm.tools.physio.verbose.level = 0;
            matlabbatch{gi,sub}.spm.tools.physio.verbose.fig_output_file = '';
            matlabbatch{gi,sub}.spm.tools.physio.verbose.use_tabs = false;
            gi = gi+1;
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


end
