function create_fieldmap(all_sub_ids)
% function create_fieldmap(all_sub_ids)
% takes magnitude and phasediff images to estimate a voxel displacement map
% (vdm) to be used to correct EPI distortions
% it provides vdm5_<P><epi> files for each 4D epi file (run) in the func directory,
% i.e. with the name of the EPI (e.g. vdm5_ansub-01_ses-01_task-x_run-01_part-mag_bold.nii,
% P: prefix of the realigned images, see epi_prefix and epi_file(...,'vdm')).
% It also writes the unwarped mean umean<P><epi> (the mean of realign 1 is unwarped).
% see https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;d2c47aa4.1606 for
% details

specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

% check for FieldMaop toolbox and get T1 image
all_tb        = spm('TBs');
all_tb_names  = cellstr(strvcat(all_tb.name));
FM_tb_present = any(contains(all_tb_names,'FieldMap'));
if ~FM_tb_present
    error('FieldMap toolbox missing. Please install the toolbox');
else
    FM_tb_dir = all_tb(contains(all_tb_names,'FieldMap')).dir;
    FM_templ  = fullfile(FM_tb_dir,'T1.nii');
end

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    have_pm    = false(0); % reset for each subject, (ses,run) matrix, which runs have a phasediff
    have_run   = false(0); % (ses,run) matrix, which runs exist
    all_niftis = {};
    all_epis   = {};
    all_meta_pd  = spm_BIDS(BIDS,'metadata','sub',sprintf('%02d',sub_id),'type','phasediff'); % get the 2 echo times from the json files
    if numel(all_meta_pd) > 1
        meta_pd = all_meta_pd{1};
    else
        meta_pd = all_meta_pd;
    end
    meta_bold  = spm_BIDS(BIDS,'metadata','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'acq',{'brain',''},'part',{'mag',''},'type',{'bold','cbv'}); % get total readout time (of the main task, brain, magnitude: one entry)
    epi        = get_ref_epi(BIDS,sub_id);
    mean_file  = epi_file(epi,vars,'mean');    % we need mean<P><epi> (for session1 and first epi images for all other sessions
    mean_dir   = spm_file(mean_file,'path');
    
    cnt        = 1;
    for ses = 1:vars.nSess
        for run = get_run_numbers(BIDS,sub_id,ses)
            epi                =  get_func_files(BIDS,sub_id,ses,run);
            run_niftis         =  epi_frames(epi,vars,'in',1); %first EPI in each run (P<epi>)
            all_niftis{ses,run}=  cellstr(run_niftis);
            all_epis{ses,run}  =  epi;                         % raw BIDS name, defines the name of the vdm file
            pm                 =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'type','phasediff');
            have_pm(ses,run)   =  ~isempty(pm); % binary matrix to quickly see which PMs we have
            have_run(ses,run)  =  true;
            cnt                =  cnt + 1;
        end
    end
    % define default options
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 1;
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1;
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = meta_bold.TotalReadoutTime.*1000; % assume is constant over ses and runs
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [meta_pd.EchoTime1.*1000 meta_pd.EchoTime2.*1000]; % dito
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
    template.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 1;
    template.spm.tools.fieldmap.calculatevdm.subj.anat = {''};
    template.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
    template.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 0;
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {FM_templ};
    % we have 2 cases: 1) a PM for each run 2) one PM per session
    if all(have_pm(have_run)) % one pm per run (of the runs that exist)
        gi  = 1;
        for ses = 1:vars.nSess
            for run = find(have_run(ses,:))
                matlabbatch{gi,sub} = template;
                pd = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'type','phasediff');
                fmap_dir = spm_file(pd,'path');
                matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = pd;
                magn_f = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'type','magnitude1');
                matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = cellstr(magn_f);
                if ses == 1 & run == 1
                    matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.session(1).epi = mean_file;
                else
                    matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.session(1).epi = all_niftis{ses,run};
                end
                func_dir   = spm_file(all_niftis{ses,run},'path');
                matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'run';
                gi = gi + 1;
                % now we have to move and rename the vdm file
                % e.g. vdm5_scsub-92_ses-01_run-01_phasediff.nii in fmap
                % and want: vdm5_<P>sub-88_ses-01_task-pdev_run-01_part-mag_bold.nii in func (epi_file(epi,vars,'vdm'))

                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = spm_file(pd,'prefix','vdm5_sc'); % move and rename vdms
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = func_dir;
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(1).pattern = '^.*$'; % the whole name is replaced ...
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(1).repl = char(spm_file(epi_file(all_epis{ses,run},vars,'vdm'),'basename')); % ... by the name of the vdm of this EPI
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
                gi = gi + 1;
            end
        end        
    elseif sum(have_pm(:,1)) == size(have_pm,1) % one pm per session
        gi  = 1;
        for ses = 1:vars.nSess
            matlabbatch{gi,sub} = template;
            pd = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',1),'type','phasediff');
            fmap_dir = spm_file(pd,'path');
            matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = pd;
            magn_f = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',1),'type','magnitude1');
            matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = cellstr(magn_f);
            ses_runs = find(have_run(ses,:)); % runs that exist in this session
            for k = 1:numel(ses_runs)
                matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.session(k).epi = all_niftis{ses,ses_runs(k)};
            end
            run = ses_runs(end); % (used below for func_dir)
            if ses == 1
                matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.session(1).epi = mean_file; % unwarp mean instead of 1st EPI bc this is what we coregister to T1
            end
            func_dir   = spm_file(all_niftis{ses,run},'path');
            matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'run';
            gi = gi + 1;
            
            % now we need to rename the vdm files in the batch system
            % we have vdm5_scsub-88_ses-01_run-01_phasediff_run1.nii, ..run2
            % etc in fmap
            % and want: vdm5_<P>sub-88_ses-01_task-pdev_run-01_part-mag_bold.nii in func
            %           vdm5_<P>sub-88_ses-01_task-pdev_run-02_part-mag_bold.nii in func
            % also rename umean... to cmean... (it is not yet coregistered but will be in the next step)
            for k = 1:numel(ses_runs) % SPM numbers the vdm files by the session index k, we want the run number
                run = ses_runs(k);
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = spm_file(pd,'prefix','vdm5_sc','suffix',sprintf('_run%d',k)); % move and rename vdms
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = func_dir;
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(1).pattern = '^.*$'; % the whole name is replaced ...
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(1).repl = char(spm_file(epi_file(all_epis{ses,run},vars,'vdm'),'basename')); % ... by the name of the vdm of this EPI
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
                gi = gi + 1;
            end
        end
    else
        error('need at least one phasediff/magnitude image per sesssion')
    end
    
end
%save('create_fieldmap','matlabbatch');
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
