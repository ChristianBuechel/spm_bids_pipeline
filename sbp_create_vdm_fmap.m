function create_fieldmap(all_sub_ids)
% function create_fieldmap(all_sub_ids)
% takes magnitude and phasediff images to estimate a voxel displacement map
% (vdm) to be used to correct EPI distortions
% it provides vdm5_asub... files for each 4D epi file (run) in the func directory
% see https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;d2c47aa4.1606 for
% details

[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    all_meta_pd    = spm_BIDS(BIDS,'metadata','sub',sprintf('%02d',sub_id),'type','phasediff'); % get the 2 echo times from the json files
    if numel(all_meta_pd) > 1
        meta_pd = all_meta_pd{1};
    else
        meta_pd = all_meta_pd;
    end
    meta_bold  = spm_BIDS(BIDS,'metadata','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'type','bold'); % get total readout time
    epi        = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    epi        = epi(1); % if there are brain and spinal, just take brain
    mean_file  = spm_file(epi,'prefix','meana');    % we need meanasub (for session1 and first epi images for all other sessions
    mean_dir   = spm_file(mean_file,'path');
    
    cnt        = 1;
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            epi                =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            epi                =  epi(1); % if there are brain and spinal, just take brain
            run_niftis         =  spm_select('ExtFPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^a'),1); %first EPI in each run
            all_niftis{ses,run}=  cellstr(run_niftis);
            pm                 =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'type','phasediff');
            have_pm(ses,run)   =  ~isempty(pm); % binary matrix to quickly see which PMs we have
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
    template.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
    template.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {'C:\Users\buechel\Documents\MATLAB\spm\toolbox\FieldMap\T1.nii'};
    % we have 2 cases: 1) a PM for each run 2) one PM per session
    if sum(have_pm(:)) == prod(size(have_pm)) % one pm per run
        gi  = 1;
        for ses = 1:vars.nSess
            for run = 1:vars.nRuns
                matlabbatch{gi,sub} = template;
                pd = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'type','phasediff');
                fmap_dir = spm_file(pd,'path');
                matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = pd;
                magn_f = strrep(pd,'phasediff','magnitude1'); % hack, bc spm_BIDS does not reveal the magnitude image name although its in there?
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
                % e.g. vdm5_scsub-92_ses-01_run-01_phasediff_run1.nii in fmap
                % and want: vdm5_asub-88_ses-01_task-pdev_run-01_bold.nii in func
                
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = fullfile(fmap_dir,sprintf('vdm5_scsub-%02d_ses-%02d_run-%02d_phasediff.nii',sub_id,ses,run)); % move and rename vdms
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = func_dir;
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(1).pattern = 'vdm5_sc';
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(1).repl = 'vdm5_a';
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(2).pattern = sprintf('_run-%02d_phasediff',run);
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(2).repl = sprintf('_task-%s_run-%02d_bold',vars.task,run);
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
            magn_f = strrep(pd,'phasediff','magnitude1');
            matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = cellstr(magn_f);
            for run = 1:vars.nRuns
                matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.session(run).epi = all_niftis{ses,run};
            end
            if ses == 1
                matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.session(1).epi = mean_file; %unwarp mean instead of 1st EPI bc this is what we coregister to T1
            end
            func_dir   = spm_file(all_niftis{ses,run},'path');
            matlabbatch{gi,sub}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'run';
            gi = gi + 1;
            
            % now we need to rename the vdm files in the batch system
            % we have vdm5_scsub-88_ses-01_run-01_phasediff_run1.nii, ..run2
            % etc in fmap
            % and want: vdm5_asub-88_ses-01_task-pdev_run-01_bold.nii in func
            %           vdm5_asub-88_ses-01_task-pdev_run-02_bold.nii in func
            % also rename umean... to cmean... (it is not yet coregistered but will be in the next step)
            for run = 1:vars.nRuns
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = fullfile(fmap_dir,sprintf('vdm5_scsub-%02d_ses-%02d_run-01_phasediff_run%d.nii',sub_id,ses,run)); % move and rename vdms
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = func_dir;
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(1).pattern = 'vdm5_sc';
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(1).repl = 'vdm5_a';
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(2).pattern = sprintf('_run-01_phasediff_run%d',run);
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(2).repl = sprintf('_task-%s_run-%02d_bold',vars.task,run);
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
