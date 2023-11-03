function clean_physio(all_sub_ids)
% clean_physio(all_sub_ids)
% residualize images with respect to physio/CSF noise
% to make realignment less prone to intesity artefacts

[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);
run_parallel = 1;

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    gi         = 1;
    epi                =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    func_dir           =  spm_file(epi,'path');
    mkdir(char(fullfile(func_dir,'res')));
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            epi                =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            rename_mat = 0;
            if exist(char(spm_file(epi,'prefix','a','ext','.mat')),'file') %assume original realign
                rename_mat = 1;
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = spm_file(epi,'prefix','a','ext','.mat');
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = func_dir;
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = 'asub';
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = 'brain_asub';
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
                gi = gi + 1;
            end
            
            matlabbatch{gi,sub}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = func_dir;
            matlabbatch{gi,sub}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'tmp';
            gi = gi + 1;
                     
            run_niftis         =  spm_select('ExtFPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^a'),Inf);
            all_nuis           = [];

            physio_noise_f     =  char(spm_file(epi,'prefix','physio_','ext','.mat'));
            if exist(physio_noise_f,'file')
                physio_noise       =  load(physio_noise_f);
                all_nuis           =  [all_nuis physio_noise.physio];
            end
            bs_csf_noise_f     =  char(spm_file(epi,'prefix','noise_bs_csf_a','ext','.mat'));
            if exist(bs_csf_noise_f,'file')
                bs_csf_noise       =  load(bs_csf_noise_f);
                all_nuis           =  [all_nuis bs_csf_noise.segment.data];
            end
            
            n_nuis             =  size(all_nuis,2);
            %do firstlevel with physio regressors only
            matlabbatch{gi,sub}.spm.stats.fmri_spec.sess.scans = cellstr(run_niftis);
            matlabbatch{gi,sub}.spm.stats.fmri_spec.dir = fullfile(func_dir,'tmp');
            matlabbatch{gi,sub}.spm.stats.fmri_spec.timing.units = 'scans';
            matlabbatch{gi,sub}.spm.stats.fmri_spec.timing.RT = vars.sliceTiming.tr;
            matlabbatch{gi,sub}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{gi,sub}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
            matlabbatch{gi,sub}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{gi,sub}.spm.stats.fmri_spec.sess.multi = {''};
            matlabbatch{gi,sub}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
            
            for nuis = 1:n_nuis
                matlabbatch{gi,sub}.spm.stats.fmri_spec.sess.regress(nuis) = struct('name', cellstr(num2str(nuis)), 'val', all_nuis(:,nuis));
            end
            
            matlabbatch{gi,sub}.spm.stats.fmri_spec.sess.multi_reg = {''};
            matlabbatch{gi,sub}.spm.stats.fmri_spec.sess.hpf = Inf; %no highpass filter
            matlabbatch{gi,sub}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{gi,sub}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{gi,sub}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{gi,sub}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{gi,sub}.spm.stats.fmri_spec.mthresh = -Inf;
            matlabbatch{gi,sub}.spm.stats.fmri_spec.mask = {''};
            matlabbatch{gi,sub}.spm.stats.fmri_spec.cvi = 'none';
            gi = gi + 1;
            
            matlabbatch{gi,sub}.spm.stats.fmri_est.spmmat = fullfile(fullfile(func_dir,'tmp'),'SPM.mat');
            matlabbatch{gi,sub}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{gi,sub}.spm.stats.fmri_est.method.Classical = 1;
            gi = gi + 1;
            
            matlabbatch{gi,sub}.spm.stats.con.spmmat = fullfile(fullfile(func_dir,'tmp'),'SPM.mat');
            matlabbatch{gi,sub}.spm.stats.con.consess{1}.fcon.name = 'mean';%the one we need to residualize with
            matlabbatch{gi,sub}.spm.stats.con.consess{1}.fcon.weights = [zeros(1,n_nuis) 1];
            matlabbatch{gi,sub}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
            
            matlabbatch{gi,sub}.spm.stats.con.consess{2}.fcon.name = 'nuis'; %F-map is good to see whether physio/CSF is correct
            matlabbatch{gi,sub}.spm.stats.con.consess{2}.fcon.weights = [eye(n_nuis) zeros(n_nuis,1)];
            matlabbatch{gi,sub}.spm.stats.con.consess{2}.fcon.sessrep = 'none';
            
            matlabbatch{gi,sub}.spm.stats.con.delete = 1;
            gi = gi + 1;
            
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.string    = char(fullfile(fullfile(func_dir,'tmp'),'SPM.mat'));
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.evaluated = 1; % the first contrast i.e. mean stays in
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'spm_write_residuals';
            gi = gi + 1;
            
            all_res = [];
            for r=1:size(run_niftis,1)
                all_res{r,1} = char(fullfile(func_dir,'tmp',sprintf('Res_%4.4d.nii',r)));
            end
            matlabbatch{gi,sub}.spm.util.cat.vols = all_res;
            matlabbatch{gi,sub}.spm.util.cat.name = char(spm_file(spm_file(epi,'filename'),'prefix','res_'));
            matlabbatch{gi,sub}.spm.util.cat.dtype = 4; %int16, SAME is 0
            matlabbatch{gi,sub}.spm.util.cat.RT = vars.sliceTiming.tr;
            gi = gi + 1;
            
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = fullfile(func_dir,'tmp',spm_file(spm_file(epi,'filename'),'prefix','res_'));
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveto = fullfile(func_dir,'res');
            gi = gi + 1;
            
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = fullfile(fullfile(func_dir,'tmp'),'spmF_0002.nii') ;
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = fullfile(func_dir,'res');
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = 'spmF_0002';
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = char(spm_file(spm_file(epi,'basename'),'prefix','Physio_'));
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
            gi = gi + 1;
            
            matlabbatch{gi,sub}.cfg_basicio.file_dir.dir_ops.dir_move.dir = fullfile(func_dir,'tmp');%delete tmp dir
            matlabbatch{gi,sub}.cfg_basicio.file_dir.dir_ops.dir_move.action.delete = true;
            gi = gi + 1;
            
            if rename_mat == 1 %assume original realign
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = spm_file(epi,'prefix','brain_a','ext','.mat');
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = func_dir;
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = 'brain_asub';
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = 'asub';
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
                gi = gi + 1;
            end                     
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
    %run_spm_multiple(matlabbatch, n_procs);
else
    spm_figure('CreateWin', 'Graphics');
    save('clean','matlabbatch');
    spm_jobman('run',matlabbatch);
end

end
