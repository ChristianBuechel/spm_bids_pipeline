function realign_bs_physio(all_sub_ids)
% function realign(all_sub_ids)
% realign and reslice epi images otimized for the brainstem
% this version uses residualized images, ie those where the effect of
% physio noise has been removed by clean_physio.m

[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    
    epi        =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    func_dir   = spm_file(epi,'path');
    
    
    cnt        = 1;
    all_aniftis = [];
    gi  = 1;
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            epi =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            % copy all asub files to ./res
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = spm_file(epi,'prefix','a');
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.copyto = fullfile(func_dir,'res');
            gi = gi + 1;
            
            run_niftis         =  spm_select('ExtFPlist', fullfile(spm_file(epi,'path'),'res'), spm_file(spm_file(epi,'filename'),'prefix','^res_'),Inf);
            all_niftis{cnt,1}  =  cellstr(run_niftis);
            
            %select aniftis in func dir
            run_aniftis         =  spm_select('ExtFPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^a'),Inf);
            %change all to the ones in /res
            run_aniftis         =  spm_file(run_aniftis,'path', fullfile(spm_file(epi,'path'),'res'));
            all_aniftis         =  strvcat(all_aniftis,run_aniftis);            
            cnt                 =  cnt + 1;
        end
    end

    %estimate movement parameters for residualized EPIs
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.data = all_niftis; % res_sub*nii in ./res
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.sep = 3;
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.fwhm = 4;
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.rtm = 1;
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.interp = 3;
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
    matlabbatch{gi,sub}.spm.spatial.realign.estimate.eoptions.weight = fullfile(func_dir,'bins3cbrainstem_mask.nii');
    gi = gi + 1;
    
    %rename *.mat to asub*mat
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            epi =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = fullfile(func_dir,'res',spm_file(spm_file(epi,'filename'),'prefix','res_','ext','.mat')); %
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = fullfile(func_dir,'res');
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = 'res_';
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = 'a';
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
            gi = gi + 1;
        end
    end
    
    %and reslice asubs to pbrasub
    matlabbatch{gi,sub}.spm.spatial.realign.write.data = cellstr(all_aniftis);
    matlabbatch{gi,sub}.spm.spatial.realign.write.roptions.which = [2 1];
    matlabbatch{gi,sub}.spm.spatial.realign.write.roptions.interp = 4;
    matlabbatch{gi,sub}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
    matlabbatch{gi,sub}.spm.spatial.realign.write.roptions.mask = 1;
    matlabbatch{gi,sub}.spm.spatial.realign.write.roptions.prefix = 'pbr';
    gi = gi + 1;
    
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            epi =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');

            % now copy pbrasub files to fun dir (and the rp*txt as rp_pbasub)
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = fullfile(func_dir,'res',spm_file(spm_file(epi,'filename'),'prefix','pbra')); %
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = func_dir;
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = 'pbra';
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = 'pbra';
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
            gi = gi + 1;
            
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = fullfile(func_dir,'res',spm_file(spm_file(epi,'filename'),'prefix','rp_res_','ext','.txt')); %rename rp*txt files to brp*txt
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = func_dir;
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = 'rp_res';
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = 'rp_pbasub';
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
            gi = gi + 1;

            % delete the asub files in ./res            
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = fullfile(func_dir,'res',spm_file(spm_file(epi,'filename'),'prefix','a'));
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false; % 
            gi = gi + 1;
            
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
