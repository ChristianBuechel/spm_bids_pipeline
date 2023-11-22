function realign_bs(all_sub_ids)
% function realign(all_sub_ids)
% realign and reslice epi images
% Gets the filename of the raw epi files from get_base_dir

[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);
run_parallel = 1;

% TO DO

% estimate motion and apply to asub* images (could get rid of hi freq oscillations --> yes!)

% better: if physio exists --> residualize EPIs and save
% do all religenments with them but copy *.mat

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    
    epi        =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    
    cnt        = 1;
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            epi =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            if exist(char(spm_file(epi,'prefix','a','ext','.mat')),'file') %assume original realign
                movefile(char(spm_file(epi,'prefix','a','ext','.mat')),char(spm_file(epi,'prefix','brain_a','ext','.mat')));
            end
            if exist(char(spm_file(epi,'prefix','rp_a','ext','.txt')),'file') %assume original realign
                movefile(char(spm_file(epi,'prefix','rp_a','ext','.txt')),char(spm_file(epi,'prefix','brain_rp_a','ext','.txt')))
            end            
            if exist(char(spm_file(epi,'prefix','meana')),'file') %assume original realign
                movefile(char(spm_file(epi,'prefix','meana')),char(spm_file(epi,'prefix','brain_meana')));
            end
            
            run_niftis         =  spm_select('ExtFPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^a'),Inf);
            all_niftis{cnt,1}  =  cellstr(run_niftis);
            cnt                =  cnt + 1;
        end
    end
    func_dir   = spm_file(epi,'path');
    % fill with default options
    gi  = 1;    
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.data             = all_niftis;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.sep     = 3;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.fwhm    = 4;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.rtm     = 1;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.interp  = 3;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.weight  = fullfile(func_dir,'bins3cbrainstem_mask.nii');
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.which   = [2 1]; %all images, no mean
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.interp  = 4;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.mask    = 1; % if set to 1 takes out people who move out of FOV
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.prefix  = 'br';
    gi = gi + 1;
    
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            epi =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = spm_file(epi,'prefix','rp_a','ext','.txt'); %rename rp*txt files to brp*txt
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = func_dir;
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = 'rp_asub';
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = 'rp_basub';
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
            gi = gi + 1;
            
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = spm_file(epi,'prefix','a','ext','.mat'); %rename rp*txt files to brp*txt
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = func_dir;
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = 'asub';
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = 'basub';
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
            gi = gi + 1;
            
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = spm_file(epi,'prefix','meana'); %rename mean to bs_mean
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = func_dir;
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = 'meana';
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = 'bs_meana';
            matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
            gi = gi + 1;
                        
            if exist(char(spm_file(epi,'prefix','brain_a','ext','.mat')),'file') % restore original rp*txt files
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = spm_file(epi,'prefix','brain_a','ext','.mat');
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = func_dir;
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = 'brain_a';
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = 'a';
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
                gi = gi + 1;
            end
                        
            if exist(char(spm_file(epi,'prefix','brain_rp_a','ext','.txt')),'file') % restore original rp*txt files
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = spm_file(epi,'prefix','brain_rp_a','ext','.txt');
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = func_dir;
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = 'brain_rp_a';
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = 'rp_a';
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
                gi = gi + 1;
            end
            
            if exist(char(spm_file(epi,'prefix','brain_meana')),'file') % restore original mean file
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = spm_file(epi,'prefix','brain_meana');
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = func_dir;
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = 'brain_meana';
                matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = 'meana';
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
    save('realign_bs','matlabbatch');
    spm_jobman('run',matlabbatch);
end

end
