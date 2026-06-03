function sbp_multiband_correct(all_sub_ids)
% function sbp_nordic
% does NORDIC noise filtering on data in image space

% add paths
[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id    = all_sub_ids(sub);
       
    epi        = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    epi        = epi(1); % if there are brain and spinal, just take brain
    
    func_dir   = spm_file(epi,'path');
    
    gi = 1;
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            fprintf('ses %d run %d\n',ses,run);
            epi       =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            epi       =  epi(1); % if there are brain and spinal, just take brain
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.string    = char(epi);
            
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.inputs{3}.evaluated = arg;
            
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
            matlabbatch{gi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'multiband_correct';
            gi    = gi + 1;
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
