function lin_coreg(all_sub_ids)
% function lin_coreg
% performs a linear coregistration between umean EPI (ie epi fieldmap corrected) and T1
% changes the meana* file to be in the space of the T1
% it then applies this transformation to all epis and all vdms that were already
% realigned once


% resolve paths
[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    gi  = 1;
    sub_id     = all_sub_ids(sub);
    
    struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'type','T1w');
    struc_file = struc_file(1); %in case there are more sessions with T1w
    anat_dir   = spm_file(struc_file,'path');
    
    epi        = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    epi        = epi(1); % if there are brain and spinal, just take brain
    mean_file  = spm_file(epi,'prefix','umeana');    % we need the field map corrected mean umeanasub 
    mean_dir   = spm_file(mean_file,'path');

    all_niftis = [];
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            epi                =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            epi                =  epi(1); % if there are brain and spinal, just take brain
            run_niftis         =  spm_select('ExtFPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^a'),Inf); % all epis
            all_niftis         =  strvcat(all_niftis,run_niftis);
            run_vdms           =  spm_select('FPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^vdm5_a')); % all vdms
            all_niftis         =  strvcat(all_niftis,run_vdms);            
        end
    end    
    matlabbatch{gi,sub}.spm.spatial.coreg.estimate.ref = struc_file;
    matlabbatch{gi,sub}.spm.spatial.coreg.estimate.source = mean_file;
    matlabbatch{gi,sub}.spm.spatial.coreg.estimate.other = cellstr(all_niftis);
    matlabbatch{gi,sub}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{gi,sub}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{gi,sub}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{gi,sub}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    gi = gi + 1;
    
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = mean_file; % rename umeana file to cmeana now that it is coregistered
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = mean_dir; % same dir
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(1).pattern = 'umeana';
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(1).repl = 'cmeana';
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
    gi = gi + 1;

end

% run matlabbatch
n_procs = n_subs; % to not block too many cores on the server

if n_procs > vars.max_procs
    n_procs = vars.max_procs;
end
save('lin_coreg','matlabbatch')
if vars.parallel == 1
    run_spm_parallel(matlabbatch, n_procs);
    %run_spm_multiple(matlabbatch, n_procs);
else
    run_spm_sequential(matlabbatch);
end
