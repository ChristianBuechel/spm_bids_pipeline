function realign_1_2(all_sub_ids)
% function realign_1_2(all_sub_ids)
% realign epi images (first pass)
% This is part 1 of 2 here we simply create a mean EPI for all the upcoming
% spatial things, but also update the .mat info
% Part 2 will then use a brain mask to get rid of the eyes

[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    cnt        = 1;
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            epi                =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            epi                =  epi(1); % if there are brain and spinal, just take brain
            run_niftis         =  spm_select('ExtFPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^a'),Inf);
            all_niftis{cnt,1}  =  cellstr(run_niftis);
            cnt                =  cnt + 1;
        end
    end
    % fill with default options
    gi  = 1;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.data             = all_niftis;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.quality = 0.95; % was 0.9 in SPM12
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.sep     = 1.5; % was 4 in SPM12
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.fwhm    = 1; % was 5 in SPM12;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.rtm     = 0; % one pass is ok
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.interp  = 2;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.eoptions.weight  = '';
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.which   = [0 1]; %only mean
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.interp  = 4;
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.mask    = 1; % if set to 1 takes out voxels from people who move out of FOV
    matlabbatch{gi,sub}.spm.spatial.realign.estwrite.roptions.prefix  = 'r';    
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
