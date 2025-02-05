function slice_timing(all_sub_ids)
% function slice_timing(all_sub_ids)
% corrects slice timing in epi images
% Gets the filename of the raw epi files from get_base_dir
% Also gets the timing of the scans form get_base_dir

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
            run_niftis         =  spm_select('ExtFPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^'),Inf);
            all_niftis{cnt,1}  =  cellstr(run_niftis);
            cnt                =  cnt + 1;
        end
    end
    % fill with default options
    gi  = 1;
    matlabbatch{gi,sub}.spm.temporal.st.scans    = all_niftis;
    matlabbatch{gi,sub}.spm.temporal.st.nslices  = vars.sliceTiming.nslices;
    matlabbatch{gi,sub}.spm.temporal.st.tr       = vars.sliceTiming.tr;
    matlabbatch{gi,sub}.spm.temporal.st.ta       = 0;
    matlabbatch{gi,sub}.spm.temporal.st.so       = vars.sliceTiming.so;
    matlabbatch{gi,sub}.spm.temporal.st.refslice = vars.sliceTiming.refslice;
    matlabbatch{gi,sub}.spm.temporal.st.prefix   = 'a';
    gi = gi + 1;    
end

% run matlabbatch
n_procs = n_subs; % to not block too many cores on the server

if n_procs > vars.max_procs
    n_procs = vars.max_procs;
end

if vars.parallel == 1
    run_spm_parallel(matlabbatch, n_procs);
    %run_spm_multiple(matlabbatch, n_procs);
else
    run_spm_sequential(matlabbatch);
end

