function apply_vdm(all_sub_ids)
% function realign(all_sub_ids)
% realign and reslice epi images
% This is part 2 of 2 here we take a brain mask and actually reslice all images
% Applies the vdm (vdm5_<P><epi>) and the realignment to the P<epi> images (P: see epi_prefix)
% and writes the resliced images r<P><epi> (fixed prefix r, used by sbp_create_noise and the analysis)

specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);

    gi  = 1;
    cnt = 1;
    for ses = 1:vars.nSess
        for run = get_run_numbers(BIDS,sub_id,ses)
            epi                =  get_func_files(BIDS,sub_id,ses,run);
            run_niftis         =  epi_frames(epi,vars,'in'); % frames of P<epi>
            matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.data(cnt).scans   = cellstr(run_niftis);
            vdm_file = epi_file(epi,vars,'vdm'); % vdm5_<P><epi>
            if ~isfile(char(vdm_file))
                error('sub-%02d ses-%02d run-%02d: vdm not found (%s), was sbp_create_vdm_fmap run?',sub_id,ses,run,char(vdm_file));
            end
            matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.data(cnt).vdmfile = cellstr(vdm_file);
            cnt                =  cnt + 1;
        end
    end
    
    matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.roptions.pedir = 2; % Y direction
    matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.roptions.which = [2 1];
    matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.roptions.rinterp = 4;
    matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.roptions.wrap = [0 0 0];
    matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.roptions.mask = 1;
    matlabbatch{gi,sub}.spm.tools.fieldmap.applyvdm.roptions.prefix = 'r';
end

% save('apply_vdm','matlabbatch');
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
