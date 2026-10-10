function sbp_coreg_lin(all_sub_ids)
% function lin_coreg
% performs a linear coregistration between umean EPI (ie epi fieldmap corrected) and T1
% changes the mean<P>* file to be in the space of the T1 (cmean<P><epi>, a copy of umean<P><epi>)
% it then applies this transformation to all epis (frames of P<epi>) and all vdms (vdm5_<P><epi>)
% that were already realigned once (P: prefix of the realigned images, see epi_prefix)


% resolve paths
specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;

BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    gi  = 1;
    sub_id     = all_sub_ids(sub);
    
    struc_file = get_anat(BIDS,sub_id);
    anat_dir   = spm_file(struc_file,'path');
    
    epi        = get_ref_epi(BIDS,sub_id);
    mean_file  = epi_file(epi,vars,'umean');    % we need the field map corrected mean umean<P><epi> (P: see epi_prefix)
    cmean_file = epi_file(epi,vars,'cmean');    % ... which is copied to cmean<P><epi> and coregistered to the T1
    mean_dir   = spm_file(mean_file,'path');

    all_niftis = [];
    for ses = 1:vars.nSess
        for run = get_run_numbers(BIDS,sub_id,ses)
            epi                =  get_func_files(BIDS,sub_id,ses,run);
            run_niftis         =  epi_frames(epi,vars,'in'); % all epis (frames of P<epi>)
            all_niftis         =  strvcat(all_niftis,run_niftis);
            run_vdm            =  char(epi_file(epi,vars,'vdm')); % the vdm of this run (vdm5_<P><epi>)
            if ~isfile(run_vdm)
                error('sub-%02d ses-%02d run-%02d: vdm not found (%s), was sbp_create_vdm_fmap run?',sub_id,ses,run,run_vdm);
            end
            all_niftis         =  strvcat(all_niftis,run_vdm);
        end
    end


    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = mean_file;
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.copyto = mean_dir;
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.patrep = struct('pattern', {'^.*$'}, 'repl', {char(spm_file(cmean_file,'basename'))}); % umean<P><epi> --> cmean<P><epi>
    matlabbatch{gi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.unique = false;
    gi = gi + 1;

    matlabbatch{gi,sub}.spm.spatial.coreg.estimate.ref = struc_file;
    matlabbatch{gi,sub}.spm.spatial.coreg.estimate.source = cmean_file;
    matlabbatch{gi,sub}.spm.spatial.coreg.estimate.other = cellstr(all_niftis);
    matlabbatch{gi,sub}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{gi,sub}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{gi,sub}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{gi,sub}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    gi = gi + 1;
    
end

% run matlabbatch
n_procs = n_subs; % to not block too many cores on the server

if n_procs > vars.max_procs
    n_procs = vars.max_procs;
end
%save('lin_coreg','matlabbatch')
if vars.parallel == 1
    run_spm_parallel(matlabbatch, n_procs);
    %run_spm_multiple(matlabbatch, n_procs);
else
    run_spm_sequential(matlabbatch);
end
