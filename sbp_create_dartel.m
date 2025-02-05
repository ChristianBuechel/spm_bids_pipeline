function create_dartel(all_sub_ids)
% function create_dartel
% estimates the dartel flowfields from the T1 to the Template


% resolve paths
[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);

template_t   = vars.templateID;

n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    
    struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'type','T1w');
    struc_file = struc_file(1); % if there are more than 1
    anat_dir   = spm_file(struc_file,'path');

    rc1_file     = spm_file(struc_file,'prefix','rc1');    
    rc2_file     = spm_file(struc_file,'prefix','rc2');    

    gi    = 1;
    %matlabbatch{gi,sub}.spm.tools.dartel.warp1.images = {cellstr(rc1_file{1}),cellstr(rc2_file{1})};
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.images = {rc1_file rc2_file}';
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.rform = 0;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(1).its = 3;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(1).K = 0;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(1).template = {[path.templateDir sprintf(template_t,1)]};
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(2).its = 3;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(2).K = 0;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(2).template = {[path.templateDir sprintf(template_t,2)]};
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(3).its = 3;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(3).K = 1;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(3).template = {[path.templateDir sprintf(template_t,3)]};
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(4).its = 3;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(4).K = 2;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(4).template = {[path.templateDir sprintf(template_t,4)]};
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(5).its = 3;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(5).K = 4;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(5).template = {[path.templateDir sprintf(template_t,5)]};
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(6).its = 3;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(6).K = 6;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.param(6).template = {[path.templateDir sprintf(template_t,6)]};
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
    matlabbatch{gi,sub}.spm.tools.dartel.warp1.settings.optim.its = 3;
    gi = gi + 1;
    
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
