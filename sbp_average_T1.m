function sbp_average_T1(all_sub_ids)

% function sbp_average_T1(all_sub_ids)
% To improve the SNR of the T1w image: coregisters the first vars.nT1 T1w images
% (run-01 ... run-<nT1>) of each session to run-01 (estimate only), averages them
% and replaces run-01 by the average (the original run-01 is kept as *.bak).
% The other runs stay where they are (with their header coregistered to run-01).
% Run this BEFORE sbp_align_images (which uses T1w run-01 and moves all images).
% Safe to run twice: a session with a *.bak file is skipped.
%
% Across sessions (vars.average_T1.across_sessions = 1 in get_study_specs): if there is only one
% T1w per session (vars.nT1 = 1) but several sessions, T1w run-01 of all sessions are averaged instead:
% the T1w of the other sessions are coregistered to the one of session 1 (on temporary copies, the
% original images of the other sessions are NOT changed), averaged, and the average replaces
% T1w ses-01 run-01 (the original is kept as *.bak). This is the image that get_anat returns
% (ses-01 run-01). sbp_align_images still aligns each session with its own T1w.

% resolve paths
specs    = get_study_specs;
path     = specs.path;
vars     = specs.vars;

BIDS                  = spm_BIDS(path.preprocDir);
n_subs                = length(all_sub_ids);

across = isfield(vars,'average_T1') && isfield(vars.average_T1,'across_sessions') && vars.average_T1.across_sessions;
if across
    if vars.nSess < 2
        error('sbp_average_T1:nSess', 'vars.average_T1.across_sessions needs at least 2 sessions, but vars.nSess = %d.', vars.nSess);
    end
elseif vars.nT1 < 2
    error('sbp_average_T1:nT1', ['Need at least 2 T1w scans to average, but vars.nT1 = %d. With one T1w per session you can average ' ...
        'across the sessions: vars.average_T1.across_sessions = 1.'], vars.nT1);
end

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    if across
        average_across_sessions(BIDS,sub_id,vars.nSess);
        continue
    end
    for ses = 1:vars.nSess
        struc_files = cell(vars.nT1,1);
        for i_T1 = 1:vars.nT1        % should also work for more than 2 T1w scans
            % find all T1s, realign, create mean, rename 1st T1 and use mean as T1
            struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',i_T1),'type','T1w');
            if isempty(struc_file)
                error('sbp_average_T1:missingT1', 'sub-%02d ses-%02d: T1w run-%02d not found.', sub_id, ses, i_T1);
            end
            struc_files{i_T1,1} = char(struc_file(1));
        end

        % re-run guard: the .bak only exists if the 1st T1 was already replaced by the mean
        if exist(spm_file(struc_files{1},'ext','.bak'),'file')
            fprintf('sub-%02d ses-%02d: already averaged (.bak exists), skipping.\n', sub_id, ses);
            continue
        end
        % coregister the others to the 1st
        matlabbatch = [];
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {struc_files{1}};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {struc_files{2:end}};
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

        % now do mean
        matlabbatch{2}.spm.util.imcalc.input = struc_files;
        matlabbatch{2}.spm.util.imcalc.output = spm_file(struc_files{1},'prefix','mean');
        matlabbatch{2}.spm.util.imcalc.outdir = {spm_file(struc_files{1},'path')};
        matlabbatch{2}.spm.util.imcalc.expression = 'mean(X)';
        matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{2}.spm.util.imcalc.options.dmtx = 1;
        matlabbatch{2}.spm.util.imcalc.options.mask = 0;
        matlabbatch{2}.spm.util.imcalc.options.interp = -4;
        matlabbatch{2}.spm.util.imcalc.options.dtype = 512;

        spm_jobman('initcfg');
        spm('defaults', 'FMRI');
        spm_jobman('run',matlabbatch);

        % The coregistration changed the header of the other T1s (mat ~= mat0). Make that their
        % new reference position (mat0 = mat), otherwise sbp_align_images would take them for
        % images that are already aligned and would leave them behind.
        for i_T1 = 2:vars.nT1
            N      = nifti(struc_files{i_T1});
            N.mat0 = N.mat;
            create(N);
        end

        % rename first T1 file to *.bak
        movefile(struc_files{1}, spm_file(struc_files{1},'ext','.bak'));
        % rename mean to 1st T1
        movefile(spm_file(struc_files{1},'prefix','mean'),struc_files{1});
    end
end
end

%% ------------------------------------------------------------------------
function average_across_sessions(BIDS,sub_id,n_ses)
% T1w run-01 of all sessions: the ones of sessions 2..n are coregistered to session 1 (temporary
% copies), the average replaces T1w ses-01 run-01 (original kept as .bak)
files = cell(n_ses,1);
for ses = 1:n_ses
    f = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run','01','type','T1w');
    if isempty(f)
        error('sbp_average_T1:missingT1', 'sub-%02d ses-%02d: T1w run-01 not found.', sub_id, ses);
    end
    files{ses} = char(f(1));
end
if exist(spm_file(files{1},'ext','.bak'),'file')
    fprintf('sub-%02d: already averaged (.bak exists), skipping.\n', sub_id);
    return
end

% temporary copies of the other sessions (their header is changed by the coregistration, the originals stay as they are)
tmp = cell(n_ses-1,1);
for ses = 2:n_ses
    tmp{ses-1} = fullfile(spm_file(files{1},'path'),sprintf('tmpses%02d_%s',ses,spm_file(files{ses},'filename')));
    copyfile(files{ses},tmp{ses-1});
end

matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {files{1}};
matlabbatch{1}.spm.spatial.coreg.estimate.source = tmp;
matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

matlabbatch{2}.spm.util.imcalc.input = [files(1); tmp];
matlabbatch{2}.spm.util.imcalc.output = spm_file(files{1},'prefix','mean');
matlabbatch{2}.spm.util.imcalc.outdir = {spm_file(files{1},'path')};
matlabbatch{2}.spm.util.imcalc.expression = 'mean(X)';
matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{2}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{2}.spm.util.imcalc.options.mask = 0;
matlabbatch{2}.spm.util.imcalc.options.interp = -4;
matlabbatch{2}.spm.util.imcalc.options.dtype = 512;

spm_jobman('initcfg');
spm('defaults', 'FMRI');
spm_jobman('run',matlabbatch);

for ses = 1:numel(tmp)
    P = spm_imatrix(nifti(tmp{ses}).mat/nifti(files{ses+1}).mat);   % how far the session was moved to fit session 1
    fprintf('sub-%02d: ses-%02d T1w moved by %.1f %.1f %.1f mm, %.1f %.1f %.1f deg for the average\n',sub_id,ses+1,P(1:3),P(4:6)*180/pi);
    delete(tmp{ses});
end

movefile(files{1}, spm_file(files{1},'ext','.bak'));                % original of session 1
movefile(spm_file(files{1},'prefix','mean'),files{1});               % the average is T1w ses-01 run-01
N      = nifti(files{1});                                            % not an aligned image: mat0 = mat
N.mat0 = N.mat;
create(N);
end
