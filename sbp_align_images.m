function sbp_align_images(all_sub_ids)

% function sbp_align_images(all_sub_ids)
% does a quick (rigid body) coregistration of the T1 image in each session with a
% skullstripped template and then moves ALL images of that session by the same
% shift. This makes it easier to look at single subject data as they are (close)
% to template space. It is a coarse alignment, movement within a session is not
% corrected.
%
% The T1 used is T1w run-01 of the session. If you average several T1s
% (sbp_average_T1), do that BEFORE, then run-01 is the average.
%
% What is moved: every BIDS named nifti (*.nii whose name starts with sub-) in
% the session folder, i.e. anat (all T1w runs), func (brain AND spinal, all tasks,
% mag and phase) and fmap (magnitude, phasediff, epi ...), and anything new that
% follows this naming. This works for spinal studies too (the T1 covers brain and
% spinal cord): everything moves, so everything stays registered.
% What is NOT moved: images created by later processing steps (their names do not
% start with sub-), the gzip copies in path.baseDir (raw data stay as acquired),
% images that were already realigned (they have a .mat file) and BIDS named images that
% were computed from aligned images (descrip starts with "derived", e.g. the T1-like
% images of sbp_vaso_T1).
%
% Images are moved by changing mat and the ORIGINAL position stays in mat0, so
% running the function again is safe: moved images are recognised (mat ~= mat0) and
% skipped, and images that were added to an already aligned session later are
% moved by the shift of the T1.
% Only works if you have a T1 in EACH session, and the T1 has to contain the brain
% (a T1 or T2 that only covers the spinal cord cannot be aligned to the brain template).
% The coregistration starts from the scanner position. It worked for all test T1s
% (brain only and brain + spinal cord) whose field of view centre was within about 80 mm
% of the template centre; a warning is given if the distance is larger than 70 mm.

% resolve paths
specs    = get_study_specs;
path     = specs.path;
vars     = specs.vars;

BIDS     = spm_BIDS(path.preprocDir);
n_subs   = length(all_sub_ids);

smoothing      = 6;
max_start_dist = 70; % mm, see below

% load and smooth template once
Vtemplate         =  spm_vol(fullfile(path.templateDir,vars.templateStripID));
Vtemplate.uint8   =  loaduint8(Vtemplate);
voxsize           =  sqrt(sum(Vtemplate.mat(1:3,1:3).^2, 1)); % in mm
Vtemplate         =  smooth_uint8(Vtemplate,smoothing./voxsize); % smoothing

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    for ses = 1:vars.nSess
        ses_dir    = fullfile(path.preprocDir,sprintf('sub-%02d',sub_id),sprintf('ses-%02d',ses));
        struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'type','T1w','run','01');
        if isempty(struc_file)
            error('sub-%02d ses-%02d: no T1w run-01 found, cannot align',sub_id,ses);
        end
        struc_file = char(struc_file(1));
        N_T1       = nifti(struc_file);

        if is_moved(N_T1)
            % T1 was moved before: apply its shift to images that are not moved yet
            shift_mat = N_T1.mat/N_T1.mat0;
            fprintf('sub-%02d ses-%02d: T1 is already aligned, applying its shift to the remaining images\n',sub_id,ses);
        else
            VT1               =  spm_vol(struc_file);
            % header only check: the coarse coregistration starts from the scanner position and
            % finds the template if the images start within roughly 60-80 mm of each other
            c_T1  = VT1.mat*[VT1.dim(1:3)/2 1]';
            c_tpl = Vtemplate.mat*[Vtemplate.dim(1:3)/2 1]';
            d_fov = norm(c_T1(1:3)-c_tpl(1:3));
            if d_fov > max_start_dist
                warning(['sub-%02d ses-%02d: the centre of the T1 field of view is %.0f mm from the template centre, the coarse coregistration ' ...
                    'may fail. Check the result (e.g. spm_check_registration) and, if needed, set the origin of the T1 by hand (SPM Display) first.'],sub_id,ses,d_fov);
            end
            VT1.uint8         =  loaduint8(VT1);
            voxsize           =  sqrt(sum(VT1.mat(1:3,1:3).^2, 1)); % in mm
            VT1               =  smooth_uint8(VT1,smoothing./voxsize); % smoothing

            flags              = struct('cost_fun','nmi','sep',[4 2]); % both are T1 weighted, but ncc is worse than nmi!
            x                  =  spm_coreg(Vtemplate,VT1,flags);
            shift_mat          =  inv(spm_matrix(x(:)'));
            P                  =  spm_imatrix(shift_mat);
            fprintf('sub-%02d ses-%02d: shift %.1f %.1f %.1f mm, rotation %.1f %.1f %.1f deg\n',sub_id,ses,P(1:3),P(4:6)*180/pi);
        end

        images  = session_images(ses_dir);
        n_moved = 0; skipped = {};
        for i = 1:numel(images)
            if move_image(images{i},shift_mat)
                n_moved = n_moved + 1;
            else
                skipped{end+1} = spm_file(images{i},'filename'); %#ok<AGROW>
            end
        end
        fprintf('   %d images moved, %d skipped (already moved or realigned)\n',n_moved,numel(skipped));
        for i = 1:numel(skipped)
            fprintf('      skipped: %s\n',skipped{i});
        end
    end
end
end

%% ------------------------------------------------------------------------
function files = session_images(ses_dir)
% all BIDS named niftis (sub-*.nii) below the session folder
d     = dir(fullfile(ses_dir,'**','*.nii'));
d     = d(~[d.isdir]);
d     = d(startsWith({d.name},'sub-'));
files = fullfile({d.folder},{d.name})';
end

function moved = is_moved(N)
% images that were moved have mat ~= mat0 (mat0 keeps the original position)
moved = sum(sum((N.mat-N.mat0).^2)) > 1e-3;
end

function did_move = move_image(in_image,shift_mat)
% change mat with create (as for the 4D EPIs, this does not need a .mat file and
% keeps the original position in mat0)
N = nifti(in_image);
if isfield(N.extras,'mat')
    fprintf('   %s has a .mat file, seems to be realigned already - SKIP\n',in_image);
    did_move = false;
elseif startsWith(strtrim(N.descrip),'derived')
    % images computed from images that are already aligned (e.g. sbp_vaso_T1) are in the aligned space
    did_move = false;
elseif is_moved(N)
    did_move = false;
else
    N.mat = shift_mat*N.mat;
    create(N);
    did_move = true;
end
end
