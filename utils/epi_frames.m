function frames = epi_frames(epi,vars,kind,nframes)

% function frames = epi_frames(epi,vars,kind,nframes)
% the frames (<file>,1 <file>,2 ...) of the 4D EPI image of a stage, as char matrix
% as returned by spm_select('ExtFPlist',...). Stops with a clear message if there are
% none, which is the case if an earlier step has not been run.
%   epi      raw EPI (char or 1x1 cellstr)
%   kind     'in' (default, input of realign etc.), 'res' (resliced), 'slicetiming_in' or
%            'nordic_in' (input of those steps), see epi_file
%   nframes  number of frames, default Inf (all)

if nargin < 3 || isempty(kind),    kind    = 'in'; end
if nargin < 4 || isempty(nframes), nframes = Inf;  end
epi    = char(epi);
target = char(epi_file(epi,vars,kind));
frames = spm_select('ExtFPlist',spm_file(target,'path'),['^' regexptranslate('escape',spm_file(target,'filename')) '$'],nframes);
if isempty(frames)
    steps = '{}';
    if isfield(vars,'preproc') && isfield(vars.preproc,'steps') && ~isempty(vars.preproc.steps)
        steps = ['{' strjoin(cellstr(vars.preproc.steps),', ') '}'];
    end
    error(['EPI images not found: %s\nThey are expected with the prefix "%s" (vars.preproc.steps = %s). ' ...
        'Was the step that creates them (sbp_nordic, sbp_slice_timing, sbp_realign_2_2, sbp_apply_vdm, sbp_smooth, ...) run, and does vars.preproc.steps say what was run?'], ...
        target,strrep(spm_file(target,'filename'),spm_file(epi,'filename'),''),steps);
end
end
