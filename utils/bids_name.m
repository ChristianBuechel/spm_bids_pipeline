function fname = bids_name(sub,ses,suffix,varargin)

% function fname = bids_name(sub,ses,suffix,varargin)
% builds a (light weight) BIDS file name, e.g.
%   bids_name(3,1,'bold','task','pain','run',2,'acq','brain')
%   --> sub-03_ses-01_task-pain_run-02_acq-brain_bold.nii
%
% sub, ses  numbers (ses = [] for no session entity)
% suffix    e.g. 'bold', 'cbv', 'T1w', 'epi', 'magnitude1'
% name-value pairs (empty values are skipped):
%   'task' 'dir' 'run' 'acq' 'part' (run is a number, the rest strings)
%   'ext'  extension incl. dot, default '.nii'
% Entity order is the one used throughout this pipeline:
%   sub, ses, task, dir, run, acq, part

order = {'task','dir','run','acq','part'};
opt   = struct('task','','dir','','run',[],'acq','','part','','ext','.nii');

if mod(numel(varargin),2)
    error('Name-value pairs expected');
end
for ii = 1:2:numel(varargin)
    if ~isfield(opt,varargin{ii})
        error('Unknown entity %s',varargin{ii});
    end
    opt.(varargin{ii}) = varargin{ii+1};
end

parts = {sprintf('sub-%02d',sub)};
if ~isempty(ses)
    parts{end+1} = sprintf('ses-%02d',ses);
end
for ii = 1:numel(order)
    val = opt.(order{ii});
    if isempty(val), continue; end
    if strcmp(order{ii},'run')
        parts{end+1} = sprintf('run-%02d',val); %#ok<AGROW>
    else
        parts{end+1} = sprintf('%s-%s',order{ii},val); %#ok<AGROW>
    end
end
parts{end+1} = suffix;
fname = [strjoin(parts,'_') opt.ext];
end
