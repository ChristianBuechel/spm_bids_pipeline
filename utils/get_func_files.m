function files = get_func_files(BIDS,sub,ses,run,varargin)

% function files = get_func_files(BIDS,sub,ses,run,varargin)
% returns functional images (cellstr) of a BIDS dataset (see spm_BIDS)
% sub, ses, run are numbers; use [] to not filter on ses and/or run
%
% Defaults (can be overridden by name-value pairs, values can be cellstr):
%   'task'   vars.task from get_study_specs
%   'acq'    {'brain',''}     brain only (acq-brain or no acq- tag), NOT acq-spinal
%   'part'   {'mag',''}       magnitude only (or no part- tag), NOT part-phase
%   'type'   {'bold','cbv'}
%   'expect' 1 if run is specified, else NaN (= no check)
%            number of files that are expected, error otherwise. Use NaN to switch off
%
% Examples:
%   epi  = get_func_files(BIDS,2,1,3);                        % sub 2, ses 1, run 3
%   epi  = get_func_files(BIDS,2,[],[],'expect',NaN);         % all func images of sub 2
%   cplx = get_func_files(BIDS,2,1,3,'part',{'mag','phase'},'expect',2);
%   spin = get_func_files(BIDS,2,1,3,'acq','spinal');

opt.task   = [];
opt.acq    = {'brain',''};
opt.part   = {'mag',''};
opt.type   = {'bold','cbv'};
opt.expect = [];

if mod(numel(varargin),2)
    error('Name-value pairs expected');
end
for ii = 1:2:numel(varargin)
    if ~isfield(opt,varargin{ii})
        error('Unknown option %s',varargin{ii});
    end
    opt.(varargin{ii}) = varargin{ii+1};
end

if isempty(opt.task)
    specs    = get_study_specs;
    opt.task = specs.vars.task;
end
if isempty(opt.expect)
    if isempty(run), opt.expect = NaN; else, opt.expect = 1; end
end

q = {'sub',sprintf('%02d',sub),'task',opt.task,'acq',opt.acq,'part',opt.part,'type',opt.type};
if ~isempty(ses), q = [q {'ses',sprintf('%02d',ses)}]; end
if ~isempty(run), q = [q {'run',sprintf('%02d',run)}]; end

files = spm_BIDS(BIDS,'data',q{:});
if isempty(files), files = {}; end
files = cellstr(files);

if ~isnan(opt.expect) && numel(files) ~= opt.expect
    error('get_func_files: expected %d file(s) but found %d for sub %d ses %s run %s',...
        opt.expect,numel(files),sub,mat2str(ses),mat2str(run));
end
end
