function files = get_anat(BIDS,sub,ses,run,varargin)

% function files = get_anat(BIDS,sub,ses,run,varargin)
% the T1w image of a subject (1x1 cell, like get_ref_epi for the EPI). By default
% session 1, run 1: that is the average of the T1s if sbp_average_T1 was run, and the
% image sbp_align_images uses. Run numbers are counted after the exclusions made at
% the import (participants.tsv), so run-01 exists if there is a T1 at all.
%   sub, ses, run  numbers (default ses = 1, run = 1)
%   'type'         image type (default 'T1w')
% Stops with a message if there is not exactly one such image.
% Derived files are named with fixed prefixes (c1, rc1, sc1, wc1, ...), e.g.
% spm_file(get_anat(BIDS,sub),'prefix','c1')

if nargin < 3 || isempty(ses), ses = 1; end
if nargin < 4 || isempty(run), run = 1; end
type = 'T1w';
for ii = 1:2:numel(varargin)
    if strcmp(varargin{ii},'type'), type = varargin{ii+1}; else, error('Unknown option %s',varargin{ii}); end
end

files = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'type',type);
if isempty(files), files = {}; end
files = cellstr(files);
if numel(files) ~= 1
    error('get_anat: expected 1 %s image for sub-%02d ses-%02d run-%02d, found %d',type,sub,ses,run,numel(files));
end
end
