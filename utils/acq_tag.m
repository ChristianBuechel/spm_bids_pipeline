function tag = acq_tag(label,acq_type)

% function tag = acq_tag(label,acq_type)
% value of the BIDS acq- entity for an EPI volume/series labelled 'brain' or
% 'spinal' (see classify_epi_acq), given vars.acq_type of the study:
%   {'brain'}            --> ''       brain only studies carry NO acq- tag
%   {'spinal'}           --> 'spinal' spinal only studies are ALWAYS tagged
%                                     (they get an independent pipeline)
%   {'spinal','brain'}   --> 'spinal' or 'brain'
% The default acq filter of get_func_files ({'brain',''}) therefore selects
% brain data and never spinal data; use 'acq','spinal' for the latter.

if nargin < 2 || isempty(acq_type)
    acq_type = {'brain'};
end
acq_type = cellstr(acq_type);

if ~ismember(label,{'brain','spinal'})
    error('Unknown label %s (brain or spinal expected)',label);
end
if strcmp(label,'brain') && isequal(acq_type,{'brain'})
    tag = '';
else
    tag = label;
end
end
