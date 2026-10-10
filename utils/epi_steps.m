function steps = epi_steps(vars)

% function steps = epi_steps(vars)
% the preprocessing steps that were run on the EPIs, in the order they were run, as cellstr
% (vars.preproc.steps from get_study_specs, default {'slicetiming'} if it is missing).
% Valid steps (exactly these names): 'nordic', 'slicetiming', 'realign', 'smooth'.
%   before 'realign': nordic, slicetiming (any order)      after 'realign': smooth
% 'realign' is always part of the returned list: if it is not listed, it is put after the steps
% that are listed (before 'smooth' if that is listed), so {'nordic'} and {'nordic','realign'} are
% the same. Examples (returned list -> prefix of the images used for the analysis, see epi_prefix):
%   {'nordic','realign','smooth'} -> 'srn',  {'slicetiming','realign'} -> 'ra',  {'nordic','slicetiming','realign'} -> 'ran'
% See epi_prefix.

if isfield(vars,'preproc') && isstruct(vars.preproc) && isfield(vars.preproc,'steps')
    steps = vars.preproc.steps;
    if isempty(steps), steps = {}; else, steps = cellstr(steps); end
else
    steps = {'slicetiming'};
end
steps = lower(strtrim(steps(:)'));

valid = {'nordic','slicetiming','realign','smooth'};
bad   = steps(~ismember(steps,valid));
if ~isempty(bad)
    error('vars.preproc.steps: unknown step "%s" (valid: %s)',bad{1},strjoin(valid,', '));
end
if numel(unique(steps)) ~= numel(steps)
    error('vars.preproc.steps: a step is listed twice');
end

if ~ismember('realign',steps)                 % realignment is implicit
    k = find(strcmp(steps,'smooth'),1);
    if isempty(k), steps{end+1} = 'realign'; else, steps = [steps(1:k-1) {'realign'} steps(k:end)]; end
end

kr = find(strcmp(steps,'realign'));            % order: nordic / slicetiming  <  realign  <  smooth
for s = steps(1:kr-1)
    if ~ismember(s{1},{'nordic','slicetiming'})
        error('vars.preproc.steps: "%s" has to come after "realign"',s{1});
    end
end
for s = steps(kr+1:end)
    if ~ismember(s{1},{'smooth'})
        error('vars.preproc.steps: "%s" has to come before "realign"',s{1});
    end
end
end
