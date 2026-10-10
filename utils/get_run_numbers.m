function runs = get_run_numbers(BIDS,sub,ses,varargin)

% function runs = get_run_numbers(BIDS,sub,ses,varargin)
% run numbers (row vector, ascending) of the functional images that EXIST for
% a subject and session, e.g. for  for run = get_run_numbers(BIDS,sub_id,ses)
% instead of  for run = 1:vars.nRuns  (subjects can have fewer runs).
% Same defaults and options as get_func_files (brain, task = vars.task, ...)
% Returns [] if there are none. Which runs are used in the analysis is
% decided there (see analysis.exclude), all existing runs are preprocessed.

files = get_func_files(BIDS,sub,ses,[],varargin{:},'expect',NaN);
runs  = zeros(1,0);
for ii = 1:numel(files)
    t = regexp(spm_file(files{ii},'filename'),'_run-(\d+)','tokens','once');
    if ~isempty(t)
        runs(end+1) = str2double(t{1}); %#ok<AGROW>
    end
end
runs = unique(runs);
end
