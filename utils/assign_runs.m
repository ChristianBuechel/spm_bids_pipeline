function series = assign_runs(series)

% function series = assign_runs(series)
% assigns BIDS run numbers (and part- labels) to the series of ONE session
% No I/O, only bookkeeping, so it can be tested offline.
%
% input: struct array, one element per DICOM series, with fields
%   ser      series number on the scanner (defines acquisition order)
%   mod      'func', 'anat' or 'fmap'  (BIDS folder)
%   type     import type, e.g. 'bold', 'cbv', 'T1w', 'magnitude', 'phasediff'
%   n        number of images in the series
%   subtype  'mag', 'phase' or '' (from DICOM ImageType, '' if not applicable)
%   task     (optional) task label, e.g. 'vasoref'; '' for the default task.
%            Run numbers are counted per task, as BIDS defines them
% output: same array (sorted by ser) with the added fields
%   run      run number, counted separately for each mod (func, anat, fmap)
%            and task
%   part     'mag', 'phase' or ''
%
% Pairing: a magnitude series immediately followed (in acquisition order, within
% the same mod) by a phase series of the same type and the same number of
% images is ONE acquisition, i.e. both get the same run number and part-mag /
% part-phase. This replaces the post-hoc renaming of sbp_order_complex.
% A magnitude series without partner is a normal run (no part). A phase series
% without partner keeps part-phase and a warning is given.
% fmap: a 'magnitude' series immediately followed by a 'phasediff' series
% share a run number.

[~,o]  = sort([series.ser]);
series = series(o);
[series.run]  = deal(0);
[series.part] = deal('');
if ~isfield(series,'task')
    [series.task] = deal('');
end

% one counter per modality AND task
key  = strcat({series.mod},'|',{series.task});
keys = unique(key,'stable');
for m = 1:numel(keys)
    idx = find(strcmp(key,keys{m}));
    mod = series(idx(1)).mod;
    run = 0;
    k   = 1;
    while k <= numel(idx)
        i = idx(k);
        run = run + 1;
        series(i).run = run;
        paired = false;
        if k < numel(idx)
            j = idx(k+1);
            switch mod
                case 'func'
                    paired = strcmp(series(i).subtype,'mag') && strcmp(series(j).subtype,'phase') && ...
                        strcmp(series(i).type,series(j).type) && series(i).n == series(j).n;
                    if paired
                        series(i).part = 'mag';
                        series(j).part = 'phase';
                    end
                case 'fmap'
                    paired = strcmp(series(i).type,'magnitude') && strcmp(series(j).type,'phasediff');
            end
            if paired
                series(j).run = run;
            end
        end
        if ~paired && strcmp(series(i).subtype,'phase')
            warning('Series %d (%s) is a phase image without preceding magnitude image',series(i).ser,series(i).type);
            series(i).part = 'phase';
        end
        k = k + 1 + paired;
    end
end
end
