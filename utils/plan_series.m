function [planned,skipped] = plan_series(series,data,skip_ser)

% function [planned,skipped] = plan_series(series,data,skip_ser)
% decides which of the series of ONE session (see parse_dicq) are imported
%   series    struct array from parse_dicq
%   data      import.data from get_study_specs, fields dir, type, seq, cond and
%             optionally task (task label if not the default task)
%   skip_ser  series numbers to leave out (import.prisma_skip), optional
% A series is taken by the FIRST entry of data for which its name contains seq
% and cond (a string using n = number of images and ser = series number) is
% true. Series that match a seq but fail every cond, or are on the skip list,
% are returned in skipped (with the reason) so that nothing disappears silently.
%
% planned has the fields ser, name, n, path, mod (= dir), type, task, entry
% and an empty subtype (to be filled from the DICOM headers, then assign_runs)

if nargin < 3 || isempty(skip_ser) || (iscell(skip_ser) && isempty(skip_ser))
    skip_ser = [];
end
if iscell(skip_ser), skip_ser = [skip_ser{:}]; end

planned = struct('ser',{},'name',{},'n',{},'path',{},'mod',{},'type',{},'task',{},'entry',{},'subtype',{});
skipped = struct('ser',{},'name',{},'n',{},'reason',{});

for ii = 1:numel(series)
    s       = series(ii);
    matched = find(arrayfun(@(d) contains(s.name,d.seq),data));
    if isempty(matched), continue; end
    if ismember(s.ser,skip_ser)
        skipped(end+1) = struct('ser',s.ser,'name',s.name,'n',s.n,'reason','on skip list'); %#ok<AGROW>
        continue
    end
    take = 0;
    for tt = matched
        if eval_cond(data(tt).cond,s.n,s.ser)
            take = tt; break
        end
    end
    if take == 0
        skipped(end+1) = struct('ser',s.ser,'name',s.name,'n',s.n, ...
            'reason',sprintf('matches seq of entry %d but not cond (%s), n = %d',matched(1),data(matched(1)).cond,s.n)); %#ok<AGROW>
        continue
    end
    task = '';
    if isfield(data,'task') && ~isempty(data(take).task), task = data(take).task; end
    planned(end+1) = struct('ser',s.ser,'name',s.name,'n',s.n,'path',s.path, ...
        'mod',data(take).dir,'type',data(take).type,'task',task,'entry',take,'subtype',''); %#ok<AGROW>
end
end

function ok = eval_cond(cond,n,ser) %#ok<INUSD>
if isempty(cond), ok = true; return; end
ok = logical(eval(cond));
end
