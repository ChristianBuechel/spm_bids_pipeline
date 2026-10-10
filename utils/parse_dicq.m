function series = parse_dicq(txt)

% function series = parse_dicq(txt)
% parses the output of 'dicq -f --series --exam=...' (char or cellstr, one
% element per line is fine) into a struct array with the fields
%   ser    series number
%   name   series description exactly as listed (incl. padding blanks, so that
%          import.data.seq entries with trailing blank match as before)
%   n      number of images
%   path   DICOM directory on the server

if iscell(txt)
    lines = txt(:);
else
    lines = strsplit(char(txt),{'\r\n','\n'},'CollapseDelimiters',false)';
end
rx = '^\s*-\s*Series:\s*(\d+)\s*\{(.*)\}\s*\[[^\]]*\|[^\]]*\|\s*(\d+)\s*\]\s*(\S*)';
series = struct('ser',{},'name',{},'n',{},'path',{});
for ii = 1:numel(lines)
    t = regexp(lines{ii},rx,'tokens','once');
    if isempty(t), continue; end
    series(end+1) = struct('ser',str2double(t{1}),'name',t{2},'n',str2double(t{3}),'path',t{4}); %#ok<AGROW>
end
end
