function write_protocol(file,sub_str,ses_str,exam_id,planned,skipped,out_files,timings,err_msg)

% function write_protocol(file,sub_str,ses_str,exam_id,planned,skipped,out_files,timings,err_msg)
% short text file written by sbp_import_data: which series became which BIDS
% files, which series were not imported (and why), and how long the steps took
%   out_files  cell, one cellstr of BIDS file names (text) for each series in planned
%   timings    (optional) struct with fields (seconds) batch1, batch2
%   err_msg    (optional) if not empty the session failed: written as ERROR line

fid = fopen(file,'w');
fprintf(fid,'Import protocol %s %s (%s), written %s\n\n',sub_str,ses_str,exam_id,char(datetime('now','Format','yyyy-MM-dd HH:mm')));
if nargin > 8 && ~isempty(err_msg)
    fprintf(fid,'ERROR: this session was NOT imported completely: %s\n\n',err_msg);
end
fprintf(fid,'series      n  DICOM series name -> BIDS file(s)\n');
for k = 1:numel(planned)
    files = {};
    if numel(out_files) >= k, files = out_files{k}; end
    if isempty(files), files = {'(nothing written)'}; end
    name = strtrim(planned(k).name);
    for j = 1:numel(files)
        if j == 1
            fprintf(fid,'%6d %6d  %s -> %s\n',planned(k).ser,planned(k).n,name,files{j});
        else
            fprintf(fid,'%6s %6s  %*s    %s\n','','',numel(name),'',files{j});
        end
    end
end
if ~isempty(skipped)
    fprintf(fid,'\nNOT imported:\n');
    for k = 1:numel(skipped)
        fprintf(fid,'%6d %6d  %s -> %s\n',skipped(k).ser,skipped(k).n,strtrim(skipped(k).name),skipped(k).reason);
    end
end
if nargin > 7 && ~isempty(timings)
    fprintf(fid,'\nTime [s] (whole call, all sessions): download + DICOM conversion %.0f, 4D files and copies %.0f\n',timings.batch1,timings.batch2);
end
fclose(fid);
end
