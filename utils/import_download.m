function import_download(job)

% function import_download(job)
% worker for the SPM batch system (cfg_basicio call_matlab) used by sbp_import_data:
% runs the SFTP batch file of ONE session and checks that DICOM files arrived.
% The (very long) transcript of the SFTP program is not shown, it is saved in
% <ftp_file without .txt>.log, and only a short summary is printed.
%   job.machina   'lin_ser', 'win_loc' or 'lin_loc'
%   job.conn      struct with user and server
%   job.pw        password (empty on the server)
%   job.ftp_file  SFTP batch file (lcd / cd / mget MR.* for each series)
%   job.dirs      cellstr, local temp dirs that should contain MR.* files afterwards

[~,name] = fileparts(job.ftp_file);
fprintf('download %s: started (%d series)\n',name,numel(job.dirs));
t0  = tic;
out = sftp_cmd(job.machina,job.conn,job.pw,job.ftp_file,true);
log_file = regexprep(job.ftp_file,'\.txt$','.log');
fid = fopen(log_file,'w');
if fid > 0
    fprintf(fid,'%s',out);
    fclose(fid);
end
n_files = 0;
for i = 1:numel(job.dirs)
    n_i = size(spm_select('List',job.dirs{i},'^MR.*'),1);
    if n_i == 0
        warning('import_download: no DICOM files (MR.*) in %s after the download (see %s)',job.dirs{i},log_file);
    end
    n_files = n_files + n_i;
end
fprintf('download %s: done, %d files in %.0f s\n',name,n_files,toc(t0));
end
