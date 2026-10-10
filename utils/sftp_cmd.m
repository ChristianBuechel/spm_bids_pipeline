function result = sftp_cmd(machina,conn,pw,ftp_file,quiet)

% function result = sftp_cmd(machina,conn,pw,ftp_file,quiet)
% run an SFTP batch file against the DICOM server (used by sbp_import_data)
% the output is shown unless quiet is true, and it is returned in result
%   machina  'lin_ser', 'win_loc' or 'lin_loc' (see sbp_import_data)
%   conn     struct with the fields user and server

if nargin < 5, quiet = false; end
if quiet, echo_flag = {}; else, echo_flag = {'-echo'}; end
switch machina
    case 'lin_ser'
        [~,result] = system(sprintf('sftp -b %s %s@%s',ftp_file,conn.user,conn.server),echo_flag{:});
    case 'win_loc'
        [~,result] = system(sprintf('psftp %s@%s -pw %s -b %s',conn.user,conn.server,pw,ftp_file),echo_flag{:});
    case 'lin_loc'
        [~,result] = system(sprintf('SSHPASS=%s sshpass -e sftp -o BatchMode=no -b %s -o PubkeyAuthentication=no %s@%s',pw,ftp_file,conn.user,conn.server),echo_flag{:});
    otherwise
        error('unknown machine type %s',machina);
end
end
