function [status,result] = remote_cmd(machina,conn,pw,cmd)

% function [status,result] = remote_cmd(machina,conn,pw,cmd)
% run a command on the DICOM server (used by sbp_import_data)
%   machina  'lin_ser' (we are on the server), 'win_loc' or 'lin_loc' (see sbp_import_data)
%   conn     struct with the fields user and server
%   pw       password (not needed on the server)

switch machina
    case 'lin_ser'
        [status, result] = system(sprintf('ssh %s@%s %s',conn.user,conn.server,cmd));
    case 'win_loc'
        [status, result] = system(sprintf('plink -ssh %s@%s -pw %s %s',conn.user,conn.server,pw,cmd));
    case 'lin_loc'
        [status, result] = system(sprintf('SSHPASS=%s sshpass -e ssh %s@%s %s',pw,conn.user,conn.server,cmd));
    otherwise
        error('unknown machine type %s',machina);
end
end
