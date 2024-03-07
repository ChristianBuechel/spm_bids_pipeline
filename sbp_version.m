function [status, infos] = sbp_version()
%
% return infos about the current version of the toolbox. It requires the
% toolbox to be installed using git (not by downlading the zip) to work properly:
%       git clone https://github.com/ChristianBuechel/spm_bids_pipeline.git
% OUT:
%   - infos.path: path of the toolbox in use
%   * under git versioning, the function also returns:
%   - infos.version: the current branch and revision of the toolbox (branch/commit)
%   - infos.status: a summary of the status of the current revision compared to
%              the one found on github.
% adopted from VBA-toolbox VBA_version.m

str = which('sbp_version');
infos.path = fileparts(str);
cd (infos.path);
try % if the toolbox was installed through git
    gitInf = getGitInfo ;
    infos.version = [gitInf.branch '/' gitInf.hash];
    infos.git = true;
    
    try % try to see if new commits are online
        request = sprintf('https://api.github.com/repos/ChristianBuechel/spm_bids_pipeline/compare/%s...%s',gitInf.branch,gitInf.hash);
        tracker=webread(request);
        
        switch tracker.status
            case 'identical'
                status = 'The toolbox is up to date.';
            case 'behind'
                status = sprintf('The toolbox is %d revision(s) behind the online version.',tracker.behind_by);
            case 'ahead'
                status = sprintf('The toolbox is %d revision(s) ahead the online version.',tracker.ahead_by);
        end        
    end
    
catch
    infos.version = 'unkown (not under git control)';
    infos.git = false;
end

