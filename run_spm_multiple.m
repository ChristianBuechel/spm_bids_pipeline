function run_spm_multiple(matlabbatches,n_cores)

loop_procs = split_vect(matlabbatches,n_cores);

spm_path          = fileparts(which('spm')); %get spm path
mat_name          = which(mfilename);
[~,mat_name,~]    = fileparts(mat_name);

for np = 1:n_cores
    fname = [mat_name '_' num2str(np) '.mat'];
    matlabbatch = vec(loop_procs{np})';
    save(fname,'matlabbatch');
    [~, name_stem]      = fileparts(fname);
    name_stem_m         = [name_stem '.m'];
    log_name            = strcat(name_stem, '.log');
    lo_cmd  = sprintf('clear matlabbatch;\nload(''%s'');\n',fname);
    ex_cmd  = sprintf('addpath(''%s'');\nspm(''defaults'',''FMRI'');\nspm_jobman(''initcfg'');\nspm_jobman(''run'',matlabbatch);\n',spm_path);
    end_cmd = sprintf('delete(''%s'');\ndelete(''%s'');\n',fname,name_stem_m);
    str     = strvcat(sprintf('function %s', name_stem),lo_cmd, ex_cmd, end_cmd, 'exit');
    spm_save(name_stem, str); %spm_save does not do .m files ...
    movefile(name_stem, name_stem_m); %rename
    
    if isunix
        matlab_exe  = [matlabroot filesep 'bin' filesep 'matlab']; % should be OK for unix
        cmd         = sprintf('xterm -e ''%s -nodesktop -nosplash -logfile %s -r "%s" '' &',matlab_exe, log_name, name_stem);
    elseif ispc
        matlab_exe  = 'matlab';
        cmd         = sprintf('start %s -nodesktop -nosplash -logfile %s -r "%s" exit',matlab_exe, log_name, name_stem);
    end
    
    system(cmd);
    
end


end

function vx = vec(X)
% transform array X into a column vector

% avoid weirdly sized empty vector
if isempty(X)
    vx = [];
    return;
end

% return vectorized input
vx = full(X(:));

end
