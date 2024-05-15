function run_spm_sequential(matlabbatch)
% make sure volunteers are done sequentially if we do not want paralell
% computing

loop_procs = split_vect(matlabbatch,size(matlabbatch,2));

for job = 1:size(loop_procs,2)
    matlabbatch = vec(loop_procs{job})';
    do_one_batch(matlabbatch);
end
end
function do_one_batch(input)
spm('defaults','FMRI');
spm_jobman('initcfg');
spm_jobman('run', input);
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