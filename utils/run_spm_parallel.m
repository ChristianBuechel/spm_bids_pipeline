function run_spm_parallel(matlabbatch,n_cores)
% basically lukas function to parallelize SPM preprocessing
% relies on parallel computing toolbox
% start parallel pool if there is none, splot matlabbatch and run parallel

if ~license('test', 'Distrib_Computing_Toolbox')
    run_spm_multiple(matlabbatch,n_cores); % starts multiple matlabs via system calls
else
    start_pool_conditionally(n_cores);
    n_cores = min(n_cores,size(matlabbatch,2));
    loop_procs = split_vect(matlabbatch,n_cores);
    parfor worker = 1:n_cores
        matlabbatch = vec(loop_procs{worker})';
        %save(sprintf('%d_test.mat',worker),'matlabbatch');
        do_one_batch(matlabbatch);
    end
end
end
function do_one_batch(input)
input = input(~cellfun(@isempty,input)); % subjects can have different numbers of jobs --> empty cells, which spm_jobman does not accept
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