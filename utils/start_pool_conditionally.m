function [pool] = start_pool_conditionally(n_workers)
% check if parallel pool already runs if not start a new one

pool = gcp('nocreate');

if isempty(pool)
    pool = parpool(n_workers);
elseif pool.NumWorkers ~= n_workers
    fprintf('Wrong number of workers. Starting again.\n');
    delete(pool);
    pool = parpool(n_workers);
else
    fprintf('Pool with %d workers already runing.\n', n_workers);
end

end