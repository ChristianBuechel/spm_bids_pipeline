function remove = find_aborted_runs(ind,sess,thresh)

%find runs that have fewer than expected scans
abort = find(sess<thresh);
remove = [];

if ~isempty(abort)
    remove = [];
    disp('Deleting aborted runs')
    for ar = 1:length(abort)
        remove = [remove ind(abort(ar))+1:ind(abort(ar))+sess(abort(ar))];
    end
end
