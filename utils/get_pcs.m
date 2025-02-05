function y = get_pcs(y,ex_var)
% get principle components from time-series
% if ex_var > 1 this is the number of PCs we return
% if ex_var < 1, we take the first n PCs so that ex_var of the variance is explained
% one caveat is that in unsmoothed data, the correlation of neighboring
% voxels is not too high, which means if ex_var is eg 0.9 we get a LOT of
% PCs
% if we apply spm_filter like in 1st level analysis we can reduce the # of
% PCs but only a bit

 y      = spm_detrend(y,1);  % remove mean for PCA
[u,s,u] = svd(y*y');
if ex_var > 1
    n_pc = ex_var;
else
    n_pc = min(find(cumsum(diag(s))./sum(diag(s))>ex_var));
end
s       = diag(s);
u       = u(:,1:n_pc);
v       = y'*u/sqrt(s(1:n_pc)');
d       = sign(sum(v));
y       = u*d.*sqrt(size(y,1)-1); % normalize