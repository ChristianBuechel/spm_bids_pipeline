function y = get_pcs(y,ex_var)
y       = spm_detrend(y,1);  % remove mean for PCA
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