function n_noise = detect_noise_volumes(X)

% function n_noise = detect_noise_volumes(X)
% number of noise volumes (acquired without excitation) at the END of a 4D magnitude series.
% X is a file name (4D NIfTI) or a 4D array (magnitude, or complex: abs is used).
%
% Criterion: for each volume the mean magnitude over the brain/object voxels (top 20% of the mean
% of the first 10 volumes). Signal volumes are very constant (robust sd < 1.5% in the test data),
% noise volumes are much lower. The last volumes are counted as noise while
%   mean < 0.9 * median   and   mean < median - 10 * robust sd   (robust sd at least 0.5% of the median)
% At most 10 volumes (or a quarter of the series) are considered.
% Test data (VASO, 0.82 mm, tSNR ~2): vasoref cbv/bold: z = -66 / -226 for the 2 noise volumes,
% signal volumes z > -3; pRFVASO (no noise volumes): z > -2, nothing detected.

if ischar(X) || isstring(X)
    V = spm_vol(char(X));  nt = numel(V);
    get = @(n) abs(double(spm_read_vols(V(n))));
else
    nt = size(X,4);
    get = @(n) abs(double(X(:,:,:,n)));
end
n_noise = 0;
if nt < 12, return, end
nmax = min(10, floor(nt/4));

ref = 0;
for n = 1:10, ref = ref + get(n); end
ref  = ref(:)/10;
s    = sort(ref);
mask = ref >= s(round(0.8*numel(s)));
clear s ref

m = zeros(nt,1);
for n = 1:nt
    y = get(n);
    m(n) = mean(y(mask));
end
med = median(m);
sd  = max(1.4826*median(abs(m - med)), 0.005*med);
isn = m < 0.9*med & (m - med) < -10*sd;
while n_noise < nmax && isn(nt - n_noise)
    n_noise = n_noise + 1;
end
end
