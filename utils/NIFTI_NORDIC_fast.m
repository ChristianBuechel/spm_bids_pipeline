function NIFTI_NORDIC_fast(fn_magn_in, fn_phase_in, prefix, ARG)
% NIFTI_NORDIC_fast  NORDIC / MP-PCA denoising of 4D magnitude or complex NIfTI data.
%
%   NIFTI_NORDIC_fast(fn_magn_in, fn_phase_in, prefix, ARG)
%
%   fn_magn_in / fn_phase_in : NIfTI file names (.nii or .nii.gz, with or
%                              without path; an SPM frame suffix ',1' is ignored)
%   prefix                   : prepended to the input file name(s), e.g. 'n_'
%                              /data/sub01/run1.nii -> /data/sub01/n_run1.nii
%   Outputs are written next to the magnitude input unless ARG.DIROUT is set.
%
%   fMRI:            ARG.temporal_phase = 1; ARG.phase_filter_width = 10;
%   dMRI:            ARG.temporal_phase = 3; ARG.phase_filter_width = 3;
%   magnitude only:  ARG.magnitude_only = 1;  (fn_phase_in ignored, may be [])
%
% OPTIONS (all optional, defaults in [])
%   ARG.DIROUT                 output directory           [folder of fn_magn_in]
%   ARG.noise_volume_last      noise volume, counted from the end            [0]
%   ARG.noise_volume_split     n or 'auto': the last n volumes are noise volumes (no excitation); 'auto'
%                              detects them (detect_noise_volumes), a number is checked against it. They are
%                              NOT part of the data: they are removed, used for the noise level
%                              and saved as noise_<magn>.nii (+ noise_<phase>.nii) and noise_<magn>.json
%                              (noise level in physical units, for other series). Output has nt-n volumes  [0]
%   ARG.noise_level_file       noise_<magn>.json of another series (same sequence/coil/GRAPPA, same
%                              session): its noise level is used if this series has no noise volumes  ['']
%   ARG.factor_error           >1 higher, <1 lower noise floor               [1]
%   ARG.full_dynamic_range     0 keep input scale, 1 use integer range       [0]
%   ARG.temporal_phase         0 off, 1 standard, 2/3 extra spike handling   [1]
%   ARG.phase_filter_width     exponent of the Tukey (Hann) phase filter     [3]
%   ARG.NORDIC / ARG.MP        NORDIC=1 (default) | MP=1 MP-PCA | MP=2 no g-factor
%   ARG.kernel_size_gfactor    [k1 k2 k3 nvol]                               [14 14 1 90]
%   ARG.kernel_size_PCA        [k1 k2 k3]                    [11:1 voxels:volumes]
%   ARG.NORDIC_patch_overlap / ARG.gfactor_patch_overlap                     [2]
%   ARG.magnitude_only         1 = magnitude data only                       [0]
%   ARG.phase_slice_average_for_kspace_centering                             [0]
%   ARG.use_magn_for_gfactor   if field exists: g-factor from |data|
%   ARG.data_has_zero_elements 1 = fill all-zero voxels with noise           [0]
%   ARG.save_gfactor_map       1 write gfactor_<prefix><magn>.nii, 2 + stop  []
%   ARG.save_add_info          1 = save ARG structs to <prefix><magn>.mat    [0]
%   ARG.save_residual_matlab   1 = save residual to RESIDUAL_<prefix><magn>.mat [0]
%   ARG.make_complex_nii       1 = also write <prefix><phase>.nii            [0]
%   ARG.write_gzipped_niftis   1 = write .nii.gz                             [0]
%   ARG.output_datatype        'input' | 'int16' | 'uint16' | 'single'       ['input']
%   ARG.svd_method             'svd' | 'gram' (faster, see llr_slab)         ['svd']
%
% Requires R2017b+ (niftiread/niftiwrite). Signal Processing Toolbox no longer needed.
% Based on NIFTI_NORDIC.m, Copyright Board of Regents, University of Minnesota.

%% ---------------------------------------------------------------- options
if nargin < 4 || isempty(ARG), ARG = struct(); end
% original semantics: existence of the field switched complex output on
make_complex = isfield(ARG,'make_complex_nii') && ...
    (isempty(ARG.make_complex_nii) || ~isequal(ARG.make_complex_nii, 0));
ARG = set_defaults(ARG);
is_cplx = ARG.magnitude_only ~= 1;
if ~is_cplx, ARG.temporal_phase = 0; make_complex = false; end
nvl = ARG.noise_volume_last;

%% ---------------------------------------------------------------- load
fn_magn_in = clean_name(fn_magn_in);
[outdir, base_magn] = out_name(fn_magn_in, prefix, ARG);
if is_cplx
    fn_phase_in = clean_name(fn_phase_in);
    [~, base_phase] = out_name(fn_phase_in, prefix, ARG);
end
info = niftiinfo(fn_magn_in);
M = niftiread(fn_magn_in);                 % stored values; scaling stays in header
matdim = size(M); matdim(end+1:4) = 1; nt = matdim(4);

if is_cplx
    info_phase = niftiinfo(fn_phase_in);
    P = niftiread(fn_phase_in);
    if ~isequal(size(P), size(M)), error('NORDIC:size','Magnitude/phase size differ.'); end
    pmax = single(max(P(:))); pmin = single(min(P(:)));
    range_norm   = pmax - pmin;
    range_center = (pmax + pmin) / range_norm / 2;
    % build complex volume-wise from native (int16) arrays: no full single copies
    K = zeros(matdim, 'like', complex(single(0)));
    for n = 1:nt
        ph = (single(P(:,:,:,n)) ./ range_norm - range_center) * 2*pi;
        K(:,:,:,n) = abs(single(M(:,:,:,n))) .* exp(1i*ph);
    end
    clear P ph
else
    K = abs(single(M));
end
clear M

v1 = abs(K(:,:,:,1)); ARG.ABSOLUTE_SCALE = min(v1(v1 ~= 0)); clear v1
if isempty(ARG.ABSOLUTE_SCALE), error('NORDIC:empty','First volume is all zero.'); end
for n = 1:nt, K(:,:,:,n) = K(:,:,:,n) ./ ARG.ABSOLUTE_SCALE; end

Kn = [];                                   % noise volumes (normalised like the data), kept apart
ns = ARG.noise_volume_split;
if ischar(ns) || isstring(ns)              % 'auto': NORDIC decides which of the last volumes are noise
    if ~strcmpi(ns,'auto'), error('NORDIC:noise','noise_volume_split: a number or ''auto'''); end
    ns = detect_noise_volumes(K);
    fprintf('noise volumes at the end of the series (automatic): %d\n',ns);
    ARG.noise_volume_split = ns;
elseif ns > 0                              % number given: check it
    nd = detect_noise_volumes(K);
    if nd ~= ns
        warning('NORDIC:noise','%d noise volumes requested, but %d detected at the end of %s.',ns,nd,fn_magn_in);
    end
end
if ns > 0
    if nvl > 0, error('NORDIC:noise','noise_volume_split and noise_volume_last exclude each other.'); end
    if nt - ns < 6, error('NORDIC:noise','Only %d volumes left after splitting off %d noise volumes.',nt-ns,ns); end
    Kn = K(:,:,:,nt-ns+1:nt);
    K  = K(:,:,:,1:nt-ns);
    nt = nt - ns;  matdim(4) = nt;
end
slope = info.MultiplicativeScaling;  if slope == 0, slope = 1; end
if nt < 6, disp('Too few volumes'); end

%% ---------------------------------------------------------------- phase maps
use_meanphase = is_cplx && ARG.phase_slice_average_for_kspace_centering ~= 0;
if use_meanphase
    acc = zeros(matdim(1:3), 'like', K);
    for n = 1:nt-nvl, acc = acc + K(:,:,:,n); end
    meanphase = angle(acc/(nt-nvl) * ARG.phase_slice_average_for_kspace_centering);
    clear acc
    % as in the original (KKSP2 typo): mean phase is only removed in the NORDIC stage
end

DD_phase = [];                             % stored as angle: single, 4 B/voxel
if ARG.temporal_phase > 0
    disp('estimating slice-dependent phases ...')
    wx = tukey1(matdim(1)) .^ ARG.phase_filter_width;
    wy = reshape(tukey1(matdim(2)) .^ ARG.phase_filter_width, 1, []);
    DD_phase = zeros(matdim, 'single');
    for n = 1:nt
        raw = K(:,:,:,n);
        tmp = lowpass_2d(raw, wx, wy);
        if ARG.temporal_phase == 2
            mask = abs(angle(raw ./ tmp)) > 1;
            tmp(mask) = raw(mask);
        end
        DD_phase(:,:,:,n) = angle(tmp);
    end
    clear raw tmp mask
end

%% ---------------------------------------------------------------- g-factor
kdef = [14 14 1];
if isempty(ARG.kernel_size_gfactor)
    kg = kdef; nG = 90;
else
    kg = ARG.kernel_size_gfactor(1:min(3,end));
    kg(end+1:3) = kdef(numel(kg)+1:3);
    if numel(ARG.kernel_size_gfactor) >= 4, nG = ARG.kernel_size_gfactor(4); else, nG = 90; end
end
nG = min(nG, nt);

Kg = K(:,:,:,1:nG);                        % only the volumes the g-factor needs
if ~isempty(DD_phase)
    for n = 1:nG, Kg(:,:,:,n) = Kg(:,:,:,n) .* exp(-1i*DD_phase(:,:,:,n)); end
end
if isfield(ARG,'use_magn_for_gfactor'), Kg = abs(Kg); end
Kg(~isfinite(Kg)) = 0;

ARG2 = ARG;
ARG2.kernel_size = clamp_kernel(kg, matdim);
ARG2.patch_average_sub = ARG.gfactor_patch_overlap;
ARG2.LLR_scale = 0;  ARG2.NVR_threshold = 1;  ARG2.soft_thrs = 10;   % MP-PCA

disp('estimating g-factor ...')
[~, st] = llr_denoise(Kg, ARG2, false);  clear Kg
ARG2.NOISE               = sqrt(st.NOISE ./ st.weight);
ARG2.Component_threshold = st.Component_threshold ./ st.weight;
ARG2.energy_removed      = st.energy_removed ./ st.weight;
ARG2.SNR_weight          = st.SNR_weight ./ st.weight;
gfactor = ARG2.NOISE;

if nG < 6
    gfactor(isnan(gfactor)) = 0;
    gfactor(gfactor == 0) = median(gfactor(gfactor ~= 0));
end
if any(gfactor(:) == 0)
    gfactor(isnan(gfactor)) = 0;
    gfactor(gfactor < 1) = median(gfactor(gfactor ~= 0));
    ARG.data_has_zero_elements = 1;
end
if ~any(isfinite(gfactor(:)) & gfactor(:) > 0)
    error('NORDIC:gfactor','g-factor estimation failed (all zero/NaN).');   % fail loudly
end
if ARG.MP == 2, gfactor = ones(size(gfactor), 'single'); end

if ~isempty(ARG.save_gfactor_map) && any(ARG.save_gfactor_map == [1 2])
    write_gfactor(gfactor, info, ARG, fullfile(outdir, ['gfactor_' base_magn '.nii']));
    if ARG.save_gfactor_map == 2, return; end
end

%% ---------------------------------------------------------------- prepare (in place)
noise_vol = [];
for n = 1:nt
    x = K(:,:,:,n);
    if use_meanphase, x = x .* exp(-1i*meanphase); end
    x = x ./ gfactor;
    if nvl > 0 && n == nt+1-nvl, noise_vol = x; end
    if ~isempty(DD_phase)
        dd = DD_phase(:,:,:,n);
        if ARG.temporal_phase == 3
            mask = abs(angle(x .* exp(-1i*dd))) > 1 & abs(x) > sqrt(2);
            dd(mask) = angle(x(mask));
            DD_phase(:,:,:,n) = dd;
        end
        x = x .* exp(-1i*dd);
    end
    x(~isfinite(x)) = 0;
    K(:,:,:,n) = x;
end
clear x dd mask

own_noise = false;
if ns > 0                                  % noise volumes of this series (same g-factor normalisation as the data)
    z = [];
    for q = 1:ns
        x = Kn(:,:,:,q) ./ gfactor;  x(~isfinite(x)) = 0;
        z = [z; x(x ~= 0)]; %#ok<AGROW>
    end
    ARG.measured_noise = std(z);  own_noise = true;
    clear z x
elseif nvl > 0
    noise_vol(~isfinite(noise_vol)) = 0;
    ARG.measured_noise = std(noise_vol(noise_vol ~= 0));  own_noise = true;
elseif ~isempty(ARG.noise_level_file)      % noise level of another series, physical units -> units of this series
    if isfile(ARG.noise_level_file)
        L = jsondecode(fileread(ARG.noise_level_file));
        ARG.measured_noise = L.NoiseLevelPhysical / (ARG.ABSOLUTE_SCALE * slope);
        fprintf('noise level from %s\n',ARG.noise_level_file);
    else
        warning('NORDIC:noisefile','Noise level file %s not found: assuming unit noise.',ARG.noise_level_file);
        ARG.measured_noise = 1;
    end
else
    ARG.measured_noise = 1;
end
if own_noise, ARG.noise_level_phys = ARG.measured_noise * ARG.ABSOLUTE_SCALE * slope; end   % before the complex factor
clear noise_vol
if ~isfield(ARG,'use_magn_for_gfactor') && is_cplx
    ARG.measured_noise = ARG.measured_noise / sqrt(2);
end

if ARG.data_has_zero_elements == 1
    s = zeros(matdim(1:3), 'single');
    for n = 1:nt, s = s + abs(K(:,:,:,n)); end
    MASK = (s == 0); nz0 = nnz(MASK); clear s
    if nz0 > 0
        for n = 1:nt
            x = K(:,:,:,n);
            if is_cplx, x(MASK) = (randn(nz0,1) + 1i*randn(nz0,1)) / sqrt(2);
            else,       x(MASK) = randn(nz0,1);        % keeps magnitude data real
            end
            K(:,:,:,n) = x;
        end
    end
    clear MASK x
end

%% ---------------------------------------------------------------- NORDIC
if isempty(ARG.kernel_size_PCA), ks = repmat(round((nt*11)^(1/3)), 1, 3);
else, ks = ARG.kernel_size_PCA; end
if matdim(3) <= ks(3)
    ks = repmat(round((nt*11/matdim(3))^(1/2)), 1, 2);  ks(3) = matdim(3);
end
ARG.kernel_size = clamp_kernel(ks, matdim);
ARG.patch_average_sub = ARG.NORDIC_patch_overlap;
ARG.LLR_scale = 1;

thr = 0;
for k = 1:10
    s = svd(randn(prod(ARG.kernel_size), nt));      % values only
    thr = thr + s(1);
end
if is_cplx, ARG.NVR_threshold = thr/10*sqrt(2)*ARG.measured_noise*ARG.factor_error;
else,       ARG.NVR_threshold = thr/10*ARG.measured_noise*ARG.factor_error; end
ARG.soft_thrs = [];
if ARG.MP > 0, ARG.soft_thrs = 10; end
if isfield(ARG,'soft_thrs_in'), ARG.soft_thrs = ARG.soft_thrs_in; end

if ARG.save_residual_matlab == 1, Residual = K; end  % extra copy only on request

disp('starting NORDIC ...')
[K, st] = llr_denoise(K, ARG, true);       % K now holds the denoised data
ARG.NOISE               = sqrt(st.NOISE ./ st.weight);
ARG.Component_threshold = st.Component_threshold ./ st.weight;
ARG.energy_removed      = st.energy_removed ./ st.weight;
ARG.SNR_weight          = st.SNR_weight ./ st.weight;
clear st
disp('completing NORDIC ...')

if ARG.save_residual_matlab == 1
    Residual = Residual - K; %#ok<NASGU>
    save(fullfile(outdir, ['RESIDUAL_' base_magn '.mat']), 'Residual', '-v7.3');
    clear Residual
end

%% ---------------------------------------------------------------- undo normalisation
for n = 1:nt
    x = K(:,:,:,n) .* gfactor;
    if use_meanphase,      x = x .* exp(1i*meanphase); end
    if ~isempty(DD_phase), x = x .* exp(1i*DD_phase(:,:,:,n)); end
    x = x .* ARG.ABSOLUTE_SCALE;
    x(~isfinite(x)) = 0;
    K(:,:,:,n) = x;
end
clear x DD_phase gfactor

%% ---------------------------------------------------------------- write
write_magnitude(K, info, ARG, fullfile(outdir, [base_magn '.nii']));
if make_complex
    write_phase(K, info_phase, ARG, range_norm, range_center, fullfile(outdir, [base_phase '.nii']));
end
if ns > 0
    ph = [];
    if is_cplx, ph = struct('info',info_phase,'range_norm',range_norm,'range_center',range_center,'fn',fn_phase_in); end
    write_noise(Kn, info, ph, ARG, outdir, fn_magn_in);
end
if ARG.save_add_info == 1
    save(fullfile(outdir, [base_magn '.mat']), 'ARG2', 'ARG', '-v7.3');
end
end

%% ========================================================================
function ARG = set_defaults(ARG)
def = struct('DIROUT','', 'noise_volume_last',0, 'noise_volume_split',0, 'noise_level_file','', 'factor_error',1, ...
    'full_dynamic_range',0, 'temporal_phase',1, 'phase_filter_width',3, ...
    'NORDIC_patch_overlap',2, 'gfactor_patch_overlap',2, ...
    'kernel_size_gfactor',[], 'kernel_size_PCA',[], ...
    'phase_slice_average_for_kspace_centering',0, 'magnitude_only',0, ...
    'save_gfactor_map',[], 'save_add_info',0, 'save_residual_matlab',0, ...
    'data_has_zero_elements',0, 'write_gzipped_niftis',0, ...
    'output_datatype','input', 'svd_method','svd');
f = fieldnames(def);
for k = 1:numel(f)
    if ~isfield(ARG, f{k}), ARG.(f{k}) = def.(f{k}); end
end
if isempty(ARG.magnitude_only), ARG.magnitude_only = 0; end
if ~isfield(ARG,'NORDIC') && ~isfield(ARG,'MP'), ARG.NORDIC = 1; ARG.MP = 0;
elseif ~isfield(ARG,'NORDIC'), ARG.NORDIC = double(ARG.MP ~= 1);
elseif ~isfield(ARG,'MP'),     ARG.MP = double(ARG.NORDIC ~= 1);
end
end

function fn = clean_name(fn)
% char, no surrounding blanks, no SPM frame suffix (',1')
fn = regexprep(strtrim(char(fn)), ',\d+$', '');
end

function [outdir, base] = out_name(fn, prefix, ARG)
% output folder and output file name WITHOUT extension: <prefix><name of fn>
if isempty(prefix)
    error('NORDIC:prefix', 'An empty prefix would overwrite the input file.');
end
[p, n, e] = fileparts(fn);
if strcmpi(e, '.gz'), [~, n] = fileparts(n); end        % name.nii.gz -> name
base = [char(prefix) n];
if ~isempty(ARG.DIROUT), outdir = ARG.DIROUT; else, outdir = p; end
if isempty(outdir), outdir = pwd; end
if ~exist(outdir, 'dir'), mkdir(outdir); end
end

function ks = clamp_kernel(ks, matdim)
ks = min(reshape(ks(1:3),1,3), matdim(1:3));
end

function w = tukey1(N)          % == tukeywin(N,1) == hann(N)
if N == 1, w = 1; else, w = 0.5*(1 - cos(2*pi*(0:N-1)'/(N-1))); end
end

function x = lowpass_2d(x, wx, wy)   % same ops as original, vectorised over slices
for d = 1:2, x = ifftshift(ifft(ifftshift(x,d),[],d),d); end
x = x .* wy;  x = x .* wx;
for d = 1:2, x = fftshift(fft(fftshift(x,d),[],d),d); end
end

function [K, st] = llr_denoise(K, ARG, write_back)
% Slab-wise LLR along dim 1. With write_back, finished rows (never read again)
% are written into K, so only a slab-sized accumulator is needed.
[nx, ny, nz, nt] = size(K);
k1 = ARG.kernel_size(1);
z = zeros(nx, ny, nz, 'single');
st = struct('weight',z, 'NOISE',z, 'Component_threshold',z, 'energy_removed',z, 'SNR_weight',z);
clear z
step   = max(1, floor(k1/ARG.patch_average_sub));
starts = unique([1:step:nx-k1+1, nx-k1+1]);    % original stopped at nx-k1 (last row lost)
if write_back, acc = zeros([k1 ny nz nt], 'like', K); base = 1; end
for s = starts
    if write_back && s > base
        d = s - base;  rows = base:s-1;
        K(rows,:,:,:) = acc(1:d,:,:,:) ./ st.weight(rows,:,:);
        acc(1:k1-d,:,:,:)    = acc(d+1:k1,:,:,:);
        acc(k1-d+1:k1,:,:,:) = 0;
        base = s;
    end
    r = s:s+k1-1;
    [upd, st] = llr_slab(K(r,:,:,:), ARG, st, r, write_back);
    if write_back, acc = acc + upd; end
end
if write_back
    rows = base:nx;
    K(rows,:,:,:) = acc(1:numel(rows),:,:,:) ./ st.weight(rows,:,:);
end
end

function [upd, st] = llr_slab(Ks, ARG, st, r, need_recon)
% Thresholding identical to subfunction_loop_for_NVR_avg_update.
% svd_method 'gram': eig of double-precision X'*X (nt x nt) instead of svd(X);
% usually several times faster for k^3 >> nt. Validate against 'svd' first.
[~, ny, nz, nt] = size(Ks);
w2 = ARG.kernel_size(2);  w3 = ARG.kernel_size(3);
lambda = ARG.LLR_scale * ARG.NVR_threshold;
soft = ARG.soft_thrs;
fixed_scale = ~isfield(ARG,'patch_scale') || ARG.patch_scale == 1;
use_gram = strcmpi(ARG.svd_method, 'gram');
if need_recon, upd = zeros(size(Ks), 'like', Ks); else, upd = []; end
W = st.weight(r,:,:); NS = st.NOISE(r,:,:); CT = st.Component_threshold(r,:,:);
ER = st.energy_removed(r,:,:); SW = st.SNR_weight(r,:,:);
n2list = [1:max(1,floor(w2/ARG.patch_average_sub)):ny-w2+1, ny-w2+1];
n3list = [1:max(1,floor(w3/ARG.patch_average_sub)):nz-w3+1, nz-w3+1];

for n2 = n2list
    i2 = n2:n2+w2-1;
    for n3 = n3list
        i3 = n3:n3+w3-1;
        Xp = Ks(:,i2,i3,:);
        X = reshape(Xp, [], nt);  m = size(X,1);

        if ~need_recon
            S = svd(X); U = []; V = [];                     % values only
        elseif use_gram
            G = double(X)'*double(X);  G = (G+G')/2;
            [V, D] = eig(G);
            [dg, o] = sort(real(diag(D)), 'descend');
            p = min(m, nt);
            V = V(:,o(1:p));  S = single(sqrt(max(dg(1:p),0)));  U = [];
        else
            [U, S, V] = svd(X, 'econ');  S = diag(S);
        end
        S0 = S;

        idx = sum(S < lambda);  sig2 = [];
        if isempty(soft)                                    % NORDIC
            energy = sqrt(sum(S)) \ sqrt(sum(S(S < lambda)));
            S(S < lambda) = 0;  t = idx;
        elseif soft ~= 10                                   % soft threshold
            S = S - lambda*soft;  S(S < 0) = 0;  energy = 0;  t = 1;
        else                                                % MP-PCA
            MM = m - sum(sum(X,2) == 0);
            if MM > 0
                R = min(MM, nt);
                scaling   = (max(MM,nt) - (0:R-1)).' / nt;
                vals      = S.^2 / nt;
                csum      = cumsum(vals(R:-1:1));
                sigmasq_1 = (csum(R:-1:1) ./ (R:-1:1).') ./ scaling;
                sig2      = (vals(1:R) - vals(R)) ./ (4*sqrt((MM-(0:R-1)).'/nt));
                t = find(sig2 < sigmasq_1, 1);
                if isempty(t)
                    idx = 0;  energy = 0;
                else
                    idx = numel(S) - t + 1;
                    energy = sqrt(sum(S)) \ sqrt(sum(S(t:end)));
                    S(t:end) = 0;
                end
            else
                t = 1;  energy = 0;  sig2 = 0;
            end
        end
        if isempty(t), t = 1; end
        if fixed_scale, ps = 1; else, ps = numel(S) - idx; end

        if need_recon                         % reconstruct from kept components only
            if use_gram
                f = zeros(size(S), 'single');  nzs = S0 > 0;  f(nzs) = S(nzs)./S0(nzs);
                k = find(f ~= 0);  Vk = cast(V(:,k), 'like', X);
                Y = X * ((Vk .* f(k).') * Vk');
            else
                k = find(S ~= 0);
                Y = U(:,k) * (S(k) .* V(:,k)');
            end
            upd(:,i2,i3,:) = upd(:,i2,i3,:) + ps*reshape(Y, size(Xp));
        end

        W(:,i2,i3)  = W(:,i2,i3)  + ps;
        CT(:,i2,i3) = CT(:,i2,i3) + idx;
        ER(:,i2,i3) = ER(:,i2,i3) + energy;
        SW(:,i2,i3) = SW(:,i2,i3) + S(1)./S(max(1,t-1));
        if ~isempty(sig2), NS(:,i2,i3) = NS(:,i2,i3) + sig2(t); end
    end
end
st.weight(r,:,:) = W;  st.NOISE(r,:,:) = NS;  st.Component_threshold(r,:,:) = CT;
st.energy_removed(r,:,:) = ER;  st.SNR_weight(r,:,:) = SW;
end

function cls = out_class(opt, in_dt)
ints = {'int8','uint8','int16','uint16','int32','uint32'};
switch lower(opt)
    case 'input'
        if any(strcmp(in_dt, ints)), cls = in_dt; else, cls = 'single'; end
    case {'single','float32','float'}
        cls = 'single';
    otherwise
        if ~any(strcmp(opt, ints)), error('NORDIC:dtype','Unknown output_datatype %s', opt); end
        cls = opt;
end
end

function b = nbits(cls)
switch cls
    case {'int8','uint8'},   b = 8;
    case {'int16','uint16'}, b = 16;
    otherwise,               b = 32;
end
end

function write_magnitude(K, info, ARG, fname)
nt = size(K,4);
cls = out_class(ARG.output_datatype, info.Datatype);
is_int = ~strcmp(cls, 'single');
mx = 0;
for n = 1:nt, mx = max(mx, max(reshape(abs(K(:,:,:,n)),[],1))); end
g = 1;
if is_int && ARG.full_dynamic_range ~= 0
    stp = max(1, floor(numel(K)/2e7));                 % subsample instead of full sort
    a = sort(abs(K(1:stp:end)));
    target = 32000 * (1 + strcmp(cls,'uint16'));
    g = 2^floor(log2(target / (2*a(round(0.99*end)))));
    clear a
end
if is_int && mx*g > double(intmax(cls))
    warning('NORDIC:range', 'Max %.0f exceeds %s range: rescaling + scl_slope instead of clipping.', mx*g, cls);
    g = double(intmax(cls)) / mx;
end
out = zeros(size(K), cls);
for n = 1:nt, out(:,:,:,n) = cast(abs(K(:,:,:,n))*g, cls); end   % rounds
info.Datatype = cls;  info.BitsPerPixel = nbits(cls);
if g ~= 1                                   % keep physical intensities = input
    s0 = info.MultiplicativeScaling;  if s0 == 0, s0 = 1; end
    info.MultiplicativeScaling = s0 / g;
end
info = fit_info(info, size(out));
niftiwrite(out, fname, info, 'Compressed', logical(ARG.write_gzipped_niftis));
end

function write_phase(K, info_phase, ARG, range_norm, range_center, fname)
nt = size(K,4);
cls = out_class('input', info_phase.Datatype);
out = zeros(size(K), cls);
for n = 1:nt
    out(:,:,:,n) = cast((angle(K(:,:,:,n))/(2*pi) + range_center) * range_norm, cls);
end
info_phase.Datatype = cls;  info_phase.BitsPerPixel = nbits(cls);
info_phase = fit_info(info_phase, size(out));
niftiwrite(out, fname, info_phase, 'Compressed', logical(ARG.write_gzipped_niftis));
end

function info = fit_info(info, sz)
% header for an array of size sz (a different number of volumes than the input)
sz = sz(1:find(sz > 1, 1, 'last'));  sz(end+1:3) = 1;
info.ImageSize = sz;
info.PixelDimensions = info.PixelDimensions(1:numel(sz));
end

function write_noise(Kn, info, ph, ARG, outdir, fn_magn)
% the noise volumes as separate images (stored units of the input), plus the noise level (json)
[~, bm] = out_name(fn_magn, 'noise_', ARG);
Kn = Kn * ARG.ABSOLUTE_SCALE;
cls = out_class('input', info.Datatype);
out = zeros(size(Kn), cls);
for n = 1:size(Kn,4), out(:,:,:,n) = cast(abs(Kn(:,:,:,n)), cls); end
i1 = fit_info(info, size(Kn));  i1.Datatype = cls;  i1.BitsPerPixel = nbits(cls);
niftiwrite(out, fullfile(outdir, [bm '.nii']), i1, 'Compressed', logical(ARG.write_gzipped_niftis));
if ~isempty(ph)
    [~, bp] = out_name(ph.fn, 'noise_', ARG);
    clsp = out_class('input', ph.info.Datatype);
    outp = zeros(size(Kn), clsp);
    for n = 1:size(Kn,4)
        outp(:,:,:,n) = cast((angle(Kn(:,:,:,n))/(2*pi) + ph.range_center) * ph.range_norm, clsp);
    end
    ip = fit_info(ph.info, size(Kn));  ip.Datatype = clsp;  ip.BitsPerPixel = nbits(clsp);
    niftiwrite(outp, fullfile(outdir, [bp '.nii']), ip, 'Compressed', logical(ARG.write_gzipped_niftis));
end
L = struct('NoiseLevelPhysical',ARG.noise_level_phys, 'NoiseVolumes',size(Kn,4), 'Source',fn_magn, ...
    'Note','std of the g-factor normalised noise volumes (before the complex factor 1/sqrt(2)), physical units (stored value * scl_slope)');
fid = fopen(fullfile(outdir, [bm '.json']), 'w');  fwrite(fid, jsonencode(L));  fclose(fid);
end

function write_gfactor(g, info, ARG, fname)
% float with input geometry (original: no header, and integer -> g~1 became 0/1/2)
g = single(abs(g));  g(~isfinite(g)) = 0;
i3 = info;
i3.ImageSize = size(g);  i3.PixelDimensions = info.PixelDimensions(1:ndims(g));
i3.Datatype = 'single';  i3.BitsPerPixel = 32;
i3.MultiplicativeScaling = 1;  i3.AdditiveOffset = 0;
niftiwrite(g, fname, i3, ...
    'Compressed', logical(ARG.write_gzipped_niftis));
end
