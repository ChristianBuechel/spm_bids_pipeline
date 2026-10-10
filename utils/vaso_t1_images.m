function vaso_t1_images(job)

% function vaso_t1_images(job)
% worker of sbp_vaso_T1: builds T1-like images from the (denoised, realigned, coregistered)
% cbv (blood nulled) and bold (not nulled) series of the vaso reference runs.
% All volumes are resampled (linear) to the grid of the reference image (mean bold) with
% the header of each volume (sidecar .mat of the realignment/coregistration).
%
%   job.mag_c, job.ph_c   cbv magnitude and phase 4D images (phase in Siemens units, +-4096 = +-pi)
%   job.mag_b, job.ph_b   bold magnitude and phase 4D images
%   job.ref               reference image (mean bold, defines the grid)
%   job.out               struct with the output file names: median, meanstd, signed, signedpos
%   job.thresh            voxels with a median bold signal below thresh * (99th percentile) are set to 0 (0.05)
%   job.phase_fwhm        smoothing [mm] of the complex reference (temporal mean of the bold series), default 0
%   job.phase_scale       phase value to radians, default pi/4096
%   job.pos_scale         scaling of the positive version of the signed image, default 1000
%   job.n_noise           number of noise volumes at the end of each series (not used), default 0
%
% Images (float, same header as the reference, descrip starts with "derived"):
%   median    median(cbv) / median(bold)  (values >= 1 set to 0)
%   meanstd   mean / std over all volumes of both series
%   signed    the sign of the nulled images is restored with the phase: per cbv volume
%             Re( M_cbv * exp(i*(phase_cbv - phase_ref)) ), phase_ref = phase of the temporal
%             mean of the complex bold series; median over time / median(bold). CSF is negative.
%   signedpos the signed image with the range between the pos_pct percentiles ([1 99]) of the
%             signed values mapped to 0..pos_scale: positive, for tools that cannot handle negative values

if ~isfield(job,'n_noise') || isempty(job.n_noise), job.n_noise = 0; end
if ~isfield(job,'thresh') || isempty(job.thresh), job.thresh = 0.05; end
if ~isfield(job,'phase_fwhm') || isempty(job.phase_fwhm), job.phase_fwhm = 0; end
if ~isfield(job,'phase_scale') || isempty(job.phase_scale), job.phase_scale = pi/4096; end
if ~isfield(job,'pos_scale') || isempty(job.pos_scale), job.pos_scale = 1000; end
if ~isfield(job,'pos_pct') || isempty(job.pos_pct), job.pos_pct = [1 99]; end

Vref = spm_vol(job.ref);
Vref = Vref(1);                                    % (a 3D image, normally the mean bold)
dim  = Vref.dim(1:3);
N    = prod(dim);
[X,Y,Z] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
W    = Vref.mat*[X(:)';Y(:)';Z(:)';ones(1,N)];     % world coordinates of the reference voxels
clear X Y Z

Vmc = spm_vol(job.mag_c); Vpc = spm_vol(job.ph_c);
Vmb = spm_vol(job.mag_b); Vpb = spm_vol(job.ph_b);
if numel(Vpc) ~= numel(Vmc) || numel(Vpb) ~= numel(Vmb)
    error('vaso_t1_images: magnitude and phase images have different numbers of volumes');
end
nc  = numel(Vmc) - job.n_noise;                    % the noise volumes at the end are not used
nb  = numel(Vmb) - job.n_noise;

Mb      = zeros(N,nb,'single');                    % magnitude stacks (for the medians)
Mc      = zeros(N,nc,'single');
Sc      = zeros(N,nc,'single');                    % signed cbv
sum_all = zeros(N,1);  sq_all = zeros(N,1);        % running sums for mean/std over both series
n_all   = nc + nb;
sumRe   = zeros(N,1,'single'); sumIm = zeros(N,1,'single');

%% pass 1: bold series (magnitude stack, complex temporal mean)
for t = 1:nb
    [m,re,im] = resample_frame(Vmb(t),Vpb(t),W,dim,job.phase_scale);
    Mb(:,t)   = m;
    sumRe     = sumRe + re;  sumIm = sumIm + im;
    sum_all   = sum_all + double(m);  sq_all = sq_all + double(m).^2;
end
ref = complex(sumRe,sumIm)/nb;
if job.phase_fwhm > 0                              % optional smoothing of the complex reference
    vox   = sqrt(sum(Vref.mat(1:3,1:3).^2));
    ref3  = reshape(ref,dim);
    re3   = zeros(dim,'single'); im3 = zeros(dim,'single');
    spm_smooth(real(ref3),re3,job.phase_fwhm./vox);
    spm_smooth(imag(ref3),im3,job.phase_fwhm./vox);
    ref   = complex(re3(:),im3(:));
end
a = abs(ref);  u = zeros(N,1,'single');  ok = a > 0;
u(ok) = ref(ok)./a(ok);                            % unit phasor of the reference phase

%% pass 2: cbv series (magnitude stack, signed values)
for t = 1:nc
    [m,re,im] = resample_frame(Vmc(t),Vpc(t),W,dim,job.phase_scale);
    Mc(:,t)   = m;
    Sc(:,t)   = re.*real(u) + im.*imag(u);         % Re( Z * conj(u) ) = M cos(phase - phase_ref)
    sum_all   = sum_all + double(m);  sq_all = sq_all + double(m).^2;
end

%% images
med_c = median(Mc,2);  med_b = median(Mb,2);  med_s = median(Sc,2);
clear Mc Mb Sc
sb    = sort(med_b);                               % 99th percentile without the statistics toolbox
thr   = job.thresh*sb(max(1,round(0.99*numel(sb))));
clear sb
good  = med_b > thr;

r1 = zeros(N,1,'single');  r1(good) = med_c(good)./med_b(good);  r1(r1 >= 1) = 0;
mu = sum_all/n_all;  sd = sqrt(max(sq_all/n_all - mu.^2,0)*n_all/(n_all-1));
r2 = zeros(N,1,'single');  g2 = sd > 0 & good;  r2(g2) = single(mu(g2)./sd(g2));
r3 = zeros(N,1,'single');  r3(good) = med_s(good)./med_b(good);
% positive version: the range of the signed values (percentiles in the good voxels) is mapped to 0..pos_scale,
% so that the contrast uses the whole range (a fixed range of -1..1 would use only a small part of it)
sg   = sort(r3(good));
lo   = sg(max(1,round(job.pos_pct(1)/100*numel(sg))));
hi   = sg(min(numel(sg),max(1,round(job.pos_pct(2)/100*numel(sg)))));
clear sg
r4 = zeros(N,1,'single');  r4(good) = (min(max(r3(good),lo),hi)-lo)/(hi-lo)*job.pos_scale;

write_img(job.out.median,  r1,dim,Vref,'median(cbv)/median(bold)');
write_img(job.out.meanstd, r2,dim,Vref,'mean/std of cbv and bold');
write_img(job.out.signed,  r3,dim,Vref,'phase signed median(cbv)/median(bold)');
write_img(job.out.signedpos,r4,dim,Vref,'positive version of the phase signed image');
end

%% ------------------------------------------------------------------------
function [m,re,im] = resample_frame(Vm,Vp,W,dim,pscale)
% magnitude and real/imaginary part of one volume (magnitude + phase) in the reference grid
A    = spm_read_vols(Vm);
P    = spm_read_vols(Vp)*pscale;                   % radians
P    = single(P);  A = single(A);
C    = A.*cos(P);  S = A.*sin(P);
Q    = Vm.mat\W;                                   % voxel coordinates in this volume (header incl. sidecar)
grid = {1:size(A,1),1:size(A,2),1:size(A,3)};
m    = interp(grid,A,Q);
re   = interp(grid,C,Q);
im   = interp(grid,S,Q);
end

function v = interp(grid,A,Q)
F = griddedInterpolant(grid,A,'linear','none');
v = F(Q(1,:)',Q(2,:)',Q(3,:)');
v(isnan(v)) = 0;
v = single(v);
end

function write_img(fname,vec,dim,Vref,txt)
V = struct('fname',fname,'dim',dim,'dt',[spm_type('float32') 0],'mat',Vref.mat,'pinfo',[1;0;0], ...
    'descrip',['derived (sbp_vaso_T1): ' txt]);
spm_write_vol(V,reshape(vec,dim));
end
