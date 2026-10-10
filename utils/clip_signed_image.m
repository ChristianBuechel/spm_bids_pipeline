function clip_signed_image(in,out,pct)

% function clip_signed_image(in,out,pct)
% writes a copy of an image in which the values are clipped to the percentile range pct
% ([lo hi], default [0.5 99.5]) of its non-zero voxels. Voxels that are exactly 0 stay 0
% (SPM's segmentation / normalisation ignores them, see spm_preproc8).
% Used for the signed VASO T1-like image (sbp_vaso_T1): the extreme ratio values (e.g. -1.8 .. 2.1
% while 99% of the voxels lie in -0.3 .. 0.5) stretch the intensity histogram range of SPM's
% normalisation, so that the tissue contrast is described by only a few bins.
%   in, out   file names (3D image, float output)

if nargin < 3 || isempty(pct), pct = [0.5 99.5]; end
V = spm_vol(in);
Y = spm_read_vols(V);
nz = Y ~= 0 & isfinite(Y);
x  = sort(Y(nz));                                  % percentiles without the statistics toolbox
lo = x(max(1,round(pct(1)/100*numel(x))));
hi = x(min(numel(x),max(1,round(pct(2)/100*numel(x)))));
Y(nz) = min(max(Y(nz),lo),hi);
fprintf('clip_signed_image: %s clipped to %.3f .. %.3f (range was %.3f .. %.3f)\n',spm_file(in,'filename'),lo,hi,x(1),x(end));

W = struct('fname',out,'dim',V.dim,'dt',[spm_type('float32') 0],'mat',V.mat,'pinfo',[1;0;0], ...
    'descrip',sprintf('derived (clip_signed_image): clipped to %g..%g percentile, %s',pct(1),pct(2),V.descrip));
spm_write_vol(W,Y);
end
