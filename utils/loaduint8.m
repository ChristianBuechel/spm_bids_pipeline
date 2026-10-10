%==========================================================================
% function udat = loaduint8(V)
%==========================================================================
function udat = loaduint8(V)
% Load data from file indicated by V into an array of unsigned bytes.
if size(V.pinfo,2)==1 && V.pinfo(1) == 2
    mx = 255*V.pinfo(1) + V.pinfo(2);
    mn = V.pinfo(2);
else
    spm_progress_bar('Init',V.dim(3),...
        ['Computing max/min of ' spm_file(V.fname,'filename')],...
        'Planes complete');
    mx = -Inf; mn =  Inf;
    for p=1:V.dim(3)
        img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
        img = img(isfinite(img));
        mx  = max([max(img(:))+paccuracy(V,p) mx]);
        mn  = min([min(img(:)) mn]);
        spm_progress_bar('Set',p);
    end
end

% Another pass to find a maximum that allows a few hot-spots in the data.
spm_progress_bar('Init',V.dim(3),...
        ['2nd pass max/min of ' spm_file(V.fname,'filename')],...
        'Planes complete');
nh = 2048;
h  = zeros(nh,1);
for p=1:V.dim(3)
    img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
    img = img(isfinite(img));
    img = round((img+((mx-mn)/(nh-1)-mn))*((nh-1)/(mx-mn)));
    h   = h + accumarray(img,1,[nh 1]);
    spm_progress_bar('Set',p);
end
tmp = [find(cumsum(h)/sum(h)>0.9999); nh];
mx  = (mn*nh-mx+tmp(1)*(mx-mn))/(nh-1);

% Load data from file indicated by V into an array of unsigned bytes.
spm_progress_bar('Init',V.dim(3),...
    ['Loading ' spm_file(V.fname,'filename')],...
    'Planes loaded');
udat = zeros(V.dim,'uint8');

% Done the old way for compatibility with MATLAB versions older than R2011a
warning('off','MATLAB:RandStream:ActivatingLegacyGenerators')
st = rand('state');
rand('state',100);
% st = rng;           % Save old state
% rng(100,'twister'); % Replicable random numbers

for p=1:V.dim(3)
    img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
    acc = paccuracy(V,p);
    if acc==0
        udat(:,:,p) = uint8(max(min(round((img-mn)*(255/(mx-mn))),255),0));
    else
        % Add random numbers before rounding to reduce aliasing artifact
        r = rand(size(img))*acc;
        udat(:,:,p) = uint8(max(min(round((img+r-mn)*(255/(mx-mn))),255),0));
    end
    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear');
rand('state',st);
warning('on','MATLAB:RandStream:ActivatingLegacyGenerators')
% rng(st); % Return to old state

function acc = paccuracy(V,p)
if ~spm_type(V.dt(1),'intt')
    acc = 0;
else
    if size(V.pinfo,2)==1
        acc = abs(V.pinfo(1,1));
    else
        acc = abs(V.pinfo(1,p));
    end
end

