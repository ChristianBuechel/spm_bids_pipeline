function [oXYZvox, oXYZmm, ib]  = transform_back(XYZmm, Vnormalized, Vorig, Nback, prune)

if nargin < 5
    prune = 1; %prune if not specified
end

if isempty(Nback)
    oXYZmm = XYZmm;
else
    bsplin = [2 2 2]; %order for b-spline interpolation  
    XYZvox           = inv(Vnormalized.mat)*[XYZmm; ones(1,size(XYZmm,2))];
    XYZvox           = XYZvox(1:3,:);
    oXYZmm           = zeros(size(XYZvox));
    % get coords in original space
    % let's do bspline interp
    for i=1:3
        c           = spm_bsplinc(squeeze(Nback.dat(:,:,:,1,i)),[bsplin  0 0 0]); % no wrap
        oXYZmm(i,:) = spm_bsplins(c,XYZvox(1,:),XYZvox(2,:),XYZvox(3,:),[bsplin  0 0 0]);
    end
end

oXYZvox          = inv(Vorig.mat)*[oXYZmm; ones(1,size(oXYZmm,2))];
oXYZvox          = (oXYZvox(1:3,:)); %so that unique 'rows' works

if (prune == 1)
    oXYZvox          = oXYZvox'; %so that unique 'rows' works
    oXYZmm           = oXYZmm';
    XYZvox           = XYZvox';
    XYZmm            = XYZmm';
    [oXYZvox, ia, ib] = unique(round(oXYZvox),'stable','rows'); %only get rounded voxel coordinates
    oXYZvox          = oXYZvox'; %and back
    oXYZmm           = oXYZmm(ia,:)';
    %XYZvox           = XYZvox(ia,:)';
    %XYZmm            = XYZmm(ia,:)';
end