function [oXYZvox, oXYZmm]  = transform_back(XYZmm, Vnormalized, Vorig, Nback)

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

%for i_vox = 1:size(XYZvox,2)
%    oXYZmm(:,i_vox) = squeeze(Nback.dat(XYZvox(1,i_vox),XYZvox(2,i_vox),XYZvox(3,i_vox),:,:));
%end

oXYZvox          = inv(Vorig.mat)*[oXYZmm; ones(1,size(oXYZmm,2))];
oXYZvox          = (oXYZvox(1:3,:))'; %so that unique 'rows' works
oXYZmm           = oXYZmm';
XYZvox           = XYZvox';
[oXYZvox, ia, ~] = unique(round(oXYZvox),'stable','rows'); %only get rounded voxel coordinates
oXYZvox          = oXYZvox'; %and back
oXYZmm           = oXYZmm(ia,:)';
XYZvox           = XYZvox(ia,:)';
