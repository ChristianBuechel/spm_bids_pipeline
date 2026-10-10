function labels = classify_geometry(sig,acq_type)

% function labels = classify_geometry(sig,acq_type)
% decides for each volume of ONE series whether it is a 'brain' or a 'spinal'
% EPI, from its geometry. This is the core of classify_epi_acq (DICOM headers)
% and can also be used with converted NIfTIs (see nifti_signature).
%
% sig       n_vol x 3 matrix [n_slices, slice thickness, in-plane voxel size]
% acq_type  vars.acq_type from get_study_specs
%             {'brain'} (default if empty) or {'spinal'}: all volumes get this label
%             {'spinal','brain'} (either order): combined acquisition with two
%             geometries in one series. spinal = fewer slices AND larger
%             quotient (slice thickness / in-plane voxel size). If these two
%             criteria disagree (or there are not exactly two geometries) an
%             error is thrown rather than guessing.
% labels    cellstr, one label per volume

if nargin < 2 || isempty(acq_type)
    acq_type = {'brain'};
end
acq_type = cellstr(acq_type);
n_vol    = size(sig,1);

if numel(acq_type) == 1
    if ~ismember(acq_type{1},{'brain','spinal'})
        error('vars.acq_type has to be {''brain''}, {''spinal''} or {''spinal'',''brain''}');
    end
    labels = repmat(acq_type,n_vol,1);
    return
end
if ~(numel(acq_type) == 2 && all(ismember(acq_type,{'spinal','brain'})) && ~strcmp(acq_type{1},acq_type{2}))
    error('vars.acq_type has to be {''brain''}, {''spinal''} or {''spinal'',''brain''}');
end

[u,~,grp] = unique(round(sig,3),'rows');
if size(u,1) ~= 2
    error('Combined brain/spinal acquisition expected, but found %d different slice geometries',size(u,1));
end

quot = u(:,2)./u(:,3);
[~,fewer]  = min(u(:,1));
[~,larger] = max(quot);
if u(1,1) == u(2,1) || fewer ~= larger
    error('Cannot tell brain from spinal: slices %d/%d, thickness/in-plane %.2f/%.2f (spinal should have fewer slices AND larger quotient)',...
        u(1,1),u(2,1),quot(1),quot(2));
end

lab         = {'brain','brain'};
lab{fewer}  = 'spinal';
labels      = lab(grp)';
end
