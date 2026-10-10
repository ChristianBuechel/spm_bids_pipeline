function labels = classify_epi_acq(hdr,acq_type)

% function labels = classify_epi_acq(hdr,acq_type)
% decides for each volume (DICOM header from spm_dicom_headers, i.e. one cell
% per volume as returned by that function) of ONE series whether it is a
% 'brain' or a 'spinal' EPI. The label is returned as cellstr (one per volume).
% The decision itself is made by classify_geometry (see there for the rules and
% for the meaning of acq_type = vars.acq_type).

if nargin < 2 || isempty(acq_type)
    acq_type = {'brain'};
end
acq_type = cellstr(acq_type);
n_vol    = numel(hdr);

if numel(acq_type) == 1
    labels = classify_geometry(zeros(n_vol,3),acq_type);
    return
end

sig = nan(n_vol,3);
for ii = 1:n_vol
    h = hdr{ii};
    sig(ii,:) = [csa_num(h,'NumberOfImagesInMosaic') h.SliceThickness h.PixelSpacing(1)];
end
labels = classify_geometry(sig,acq_type);
end

function v = csa_num(h,name)
% numeric value of an item of the CSA image header
ix = find(strcmp({h.CSAImageHeaderInfo.name},name),1);
if isempty(ix) || isempty(h.CSAImageHeaderInfo(ix).item)
    error('CSA image header has no entry %s (not a Siemens mosaic EPI?)',name);
end
v = str2double(h.CSAImageHeaderInfo(ix).item(1).val);
end
