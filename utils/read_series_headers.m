function [subtype,ref] = read_series_headers(all_dicom,acq_type)

% function [subtype,ref] = read_series_headers(all_dicom,acq_type)
% subtype ('mag' or 'phase') from the first DICOM of a functional series, and one
% header per label (brain/spinal) in ref.brain / ref.spinal (used for the task
% json). In a combined acquisition the volumes of both geometries alternate, so
% the first ones are classified until both have been seen.
%   all_dicom  char matrix of DICOM file names (from spm_select)
%   acq_type   vars.acq_type from get_study_specs

if isempty(all_dicom)
    error('no DICOM files (MR.*) available');
end
n_read = min(size(all_dicom,1),40);
dc     = spm_dicom_headers(all_dicom(1:n_read,:));
if contains(dc{1}.ImageType,'ORIGINAL\PRIMARY\M\')
    subtype = 'mag';
elseif contains(dc{1}.ImageType,'ORIGINAL\PRIMARY\P\')
    subtype = 'phase';
else
    error('unknown image type %s (neither MAGNITUDE, nor PHASE)',dc{1}.ImageType);
end
labels = classify_epi_acq(dc,acq_type);
ref    = struct();
for ii = 1:numel(labels)
    if ~isfield(ref,labels{ii})
        ref.(labels{ii}) = dc{ii};
    end
end
end
