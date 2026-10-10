function import_convert(dicom_dir)

% function import_convert(dicom_dir)
% worker for the SPM batch system (cfg_basicio call_matlab) used by sbp_import_data:
% converts all DICOM files (MR.*) of one series to nifti (+ json) in the same folder.
% Does what the SPM batch module "DICOM Import" does (flat, nii, with meta data),
% but the file list is only needed when the job runs, i.e. after the download.

all_dicom = spm_select('FPlist',dicom_dir,'^MR.*');
if isempty(all_dicom)
    error('import_convert: no DICOM files (MR.*) in %s',dicom_dir);
end
hdr = spm_dicom_headers(all_dicom,false); % false: all fields (needed for the json files)
spm_dicom_convert(hdr,'all','flat','nii',dicom_dir,1);
end
