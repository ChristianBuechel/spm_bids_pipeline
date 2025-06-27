function sbp_create_physio(all_sub_ids)
% this is a wrapper for the three functions needed to create physiological
% noise regressors.
[path,~,~] = get_study_specs;
addpath([path.baseDir filesep 'spm_bids_pipeline' filesep 'utils'])

convert_physio_smr(all_sub_ids);
clean_physio(all_sub_ids);
create_tapas(all_sub_ids)

end
