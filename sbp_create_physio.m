function sbp_create_physio(all_sub_ids)
% this is a wrapper for the three functions needed to create physiological
% noise regressors.

convert_and_clean_physio(all_sub_ids);
create_tapas(all_sub_ids)

end
