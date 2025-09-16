function run_multiple_ppis
  
all_sub_ids = [1];

addpath(genpath('/home/keuter/scripts/spm_bids_pipeline'));

[~, ~, analysis] = get_study_specs();
    
    for i = 1:numel(analysis.ppi)
        fprintf('\nProcessing PPI %d/%d: %s\n', i, numel(analysis.ppi), analysis.ppi(i).name);
        sbp_analyses(all_sub_ids, i);
    end
end