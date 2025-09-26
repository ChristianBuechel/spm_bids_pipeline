# Changelog 2025-09-26

## `sbp_analyses.m`

- **Model Reuse**: Improved handling of existing PPI base models (no respecification or re-estimation required)
- **Verbosity**: More verbose printing of PPI settings 
- **PPI Base Model Consistency Check**: Added `check_ppi_consistency()` function to prevent inconsistent batches when some subjects have base models and others don't
- **User Interaction**: Improved directory cleanup dialog and PPI status reporting. This is especially important for potentially corrupted analyses. If in an exisiting temp dir the event files are missing or changed, the user is asked to delete the directories of all participants. Alternatively, they can opt for manually deleting the directories.
- **Multiple Analysis Support**: Added support for running analyses from struct array using the `idx` parameter

## `get_study_specs.m`

- **Multiple Analysis Definition**: Refactored to return analysis as a struct array instead of single struct, allowing multiple analyses to be run sequentially. This includes a major overhaul of the get_study_specs structure

## `sbp_ppi.m` &`get_roi_ts_fc.m` & `spm_pep_gppi.m`

-  Added deconvolution flag for time series extraction and PPI term calculation

