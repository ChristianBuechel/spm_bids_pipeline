function [path,vars,analyses,import,qc] = get_study_specs
    %% path definitions
    path.baseDir     = '/home/keuter/projects/painact/data/bids';
    path.templateDir = '/home/keuter/scripts/spm_bids_pipeline/templates/';
    path.derivDir        = fullfile(path.baseDir, 'derivatives');
    path.preprocDir      = fullfile(path.baseDir, 'derivatives', 'spm_preprocessing');
    path.firstlevelDir   = fullfile(path.baseDir, 'derivatives', 'spm_firstlevel');
    path.secondlevelDir  = fullfile(path.baseDir, 'derivatives', 'spm_secondlevel');
    path.physioDir       = fullfile(path.baseDir, 'physio');
    
    vars.max_procs   = 20; % how many parallel processes
    vars.parallel    = 0;  % parallel processing?
    
    %% Import related stuff
    import.prisma        = '/home/keuter/projects/painact/data/bids/participants.tsv';
    import.user          = 'keuter';
    import.server        = 'revelations.nin.uke.uni-hamburg.de';
    
    hh = 1;
    import.data(hh).enddir        = 'func';
    import.data(hh).type       = 'bold';
    import.data(hh).seq        = 'ninEPI_bold_v12C, fmri '; %protocol name (trailing space makes it unique)
    import.data(hh).cond       = 'n > 180'; % heuristic to get only valid runs (e.g. more than 1000 volumes)
    hh = hh + 1;
    
    import.data(hh).dir        = 'anat'; % valid BIDS dir name
    import.data(hh).type       = 'T1w'; % valid BIDS file name
    import.data(hh).seq        = 'ninFLASH_v14A_df, mprage, defa-SAT_DEFA ';
    import.data(hh).cond       = 'n == 240'; % heuristic to get only valid runs (e.g. exactly 240 slices)
    hh = hh + 1;
    
    import.data(hh).dir        = 'fmap';
    import.data(hh).type       = 'magnitude';
    import.data(hh).seq        = 'gre_field_map, 3mm, filter M ';
    import.data(hh).cond       = 'n == 132';
    hh = hh + 1;
    
    import.data(hh).dir        = 'fmap';
    import.data(hh).type       = 'phasediff';
    import.data(hh).seq        = 'gre_field_map, 3mm, filter M ';
    import.data(hh).cond       = 'n == 66';
    hh = hh + 1;
    
    % Specific runs can be excluded even if they match the criteria above see
    % sbp_import_data.m for details on how to mark this in participants.tsv
    
    import.dummies            = 0; % these scans are removed when merging 3D epifiles to a 4D file
    
    %% vars definitions
    
    % various predefined names (change only if you know what you are doing)
    vars.skullStripID    = 'skull-strip-T1.nii';
    vars.T1maskID        = 'brain_mask.nii';
    vars.templateID      = 'cb_Template_%d_Dartel.nii';
    vars.templateT1ID    = 'cb_Template_T1.nii';
    vars.groupMaskID     = 'neuromorphometrics.nii';
    vars.brainstemID     = 'Cerebellum_SUIT_05.nii'; % or 'brainstem_mask.nii'
    vars.templateStripID = 'cb_Template_SkullStrip.nii';
    
    %% General settings
    
    vars.task        = 'painact';
    vars.nRuns       = 6; %
    vars.nSess       = 1;
    
    % get info for slice timing correction MUST BE in ms
    vars.sliceTiming.so       = [1915, 1855, 1795, 1735, 1675, 1615, 1555, 1495, 1437.5, 1377.5, 1317.5, 1257.5, 1197.5, 1137.5, 1077.5, 1017.5, 957.5, 897.5, 837.5, 777.5, 717.5, 657.5, 597.5, 540, 480, 420, 360, 300, 240, 180, 120, 60, 0, ...
                                 1915, 1855, 1795, 1735, 1675, 1615, 1555, 1495, 1437.5, 1377.5, 1317.5, 1257.5, 1197.5, 1137.5, 1077.5, 1017.5, 957.5, 897.5, 837.5, 777.5, 717.5, 657.5, 597.5, 540, 480, 420, 360, 300, 240, 180, 120, 60, 0];
    vars.sliceTiming.tr       = 1.9750; % in s
    vars.sliceTiming.nslices  = 66;
    vars.sliceTiming.refslice = vars.sliceTiming.tr / 2 * 1000;
    
    %% QC settings
    %file name for mov regressors
    qc.movFilename  = 'noise_mov_rasub-%02d_ses-%02d_task-%s_run-%02d_bold.mat';
    
    %thresholds for movement inspection
    qc.threshSpike     = 0.7; %threshold for spikes
    qc.threshMov       = 3; %threshold for overall movement within 1 run or between runs
    qc.percSpike       = 0.05; %threshold in % in number of volumes that are discarded because of spikes
    
    % display options
    
    qc.maxDisImg       = 7; %number of subject images displayed at once (not all numbers make sense here, 7 means a 4x2 display including template, otherwise try 5(3x2))
    qc.contour         = 1; %1: template contour is displayed on images; 0: is not displayed
    
    % first level related qc
    
    qc.tThresh         = 1.7; %t-value threshold, only values above this threshold will be shown
    qc.tMapSmoothK     = 2; %smoothing kernel for t-map
    qc.overlayColor    = [1 0 0]; %color for overlay (here red)
    
    %% Analyses

    %% Define analyses array
    % Initialize with base settings that are common across analyses
    base_analysis = struct(...
        'all_subs', 1:58, ...
        'exclude', struct('sub',48, 'ses',1, 'run', [6]), ...
        'noise_corr', 'mov24_wm_csf_roi_physio', ...
        'cvi', 'none', ...
        'shift', 0, ...
        'skernel', 6, ...
        'hpf', 120, ...
        'bs', 0, ...
        'use_vasa', 0, ...
        'sess', [1], ...
        'prune', 1, ...
        'do_model', 1, ...
        'do_est', 1, ...
        'do_vasa', 0, ...
        'do_cons', 1, ...
        'do_correct_vasa', 0, ...
        'do_warp', 1, ...
        'do_smooth', 1, ...
        'fact_dept', 0, ...
        'fact_var', 0, ...
        'do_fact', 0, ...
        'do_fact_con', 0, ...
        'do_one_t', 0 ...
    );

    % Remove excluded subjects
    excluded_subs = [16, 43, 45, 53, 56];
    base_analysis.all_subs(ismember(base_analysis.all_subs, excluded_subs)) = [];
    
    % Group settings
    base_analysis.group_ind = ones(size(base_analysis.all_subs));
    base_analysis.group_weights = [1];
    base_analysis.group_ana_names = {'All'};

    %% Define individual analyses
    % PPI Analysis - dpINS
    analyses(1) = base_analysis;
    analyses(1).do_ppi = 1;
    analyses(1).concatenate = 1;
    analyses(1).ana = 2;  % hrf
    analyses(1).n_base = 1;
    analyses(1).name_ana = 'painact_ppi_v3';
    analyses(1).events = '';
    analyses(1).cond_names = {'CUE_Pain-Move_L', 'CUE_Pain-noMove_L', ...};  % your conditions
    analyses(1).p_mod = cell(1, numel(analyses(1).cond_names));
    analyses(1).ppi = struct(...
        'def', 'sphere', ...
        'spec', 6, ...
        'xyz', [-39 -16.5 18]', ...
        'str', 'dpINS_L_test_6mm', ...
        'name', 'dpINS_L_test', ...
        'skernel', 6, ...
        'conds', 5:8, ...
        'deconvolute', true, ...
        'tcon', struct(...
            'vec', [0 0 0 0 1 1 -1 0 1 0 0 0 0 0 0 0 0 0 0 0 0], ...
            'name', 'STM_L_Pain>Nopain' ...
        ));

    % PPI Analysis - Cerebellum
    analyses(2) = base_analysis;
    analyses(2).do_ppi = 1;
    analyses(2).concatenate = 1;
    analyses(2).ana = 2;
    analyses(2).n_base = 1;
    analyses(2).name_ana = 'painact_ppi_v3';
    analyses(2).events = '';
    analyses(2).cond_names = analyses(1).cond_names;  % same conditions
    analyses(2).p_mod = analyses(1).p_mod;  % same p_mod
    analyses(2).ppi = struct(...
        'def', 'mask', ...
        'spec', '/home/keuter/projects/painact/masks/diedrichsen_mask.nii', ...
        'str', 'cerebellum', ...
        'name', 'cerebellum', ...
        // ...rest of cerebellum PPI settings...
    );

    % HRF Analysis (if needed)
    analyses(3) = base_analysis;
    analyses(3).concatenate = 1;
    analyses(3).ana = 2;
    analyses(3).n_base = 1;
    analyses(3).name_ana = 'painact_hrf_left2';
    // ...rest of HRF analysis settings...

    % Add additional analyses as needed...
end