function [path,vars,analysis,import,qc] = get_study_specs_new

    %% path definitions
    path.baseDir     = '/home/keuter/projects/painact/data/bids';
    path.templateDir = '/home/keuter/scripts/spm_bids_pipeline/templates/';
    path.derivDir        = fullfile(path.baseDir, 'derivatives');
    path.preprocDir      = fullfile(path.baseDir, 'derivatives', 'spm_preprocessing');
    path.firstlevelDir   = fullfile(path.baseDir, 'derivatives', 'spm_firstlevel');
    path.secondlevelDir  = fullfile(path.baseDir, 'derivatives', 'spm_secondlevel');
    path.physioDir       = fullfile(path.baseDir, 'physio');
    
    vars.max_procs   = 20; % how many parallel processes
    vars.parallel    = 1;  % parallel processing?
    
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
    
    
    %% Analysis settings
    % the next section can be used to define different groups for 2nd level analysis
    % We first define base settings that are typically used across analyses
    % These can be overwritten below for specific analyses

    %whatto do
    do_hrf       = 0;
    do_hrf_param = 0;
    do_fir_cue   = 0;
    do_fir_stim  = 0;
    do_ppi       = 1;

    idx = 0; % analysis counter

    if do_hrf    

        idx = idx + 1;
        
        
        analysis(idx) = get_analysis_defaults();

        analysis(idx).concatenate       = 1; % concatenate or not ?
        analysis(idx).ana               = 2; % hrf
        analysis(idx).n_base            = 1;
        analysis(idx).name_ana          = 'painact_hrf_left2';
        analysis(idx).events            = ''; % this refers to the tsv file: e.g. sub-08_ses-01_task-lcpa_run-01_events_param.tsv
        % analysis.cond_names        = {'CUE_Pain-Move_L', 'CUE_Pain-noMove_L', 'CUE_noPain-Move_L', 'CUE_noPain-noMove_L', ...
        %                               'STM_Pain-Move_L', 'STM_Pain-noMove_L', 'STM_noPain-Move_L', 'STM_noPain-noMove_L', ...    
        %                               'CUE_Pain-Move_R', 'CUE_Pain-noMove_R', 'CUE_noPain-Move_R', 'CUE_noPain-noMove_R', ...
        %                               'STM_Pain-Move_R', 'STM_Pain-noMove_R', 'STM_noPain-Move_R', 'STM_noPain-noMove_R', ...  
        %                               'CUE_Rating', 'STM_Rating', 'Rating', 'CUE_invalid', 'STM_invalid'};
        analysis(idx).cond_names        = {'CUE_Pain-Move_L', 'CUE_Pain-noMove_L', 'CUE_noPain-Move_L', 'CUE_noPain-noMove_L', ...
                                           'STM_Pain-Move_L', 'STM_Pain-noMove_L', 'STM_noPain-Move_L', 'STM_noPain-noMove_L', ...    
                                            'CUE_Rating', 'STM_Rating', 'Rating', 'CUE_invalid', 'STM_invalid'};
        %'CUE_invalid', 'STM_invalid'
        analysis(idx).p_mod             = cell(1, numel(analysis.cond_names)); % vis and pain have 1 or 2 parametric modulators
        
    
        [analysis(idx).t_con, analysis(idx).t_con_names] = get_T_contrasts_hrf_left(analysis.cond_names);
        [analysis(idx).f_con, analysis(idx).f_con_names] = get_F_contrasts_hrf_left(analysis.cond_names);

    
    end
    
    if do_hrf_param

        idx = idx + 1;
        analysis(idx) = get_analysis_defaults();

        analysis(idx).concatenate       = 1; % concatenate or not ?
        analysis(idx).ana               = 2; % hrf
        analysis(idx).n_base            = 1;
        analysis(idx).name_ana          = 'painact_hrf_param';
        analysis(idx).events            = '_param'; % this refers to the tsv file: e.g. sub-08_ses-01_task-lcpa_run-01_events_param.tsv
        % cave: hyphens in event names dont work with pipeline, so I
        % updated the event files (not done for normal hrf analysis)
        analysis(idx).cond_names        = {'CUE_Pain_Move_L', 'CUE_Pain_noMove_L', 'CUE_noPain_Move_L', 'CUE_noPain_noMove_L', ...
                                      'STM_Pain_Move_L', 'STM_Pain_noMove_L', 'STM_noPain_Move_L', 'STM_noPain_noMove_L', ...    
                                      'CUE_Pain_Move_R', 'CUE_Pain_noMove_R', 'CUE_noPain_Move_R', 'CUE_noPain_noMove_R', ...
                                      'STM_Pain_Move_R', 'STM_Pain_noMove_R', 'STM_noPain_Move_R', 'STM_noPain_noMove_R', ...  
                                      'CUE_Rating', 'STM_Rating', 'Rating', 'CUE_invalid', 'STM_invalid'};
    
        analysis(idx).p_mod               = {{}, {}, {}, {}, ...
                                        {'STM_Pain_Move_L_param'}, {}, {'STM_noPain_Move_L_param'}, {}, ...    
                                        {}, {}, {}, {}, ...
                                        {}, {}, {}, {}, ...  
                                        {}, {}, {}, {},{}};
    
        [analysis(idx).t_con, analysis(idx).t_con_names] = get_T_contrasts_hrf_param;
        [analysis(idx).f_con, analysis(idx).f_con_names] = get_F_contrasts_hrf_param;

    end

    if do_fir_cue
        idx = idx + 1;
        analysis(idx) = get_analysis_defaults();

        analysis(idx).max_procs         = 12; % large DesMtx
        analysis(idx).concatenate       = 1;
        analysis(idx).t_con             = [];
        analysis(idx).t_con_names       = [];
        analysis(idx).ana               = 1; % fir
        analysis(idx).n_base            = 8;
        analysis(idx).bin_size          = 2; % in seconds
        analysis(idx).name_ana          = 'painact_fir_cue';
        analysis(idx).events            = ''; % this refers to the tsv file: e.g. sub-08_ses-01_task-lcpa_run-01_bold_param.tsv
        analysis(idx).cond_names        = {'CUE_Pain-Move_L', 'CUE_Pain-noMove_L', 'CUE_noPain-Move_L', 'CUE_noPain-noMove_L', ...
                                      'CUE_Pain-Move_R', 'CUE_Pain-noMove_R', 'CUE_noPain-Move_R', 'CUE_noPain-noMove_R'}; % only cue
        analysis(idx).p_mod             = cell(1, numel(analysis.cond_names));
        [analysis(idx).f_con, analysis(idx).f_con_names] = get_F_contrasts_fir(analysis.cond_names);
    end
    
    if do_fir_stim

        idx = idx + 1;
        analysis(idx) = get_analysis_defaults();

        analysis(idx).max_procs         = 12; % large DesMtx
        analysis(idx).concatenate       = 1;
        analysis(idx).t_con             = [];
        analysis(idx).t_con_names       = [];
        analysis(idx).ana               = 1; % fir
        analysis(idx).n_base            = 8;
        analysis(idx).bin_size          = 2; % in seconds
        analysis(idx).name_ana          = 'painact_fir_stim';
        analysis(idx).events            = ''; % this refers to the tsv file: e.g. sub-08_ses-01_task-lcpa_run-01_bold_param.tsv
        analysis(idx).cond_names        = {'STM_Pain-Move_L', 'STM_Pain-noMove_L', 'STM_noPain-Move_L', 'STM_noPain-noMove_L', ...
                                      'STM_Pain-Move_R', 'STM_Pain-noMove_R', 'STM_noPain-Move_R', 'STM_noPain-noMove_R'}; % only stim
        analysis(idx).p_mod             = cell(1, numel(analysis.cond_names));
        [analysis(idx).f_con, analysis(idx).f_con_names] = get_F_contrasts_fir(analysis.cond_names);
    end


    if do_ppi      

        % Get PPI specifications
        ppi_specs = get_ppi_specs();
        
        % Store base model settings for PPI
        ppi_base_settings             = get_analysis_defaults();
        ppi_base_settings.do_ppi      = 1;
        ppi_base_settings.concatenate = 1;
        ppi_base_settings.ana         = 2; % hrf
        ppi_base_settings.name_ana    = 'painact_ppi';
        ppi_base_settings.n_base      = 1;
        ppi_base_settings.events      =   '';
        ppi_base_settings.cond_names  = {'CUE_Pain-Move_L', 'CUE_Pain-noMove_L', 'CUE_noPain-Move_L', 'CUE_noPain-noMove_L', ...
                                        'STM_Pain-Move_L', 'STM_Pain-noMove_L', 'STM_noPain-Move_L', 'STM_noPain-noMove_L', ... 
                                        'CUE_Pain-Move_R', 'CUE_Pain-noMove_R', 'CUE_noPain-Move_R', 'CUE_noPain-noMove_R', ...
                                        'STM_Pain-Move_R', 'STM_Pain-noMove_R', 'STM_noPain-Move_R', 'STM_noPain-noMove_R', ...  
                                        'CUE_Rating', 'STM_Rating', 'Rating', 'CUE_invalid', 'STM_invalid'};
        ppi_base_settings.p_mod       = cell(1, numel(ppi_base_settings.cond_names));
        
        % PPI specific contrasts
        % Convention for PPI contrasts: [PSY_regs SEED PPI_regs]
        %
        % PSY_regs: psychological main effects (n columns = # conditions in ppi_specs.conds)
        % SEED:     seed region main effect (1 column)
        % PPI_regs: PPI interaction terms (n columns = # conditions)
        ppi_base_settings.t_con       = [0 0 0 0 0 1 1 -1 -1;
                                         0 0 0 0 0 1 -1 1 -1;
                                         0 0 0 0 0 -1 1 -1 1];
        ppi_base_settings.t_con_names = {'ppi_pain>nopain', 'ppi_move>nomove', 'ppi_move<nomove'};
       


        for ppi_idx = 1:length(ppi_specs)
            idx = idx + 1;

            % Create a complete copy of the base settings
            ppi_analysis = ppi_base_settings;
            
            % Add the region-specific name and PPI specification
            ppi_analysis.ppi = ppi_specs(ppi_idx);

            analysis(idx) = ppi_analysis;

        end
            % Create one analysis per PPI specification
        % for ppi_idx = 1:length(ppi_specs)
        %     idx = idx + 1;
        % 
        %     % For the first PPI, set the complete structure
        %     if ppi_idx == 1
        %         analysis(idx) = ppi_base_settings;
        %         analysis(idx).ppi = ppi_specs(ppi_idx);
        %     else
        %         % For subsequent PPIs, set fields individually to maintain structure compatibility
        %         % First, copy all fields from the previous PPI analysis
        %         analysis(idx) = analysis(idx-1);
        % 
        %         % Then update only the fields that need to change
        %         analysis(idx).ppi = ppi_specs(ppi_idx);
        %     end
        % 
        % end
     
    end


end

function defaults = get_analysis_defaults()
    
    % Returns a struct containing all default analysis settings
    % these parameters tend to be the same across analyses.
    % They can be overwritten by modifying the struct fields after calling this function.

    defaults = struct();
    
    % Subject selection
    defaults.all_subs = 1:58;
    excluded_subs = [16, 43, 45, 53, 56];
    defaults.all_subs(ismember(defaults.all_subs, excluded_subs)) = [];
    
    % Group settings
    single_group = ones(size(defaults.all_subs));
    defaults.group_ind = single_group;
    defaults.group_weights = [1];
    defaults.group_ana_names = {'All'};
    
    % Specific exclusions
    defaults.exclude(1) = struct('sub', 48, 'ses', 1, 'run', [6]);
    
    % Model settings
    defaults.noise_corr = 'mov24_wm_csf_roi_physio';
    defaults.cvi = 'none';
    defaults.shift = 0;
    defaults.skernel = 6;
    defaults.hpf = 120;
    defaults.bs = 0;
    defaults.use_vasa = 0;
    
    % Session settings
    defaults.sess = [1];
    defaults.prune = 1;
    
    % First level settings
    defaults.do_model = 1;
    defaults.do_est = 1;
    defaults.do_vasa = 0;
    defaults.do_cons = 1;
    defaults.do_correct_vasa = 0;
    defaults.do_warp = 1;
    defaults.do_smooth = 1;
    
    % Second level settings
    defaults.fact_dept = 0;
    defaults.fact_var = 0;
    defaults.do_fact = 0;
    defaults.do_fact_con = 0;
    defaults.do_one_t = 0;
end

function ppi_specs = get_ppi_specs()
    % Define all PPI specifications

     
    % PPI SETTINGS
    % analysis.ppi - Struct array defining psychophysiological interaction analyses
    
    % Fields:
    %   def         - ROI definition type: 'sphere', 'mask', 'box', 'cluster', 'all'
    %   spec        - ROI parameters depending on definition type:
    %                 - For 'sphere'/'box': radius in mm 
    %                 - For 'mask': full path to mask file
    %   xyz         - [x y z]' coordinates in mm for sphere/box ROI centers (only for def='sphere'/'box')
    %   str         - String identifier used for analysis directory naming
    %   name        - ROI name used for PPI regressor labeling
    %   skernel     - Smoothing kernel in mm applied to PPI seed time series
    %   conds       - Vector of condition indices to include in PPI (references analysis.cond_names)
    %   deconvolute - Boolean flag for HRF deconvolution of BOLD signal
    %   tcon        - Optional struct for two-pass ROI optimization:
    %                 .vec  - T-contrast vector to find peak activation within ROI
    %                 .name - Name of the contrast for documentation
    %
    % Example usage:
    %   analysis.ppi(1) = struct('def', 'sphere', 'spec', 6, 'xyz', [-39 -16.5 18]',...
    %                           'str', 'dpINS_L_6mm', 'name', 'dpINS_L',...
    %                           'skernel', 6, 'conds', 5:8, 'deconvolute', true,...
    %   
    
    % First PPI: dpINS_L
    ppi_specs(1) = struct( ...
            'def',         'sphere',    ...    
            'spec',        6,           ...    
            'xyz',         [-39 -16.5 18]', ...
            'str',         'dpINS_L_test_6mm', ...  
            'name',        'dpINS_L_test', ...       
            'skernel',     0,           ...    
            'conds',       5:8,         ...    
            'deconvolute', true,        ...                                     
            'tcon',        struct(...
                'vec',     [0 0 0 0 1 1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0], ...
                'name',    'STM_L_Pain>Nopain' ...
            ));
    
    % Second PPI: cerebellum
    ppi_specs(2) = struct( ...
            'def',         'mask',    ...   
            'spec',        '/home/keuter/projects/painact/masks/diedrichsen_mask.nii', ...   
            'xyz',         NaN, ...
            'str',         'cerebellum', ...   
            'name',        'cerebellum', ...       
            'skernel',     0,           ...    
            'conds',       5:8,         ...    
            'deconvolute', true,        ...                                     
            'tcon',        struct(...
                'vec',     [0 0 0 0 1 1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0], ...
                'name',    'STM_L_Pain>Nopain' ...
            ));
end