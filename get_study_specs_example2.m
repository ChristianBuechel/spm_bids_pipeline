function [path,vars,analysis] = get_study_specs

%% path definitions
hostname = char(getHostName(java.net.InetAddress.getLocalHost));
switch hostname
    case 'Rainbow'
        path.baseDir     = 'd:\offhum_bids\';
        path.templateDir = 'd:\offhum_bids\templates\';
        vars.max_procs   = 12;

    case 'motown'
        path.baseDir     = 'c:\Users\buechel\Data\offhum_bids\';
        path.templateDir = 'c:\Users\buechel\Data\offhum_bids\templates\';                
        vars.max_procs   = 8;

end
path.rawDir          = fullfile(path.baseDir, 'rawdata'); % that should be the folder with the zipped not-to-touch data sets
path.derivDir        = fullfile(path.baseDir, 'derivatives');
path.preprocDir      = fullfile(path.baseDir, 'derivatives', 'spm_preprocessing');
path.firstlevelDir   = fullfile(path.baseDir, 'derivatives', 'spm_firstlevel');
path.secondlevelDir  = fullfile(path.baseDir, 'derivatives', 'spm_secondlevel');

%% vars definitions

% various predefined names (do not change)
vars.skullStripID    = 'skull-strip-T1.nii';
vars.T1maskID        = 'brain_mask.nii';
vars.templateID      = 'cb_Template_%d_Dartel.nii';
vars.templateT1ID    = 'cb_Template_T1.nii';
%vars.groupMaskID     = 'brainmask.nii';
vars.groupMaskID     = 'neuromorphometrics.nii';

%% this need to be adapted to your study / computer--------------
vars.task        = 'offhum';
vars.nRuns       = 3;
vars.nSess       = 1;
% get info for slice timing correction
vars.sliceTiming.so    = [1490.00 1420.00 1347.50 1277.50 1207.50 1135.00 1065.00 992.50 922.50 852.50 780.00 710.00 640.00 567.50 497.50 427.50 355.00 285.00 212.50 142.50 72.50 0.00...
    1490.00 1420.00 1347.50 1277.50 1207.50 1135.00 1065.00 992.50 922.50 852.50 780.00 710.00 640.00 567.50 497.50 427.50 355.00 285.00 212.50 142.50 72.50 0.00...
    1490.00 1420.00 1347.50 1277.50 1207.50 1135.00 1065.00 992.50 922.50 852.50 780.00 710.00 640.00 567.50 497.50 427.50 355.00 285.00 212.50 142.50 72.50 0.00]; % in ms
vars.sliceTiming.tr       = 1.58; % in s
vars.sliceTiming.nslices  = 66;
vars.sliceTiming.refslice = 750;



analysis.all_subs  = [4:10 11 12 14:19 21 22 24 26 27 29:35 37:39]; % very all;
single_group       = ones(size(analysis.all_subs));

analysis.group_ind  = single_group; %index 1
analysis.group_weights = [1];
analysis.group_names   = {'All'};

analysis.parallel          = 0;
analysis.noise_corr        = ['mov24_physio'];
%analysis.noise_corr        = ['physio'];
analysis.cvi               = 'none'; % 'AR(1)'  'FAST' 'none'
analysis.shift             = 0; % shift all onsets
analysis.skernel           = 6; % smoothing kernel
analysis.wls               = 0;
analysis.bs                = 1;
analysis.concatenate       = 0;

analysis.cond_names        = {'control_T1','control_T2','control_T3','offset_T1','offset_T2','offset_T3',...
    'control_rate_T1','control_rate_T2','control_rate_T3','offset_rate_T1','offset_rate_T2','offset_rate_T3'};

%what to do
analysis.do_model          = 0;
analysis.do_est            = 0;
analysis.do_vasa           = 0;
analysis.do_cons           = 0;
analysis.do_correct_vasa   = 0;
analysis.do_warp           = 0;
analysis.do_smooth         = 0;

analysis.do_fact           = 0;
analysis.do_fact_con       = 0;

analysis.do_one_t          = 0;


do_hrf = 1;
do_fir = 0;

if do_hrf
    analysis.max_procs         = 12;
    
    [analysis.t_con, analysis.t_con_names] = get_hrf_cons;
    
    analysis.ana               = 2; %hrf
    analysis.n_base            = 1;
    analysis.name_ana          = 'offset_hrf';
    
elseif do_fir
    analysis.max_procs         = 8; % big design matrix ...
    analysis.t_con             = [];
    analysis.t_con_names       = {};
    analysis.ana               = 1; %fir
    analysis.n_base            = 34;
    analysis.name_ana          = 'offset_fir';
end

    function [t_con, t_con_names] = get_hrf_cons
        
        t_con = [ 0  0 0 -1 2 -1 zeros(1,6);... %akt
            1 -2 1 -1 2 -1 zeros(1,6);...
            0  1 0  0 -1 0 zeros(1,6);...
            0  0 1  0 0 -1 zeros(1,6);...
            0  0 0  1 0 -1 zeros(1,6);...
            0  0 0  -1 1 0 zeros(1,6);...
            0  0 0   0 1 -1 zeros(1,6);...
            1 -1 0  -1 1 0 zeros(1,6);...
            1 0 -1  -1 0 1 zeros(1,6)];
        t_con_names         = {'offT2>T1_T3',...
            'intT2>T1_T3',...
            'conT2>offT2',...
            'conT3>offT3',...
            'offT3<T1',...
            'offT2>T1',...
            'offT2>T3',...
            'intT2>T1',...
            'intT3>T1'};
    end
end