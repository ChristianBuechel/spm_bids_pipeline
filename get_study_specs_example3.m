function [path,vars,analysis] = get_study_specs

%% path definitions



path.baseDir     = 'd:\fearamy_bids\';
path.templateDir = 'd:\fearamy_bids\spm_bids_pipeline\templates\';

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
vars.groupMaskID     = 'neuromorphometrics.nii';
%% this need to be adapted to your study / computer--------------
vars.max_procs   = 8;
vars.task        = 'fearamy';
vars.nRuns       = 1;
vars.nSess       = 1;
% get info for slice timing correction
vars.sliceTiming.so       = [913.85, 837.69, 761.54, 685.38, 609.23, 533.08, 456.92, 380.77, 304.62,...
                             228.46, 152.31, 76.15, 0.00, 913.85, 837.69, 761.54, 685.38, 609.23, 533.08,...
                             456.92, 380.77, 304.62, 228.46, 152.31, 76.15, 0.00]; % in ms
vars.sliceTiming.tr       = 0.99; % in s
vars.sliceTiming.nslices  = 26;
vars.sliceTiming.refslice = 500;


analysis.all_subs  = [5:44]; % very all;

end