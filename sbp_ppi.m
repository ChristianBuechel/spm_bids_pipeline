function sbp_ppi(SPM, n_cond, ppi_specs, epi_2_templ, wm_file)

% Creates and estimates a PPI model based on an existing first-level analysis
% 1. Extracts PPI regressors
% 2. Creates new design matrix with PPI terms
% 3. Estimates new model in separate directory

if ~isstruct(SPM)
    spm_dir = spm_file(SPM,'fpath');
    load(fullfile(spm_dir,'SPM.mat'), 'SPM');
    SPM.swd = spm_dir;
end

%% Setup and extract PPI regressors
% xY and Uu are inputs to get_roi_ts_fc / g_ppi

xY = ppi_specs;
xY.out = Inf;
xY.Ic = 1;

% gPPI input
Uu = arrayfun(@(x) [x 1 1], ppi_specs.conds, 'UniformOutput', false);

% Extract ROI timeseries and create PPI terms
ppi = get_roi_ts_fc(SPM, char(wm_file), char(epi_2_templ), ppi_specs.skernel, xY, Uu);

%% Create new design matrix
% Includes the PPI terms and the seed regressor into the design matrix 
% and removes the experimental conditions (n_cond) from the design matrix

% save original condition names
cond_names = extractBetween(SPM.xX.name(ppi_specs.conds), " ", "*");

n_ppis = numel(ppi.xY.PPI); 

% extract regressors
[n_scans, ~]   = size(ppi.xY.PPI{1}.P);
psy_regressors = zeros(n_scans, n_ppis);
ppi_regressors = zeros(n_scans, n_ppis);

% Extract normalized regressors
for i = 1:n_ppis
    psy_regressors(:,i) = normit(ppi.xY.PPI{i}.P);
    ppi_regressors(:,i) = normit(ppi.xY.PPI{i}.o_ppi);
end

% create final design matrix
seed_ts        = normit(ppi.xY.PPI{1}.Y);         
new_regressors = [psy_regressors seed_ts ppi_regressors];
retain_idx     = n_cond+1:size(SPM.xX.X,2); % retain all regressors after the experimental conditions
SPM.xX.X       = [new_regressors, SPM.xX.X(:,retain_idx)];

% Update regressor names

roi_prefix     = ppi_specs.name;
psych_names    = strcat('PSY_', cond_names);
ppi_names      = strcat('PPI_', cond_names);
phys_name      = {sprintf('Y_%s', roi_prefix)};
SPM.xX.name    = [cellstr(psych_names), phys_name, cellstr(ppi_names), SPM.xX.name(retain_idx)];


%% Setup new analysis directory

temp_dir = SPM.swd; % this will be removed at the end
[base_dir, ana_dir] = fileparts(SPM.swd);
ppi_dir  = fullfile(base_dir, sprintf('PPI_%s_%s', ppi_specs.str, erase(ana_dir, '_tmp')));


% to avoid overwriting warnings 
if exist(ppi_dir, 'dir')
    rmdir(ppi_dir, 's')
end
mkdir(ppi_dir);

% Save and cleanup
SPM.swd  = ppi_dir;
save(fullfile(ppi_dir, 'PPI.mat'), "ppi")
rmdir(temp_dir, 's');

%% Estimate the model

spm_spm(SPM);

end


function Y = normit(X)
%normalise X to zero mean, std unity
Y = X - ones(size(X,1),1)*mean(X); %zero mean
Y = Y ./ (ones(size(Y,1),1)*std(Y)); %std 1
end