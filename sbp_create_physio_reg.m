function sbp_create_physio_reg(subIDs)
%Create and Save Physio Regressors

[path,vars]    = get_study_specs;
physioDir      = fullfile(path.baseDir,'physio');
nSubs          = length(subIDs);
nSess          = vars.nSess;
nRuns          = vars.nRuns;
taskName       = vars.task;
TR             = vars.sliceTiming.tr;
BIDS           = spm_BIDS(path.preprocDir);

for sub = 1:nSubs
    for ses = 1:nSess
        fprintf('\nCreating physio regressors for sub%02d ses%02d\n',subIDs(sub),ses);
        for run = 1:nRuns
            fprintf('Run%d ...\n', run);
            epi = spm_BIDS(BIDS,'data','sub',sprintf('%02d',subIDs(sub)),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',taskName,'type','bold');
            nScans = length(spm_select('ExtFPList',spm_file(epi{1},'path'),spm_file(epi{1},'filename'),Inf));

            % Get tsv data
            fileName = fullfile(spm_file(epi{1},'path'),sprintf('sub-%02d_ses-%02d_task-%s_run-%02d_physio',subIDs(sub),ses,taskName,run));
            try
                tsv = spm_load([fileName '.tsv.gz']);
            catch
                fprintf('No tsv physio file found for sub %02d ses %02d. Skipping ...\n',subIDs(sub),ses);
                continue;
            end

            %% Get sample rate from json header
            fid  = fopen(spm_file(fileName,'ext','json'));
            raw  = fread(fid,inf); % Reading the contents
            str  = char(raw');
            info = jsondecode(str);
            samp_int = 1/info.SamplingFrequency;
            fclose(fid); % closing the file

            [physio,breaths,beats] = get_physio(tsv,samp_int);

            if size(physio,1) ~= nScans
                warning('Physio regressors do not have the same length as number of scans! Please check why before continuing');
                return
            end

            %% save physio regs
            physioFile = fullfile(spm_file(epi{1},'path'),spm_file(spm_file(epi{1},'filename'),'prefix','physio_','ext','.mat')); % Name of the physio regressor output file
            save(physioFile,'physio');
            fprintf('finished and saved!\n');

            %save beats/breaths per minute for all subs for later quality control
            bpm.cardiac(sub,run) = beats/(size(physio,1)*TR)*60;
            bpm.resp(sub,run) = breaths/(size(physio,1)*TR)*60;
        end
    end
end
save(fullfile(physioDir,'physio_bpm_all.mat'),'bpm');
end
