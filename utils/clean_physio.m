function clean_physio(subIDs)

[path,vars,~,import]    = get_study_specs;
nSubs                   = length(subIDs);
BIDS                    = spm_BIDS(path.preprocDir);
physioDir               = fullfile(path.baseDir,'physio');
matDir                  = fullfile(physioDir,'mat');
nSess                   = vars.nSess;
nRuns                   = vars.nRuns;
taskName                = vars.task;
nDummies                = import.dummies;
TR                      = vars.sliceTiming.tr;
incomplRuns             = str2double(extract(import.data(1).cond,digitsPattern)); %rule for number of scans provided in get_study_specs for data import to find incomplete runs

for sub = 1:nSubs
    for ses = 1:nSess
        fileName = fullfile(matDir,sprintf('sub-%02d/ses-%02d',subIDs(sub),ses),sprintf('sub-%02d_ses-%02d_physio.mat',subIDs(sub),ses));

        %% First, load the mat file
        try
            data = load(fileName);
            do_preproc = 1;
            fprintf('\n\nRunning physio preprocessing for sub%02d session %d... \n',subIDs(sub),ses);
        catch
            fprintf(['No mat file found for ',sprintf('sub%02d session %d. skipping...\n',subIDs(sub),ses)]);
            do_preproc = 0;
        end

        if do_preproc

            % find which channel is for what
            phy = get_correct_channels(data);

            % Determine sample rate and tolerance
            sampInt  = phy.cardiac.interval;
            tol      = 0.005;          % 5 ms tolerance

            %data will be downsampled to 100 Hz
            sampIntNew = 1/100;

            %downsample data
            puls  = interp1(phy.cardiac.values,1:sampIntNew/sampInt:size(phy.cardiac.values,1))';
            resp  = interp1(phy.resp.values,1:sampIntNew/sampInt:size(phy.resp.values,1))';

            %trigger times are still in s!
            scan = phy.trigger.times;

            %find breaks in between runs
            [ind,runPulses] = get_runs(scan);

            %remove aborted runs that have fewer than number of scans
            %reported in get_study_specs for dicom import
            abortRuns = find_aborted_runs(ind,runPulses,incomplRuns);
            scan(abortRuns) = [];

            %% check for excess- or missing pulses based on trigger diffs
            pulse  = diff(scan);
            fprintf('\nChecking for excess or missing pulses ...\n')
            % Excess pulses
            while any(pulse<TR-tol)
                fprintf('Sub%02d: excess pulses found\n',subIDs(sub));
                wh          = find(pulse<TR-tol,1);
                scan(wh+1)  = [];
                fprintf('removing pulse %d\n',wh+1);
                pulse       = diff(scan);
            end

            % Missing pulses
            while any(discretize(pulse, [(2*TR)-tol, (2*TR)+tol]))
                fprintf('Sub%02d: missing pulses found\n',subIDs(sub));
                wh     = find((pulse >= (2*TR)-tol) & (pulse <= (2*TR)+tol),1);
                scan   = [scan(1:wh);scan(wh+1)+TR;scan(wh+1:end)];
                fprintf('adding pulse %d\n',wh+1);
                pulse  = diff(scan);
            end

            %% Check each run for number of correct pulses
            for run = 1:nRuns
                epi = spm_BIDS(BIDS,'data','sub',sprintf('%02d',subIDs(sub)),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',taskName,'type','bold');
                nScans = numel(spm_vol(epi{:,1}));

                [ind,runPulses] = get_runs(scan);
                indRun   = ind(run);
                nPulses  = runPulses(run);

                %remove dummy pulses
                if nDummies > 0
                    scan(indRun+1:indRun+nDummies) = [];
                    [ind,runPulses] = get_runs(scan);
                    indRun   = ind(run);
                    nPulses  = runPulses(run);
                end

                %identify runs that have wrong number of pulses
                %check for missing pulses in the end
                if nPulses < nScans
                    if nPulses == nScans - 1
                        fprintf('Sub%02d run%d: adding missing pulse in the end\n',subIDs(sub),run);
                    else
                        fprintf('Sub%02d run%d: more than one missing pulse in the end detected.\nAll missing pulses will be added but make sure to understand why there are several missing pulses!!\n',subIDs(sub),run);
                    end
                    scan  = [scan(indRun+1:indRun+nPulses); scan(indRun+nPulses)+TR; scan(indRun+nScans:end)];
                end

                %check for too many pulses in the end
                if nPulses > nScans
                    if nPulses == nScans + 1
                        fprintf('Sub%02d run%d: removing excess pulse in the end\n',subIDs(sub),run);
                    else
                        fprintf('Sub%02d run%d: more than one excess pulse in the end detected.\nThey will be deleted but make sure to understand why there are several excess pulses!\n',subIDs(sub),run);
                    end
                    scan(indRun+(nScans+1:nPulses)) = [];
                end

                %check if number of pulses and diff are (now) correct
                [ind, runPulses] = get_runs(scan);

                err1 = 0; err2 = 0;
                if runPulses(run) ~= nScans
                    fprintf('Number of pulses not correct in sub%02d run%d\n',subIDs(sub),run);
                    err1 = 1;
                end
                scanRun = scan(ind(run)+1:ind(run)+nScans);
                pulse   = diff(scanRun);
                if any(pulse < TR-tol) || any(pulse > TR+tol)
                    fprintf('Timing of pulses in sub%02d run%d is not correct\n',subIDs(sub),run);
                    err2 = 1;
                end
                if err1 == 0 && err2 == 0
                    fprintf('Sub%02d run%d ok\n',subIDs(sub),run);
                end

                %convert trigger times into 100 Hz and find run indeces
                scanRun = scanRun/sampIntNew;
                index = [ceil(scanRun(1)) floor(scanRun(end)+TR/sampIntNew)]; %1 more sample

                % convert trigger events into time series with same length as cardiac, resp etc
                % for bids requirements
                triggerTS = zeros((index(2)-index(1))+1,1);
                triggerTS(round(scanRun-scanRun(1))+1) = 1;

                %crop trigger, cardiac, and respiratory data
                crop.trigger      = triggerTS;
                crop.cardiac      = puls(index(1):index(2),:);
                crop.resp         = resp(index(1):index(2),:);

                %now save bids compatible tsv file with physio data
                saveName = fullfile(path.preprocDir,sprintf('sub-%02d',subIDs(sub)),sprintf('ses-%02d',ses),'func',sprintf('sub-%02d_ses-%02d_task-%s_run-%02d',subIDs(sub),ses,taskName,run));
                save_tsv_file(saveName,crop,sampIntNew);
                clear crop triggerTS
            end
        end
    end
end
