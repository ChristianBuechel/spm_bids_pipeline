function convert_and_clean_physio(all_sub_ids)

%% Path specifications
[path,vars,~,import]     = get_study_specs;
n_subs                   = length(all_sub_ids);
BIDS                     = spm_BIDS(path.preprocDir);
nSess                    = vars.nSess;
nRuns                    = vars.nRuns;
taskName                 = vars.task;
nDummies                 = import.dummies;
TR                       = vars.sliceTiming.tr;
incomplRuns              = str2double(extract(import.data(1).cond,digitsPattern)); % not always (1) ...rule for number of scans provided in get_study_specs for data import to find incomplete runs

%% Preprocessing: SMR to cleaned BIDS
%% 1. Retrieve Data in modular fashion from either run- or sessionwise smr- or matfiles.
% Main Subject Loop
for sub = 1:n_subs
    fprintf('\n\nRunning physio preprocessing for %s ...',sprintf('sub-%02d',all_sub_ids(sub)))

    % Session loop
    for ses = 1:nSess

        % First, delete all physio files
        fprintf('\nSession %d: deleting exisiting Physio files in func ....\n',ses);
        toggleWarning = warning('off','MATLAB:DELETE:FileNotFound');
        for k = 1:nRuns
            % tsv files
            tsvName = fullfile(path.preprocDir,sprintf('sub-%02d',all_sub_ids(sub)),sprintf('ses-%02d',ses),'func',sprintf('sub-%02d_ses-%02d_task-%s_run-%02d',all_sub_ids(sub),ses,taskName,k));
            delete(spm_file(tsvName,'suffix','_physio','ext','.tsv'),'FileType','text','Delimiter','\t'); % Delete the file
            % regressor txt files
            delete([path.preprocDir,sprintf('/sub-%02d/',all_sub_ids(sub)),sprintf('ses-%02d/',ses),'func/',sprintf('physio_sub-%02d_ses-%02d_task-%s_run-%02d_bold.txt',all_sub_ids(sub),ses,taskName,k)]);
            %physio info files
            delete([path.preprocDir,sprintf('/sub-%02d/',all_sub_ids(sub)),sprintf('ses-%02d/',ses),'func/',sprintf('info-physio_sub-%02d_ses-%02d_task-%s_run-%02d.mat',all_sub_ids(sub),ses,taskName,k)])
        end
        warning(toggleWarning);

        % By default, loop below runs nRun times. 
        iterations = nRuns;
  
        % But: if we have ses-wise files, it runs just once.
        if exist(fullfile(path.physioDir,'smr',sprintf('sub-%02d_ses-%02d_physio.smr',all_sub_ids(sub),ses)),'file')
            iterations = 1;
            convert_smr = 1;
            smr_file_raw = fullfile(path.physioDir,'smr',sprintf('sub-%02d_ses-%02d_physio.smr',all_sub_ids(sub),ses));
            fprintf('Processing SMR Session file ... \n');
            if exist(fullfile(path.physioDir,'smr',sprintf('sub-%02d_ses-%02d_physio.json',all_sub_ids(sub),ses)),'file')
                physioSidecar = spm_jsonread(fullfile(path.physioDir,'smr',sprintf('sub-%02d_ses-%02d_physio.json',all_sub_ids(sub),ses)));
            end
        elseif exist(fullfile(path.physioDir,'mat',sprintf('sub-%02d_ses-%02d_physio.mat',all_sub_ids(sub),ses)),'file')
            iterations = 1;
            convert_smr = 0;
            fprintf('Mat file found for sub %02d ses %02d. Continuing with session-long mat file.\n',all_sub_ids(sub),ses);
            data = load(fullfile(path.physioDir,'mat',sprintf('sub-%02d_ses-%02d_physio.mat',all_sub_ids(sub),ses)));
            if exist(fullfile(path.physioDir,'mat',sprintf('sub-%02d_ses-%02d_physio.json',all_sub_ids(sub),ses)),'file')
                physioSidecar = spm_jsonread(fullfile(path.physioDir,'mat',sprintf('sub-%02d_ses-%02d_physio.json',all_sub_ids(sub),ses)));
            end
        end

            %% Main Loop
           for iter = 1:iterations
            
               %% Import SMR file and convert to structure identical to mats
               % Check if we have runwise SMR files
               if exist(fullfile(path.physioDir,'smr',sprintf('sub-%02d_ses-%02d_run-%02d_physio.smr',all_sub_ids(sub),ses,iter)),'file')
                   smr_file_raw = fullfile(path.physioDir,'smr',sprintf('sub-%02d_ses-%02d_run-%02d_physio.smr',all_sub_ids(sub),ses,iter));
                   convert_smr = 1;
                   fprintf('Processing SMR Ses %02d Run %02d File\n',ses,iter);
               elseif iterations ~= 1
                   sprintf('Couldnt find smr file for sub %02d ses %02d run %02d. Looking for Mat Run Files...',all_sub_ids(sub),ses,iter);
                   convert_smr = 0;
               end
               
               if convert_smr
                smrimport = ImportSMR(smr_file_raw);                         % This is what is imported by readSMR (https://github.com/tjrantal/Spike-smr-reader)
                NChannels  = smrimport(end).hdr.channel;                     % Number of the final channel
                for iChan = 1:NChannels
                        fprintf('\nchannel %d %s',iChan,sprintf('sub-%02d',all_sub_ids(sub)));
                        
                        if ~isfield(smrimport(iChan),'hdr')
                                    warning('Channel %d of %s is empty, skipping...',iChan,all_sub_ids(sub))
                                continue;
                        end
                        
                        smr(iChan).title = smrimport(iChan).hdr.title;
                        smr(iChan).comment = smrimport(iChan).hdr.comment;

                        fprintf('(%s)... ',smr(iChan).title);

                        if strcmp(smrimport(iChan).hdr.channeltype, 'Continuous Waveform') % data channel

                            smr(iChan).interval = 1/smrimport(iChan).hdr.adc.SampleInterval(1); % resolution
                            smr(iChan).scale  = 0; %
                            smr(iChan).offset = 0;
                            smr(iChan).units  = smrimport(iChan).hdr.adc.Units;
                            smr(iChan).start  = smrimport(iChan).imp.tim(1);
                            smr(iChan).values = double(smrimport(iChan).imp.adc)*smrimport(iChan).hdr.adc.Scale;
                            smr(iChan).length = numel(smr(iChan).values);

                        elseif strcmp(smrimport(iChan).hdr.channeltype, 'Rising Edge')     % trigger channel

                            i64Times = smrimport(iChan).imp.tim;

                            if isempty(i64Times)
                                    warning('Channel %d of %s is empty, skipping...',iChan,all_sub_ids(sub))
                                continue;
                            end

                            smr(iChan).resolution = smrimport(iChan).hdr.tim.Scale*smrimport(iChan).hdr.tim.Units;
                            dSeconds = double(i64Times) * smr(iChan).resolution;
                            smr(iChan).resolution = smrimport(iChan).hdr.tim.Scale*smrimport(iChan).hdr.tim.Units;
                            smr(iChan).length = numel(dSeconds);
                            smr(iChan).times = dSeconds;
                        else
                            warning("Could not determine channel type for Channel %d, skipping ...",iChan)
                        end

                        % add to final data struct
                        eval(sprintf('data.sub%02d_Ch%d = smr(iChan);',all_sub_ids(sub),iChan))
                end
               end
            
            % Check if there are runwise mat files 
               if exist(fullfile(path.physioDir,'mat',sprintf('sub-%02d_ses-%02d_run-%02d_physio.mat',all_sub_ids(sub),ses,iter)),'file')
                   mat_file_run = fullfile(path.physioDir,'mat',sprintf('sub-%02d_ses-%02d_run-%02d_physio.mat',all_sub_ids(sub),ses,iter));
                   fprintf('Processing Mat Ses %02d Run %02d File\n',ses,iter);
                   data = load(mat_file_run);
               elseif iterations ~=1 && ~convert_smr
                   warning('Couldnt find mat or smr- runwise file for sub %02d ses %02d run %02d. Skipping...',all_sub_ids(sub),ses,iter);
               end
            
            %% Continue. Check for missing- or ghost pulses (and divide into runs if ses-files)
            if ~exist('data','var')
                  warning('No file(s) found for sub %02d ses %02d. Skipping ...',all_sub_ids(sub),ses);
            else
                % get correct channels
                for i = 1:numel(fieldnames(data))
                    if eval(sprintf('strcmp(data.sub%02d_Ch%d.title,"scanner")',all_sub_ids(sub),i))
                        fprintf('\nscan trigger channel is channel %d\n',i);
                        eval(sprintf('phy.trigger = data.sub%02d_Ch%d;',all_sub_ids(sub),i));
                    elseif eval(sprintf('strcmp(data.sub%02d_Ch%d.title,"PULS")',all_sub_ids(sub),i))
                        fprintf('\npulse channel is channel %d\n',i);
                        eval(sprintf('phy.cardiac = data.sub%02d_Ch%d;',all_sub_ids(sub),i));
                    elseif eval(sprintf('strcmp(data.sub%02d_Ch%d.title,"Resp")',all_sub_ids(sub),i'))
                        fprintf('\nresp channel is channel %d\n',i);
                        eval(sprintf('phy.resp = data.sub%02d_Ch%d;',all_sub_ids(sub),i));
                    end
                end
    
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
                %find runs that have fewer than expected scans
                abort = find(runPulses<incomplRuns);
                remove = [];
    
                if ~isempty(abort)
                    remove = [];
                    disp('Deleting aborted runs')
                    for ar = 1:length(abort)
                        remove = [remove ind(abort(ar))+1:ind(abort(ar))+runPulses(abort(ar))];
                    end
                end        
                scan(remove) = [];
    
                %% Now remove runs that were completed but user would like to exclude
                if exist('physioSidecar','var')
                    [ind,runPulses] = get_runs(scan);
                    exclude = find(isnan(physioSidecar.includedRuns));
                    remove = [];
                    disp('Deleting excluded runs')
                    for ex = 1:length(exclude)
                        remove = [remove ind(exclude(ex))+1:ind(exclude(ex))+runPulses(exclude(ex))];
                    end
                    scan(remove) = [];
                    [ind,runPulses] = get_runs(scan);
                    if length(ind)-1 ~= length(runPulses)
                        error('incorrect number of runs, even after correction w/ json sidecar. Please check your data!')
                    end
                else
                    [ind,runPulses] = get_runs(scan);
                    if length(ind)-1 ~= length(runPulses)
                        error(sprintf('incorrect number of runs! Please provide json sidecar to sub %02d ses %02d in your physioDir to determine runs to be excluded.', all_sub_ids(sub),ses))
                    end
                end



                %% check for excess- or missing pulses based on trigger diffs
                pulse  = diff(scan);
                fprintf('\nChecking for excess or missing pulses ...\n')
                % Excess pulses
                while any(pulse<TR-tol)
                    fprintf('Sub%02d: excess pulses found\n',all_sub_ids(sub));
                    wh          = find(pulse<TR-tol,1);
                    scan(wh+1)  = [];
                    fprintf('removing pulse %d\n',wh+1);
                    pulse       = diff(scan);
                end
    
                % Missing pulses
                while any(discretize(pulse, [(2*TR)-tol, (2*TR)+tol]))
                    fprintf('Sub%02d: missing pulses found\n',all_sub_ids(sub));
                    wh     = find((pulse >= (2*TR)-tol) & (pulse <= (2*TR)+tol),1);
                    scan   = [scan(1:wh);scan(wh+1)+TR;scan(wh+1:end)];
                    fprintf('adding pulse %d\n',wh+1);
                    pulse  = diff(scan);
                end
    
                %% Check each run for number of correct pulses
                switch iterations                           % if we have runwise SMR/Mat files, we do the loop below just once
                    case nRuns
                        runloop = 1;
                    case 1
                        runloop = nRuns;                    % if we have sess-wise SMR/Mat files, we do the loop below nRun-times
                end
    
                    for run = 1:runloop
                        %Catch correct epis
                        switch runloop
                            case 1  
                                epi = spm_BIDS(BIDS,'data','sub',sprintf('%02d',all_sub_ids(sub)),'ses',sprintf('%02d',ses),'run',sprintf('%02d',iter),'task',taskName,'type','bold');
                            case nRuns
                                epi = spm_BIDS(BIDS,'data','sub',sprintf('%02d',all_sub_ids(sub)),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',taskName,'type','bold');
                        end
        
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
                                fprintf('Sub%02d run%d: adding missing pulse in the end\n',all_sub_ids(sub),run);
                            else
                                fprintf('Sub%02d run%d: more than one missing pulse in the end detected.\nAll missing pulses will be added but make sure to understand why there are several missing pulses!!\n',all_sub_ids(sub),run);
                            end
                            scan  = [scan(1:indRun);scan(indRun+1:indRun+nPulses); scan(indRun+nPulses)+TR; scan(indRun+nScans:end)];
                        end
        
                        %check for too many pulses in the end
                        if nPulses > nScans
                            if nPulses == nScans + 1
                                fprintf('Sub%02d run%d: removing excess pulse in the end\n',all_sub_ids(sub),run);
                            else
                                fprintf('Sub%02d run%d: more than one excess pulse in the end detected.\nThey will be deleted but make sure to understand why there are several excess pulses!\n',all_sub_ids(sub),run);
                            end
                            scan(indRun+(nScans+1:nPulses)) = [];
                        end
        
                        %check if number of pulses and diff are (now) correct
                        [ind, runPulses] = get_runs(scan);
        
                        err1 = 0; err2 = 0;
                        if runPulses(run) ~= nScans
                            fprintf('Number of pulses not correct in sub%02d run%d\n',all_sub_ids(sub),run);
                            err1 = 1;
                        end
                        scanRun = scan(ind(run)+1:ind(run)+nScans);
                        pulse   = diff(scanRun);
                        if any(pulse < TR-tol) || any(pulse > TR+tol)
                            fprintf('Timing of pulses in sub%02d run%d is not correct\n',all_sub_ids(sub),run);
                            err2 = 1;
                        end
                        if err1 == 0 && err2 == 0
                            fprintf('Sub%02d run%d ok\n',all_sub_ids(sub),run);
                        end
        
                        %convert trigger times into 100 Hz and find run indeces
                        scanRun = scanRun/sampIntNew;
                        if run == runloop && scanRun(end)+TR/sampIntNew > length(puls)
                            warning('Data Series is shorter than MR collection + 1TR. Most likely, physio data collection was stopped before MR data collection was. Cropping trigger series to match length of physio data collection. Tapas will interpolate missing data, but please check your output carefully!');
                            while floor(scanRun(end)+TR/sampIntNew) > length(puls)
                                scanRun(end) = [];
                            end
                                index = [ceil(scanRun(1)) length(puls)];
                        else
                            index = [ceil(scanRun(1)) floor(scanRun(end)+TR/sampIntNew)];%+TR/sampIntNew)]; %1 more sample
                        end
        
                        % convert trigger events into time series with same length as cardiac, resp etc
                        % for bids requirements
                        triggerTS = zeros((index(2)-index(1))+1,1);
                        triggerTS(round(scanRun-scanRun(1))+1) = 1;
        
                        %crop trigger, cardiac, and respiratory data
                        crop.trigger      = triggerTS;
                        crop.cardiac      = puls(index(1):index(2),:);
                        crop.resp         = resp(index(1):index(2),:);
        
                        % final matrix
                        matrix = [crop.cardiac,crop.resp,crop.trigger];
        
                        %% 3. Convert to BIDS compatible tsv file (incl. json sidecar) for physio data
                        saveName = fullfile(path.preprocDir,sprintf('sub-%02d',all_sub_ids(sub)),sprintf('ses-%02d',ses),'func',sprintf('sub-%02d_ses-%02d_task-%s_run-%02d',all_sub_ids(sub),ses,taskName,run));
        
                        % write json header
                        bold_json = struct('SamplingFrequency',1/sampIntNew,'StartTime',0,'Columns',["cardiac","respiratory","trigger"]);
                        spm_jsonwrite(spm_file(saveName,'suffix','_physio','ext','.json'),bold_json, struct('indent','  '));
        
                        % write BIDS. keep only gz tsv
                        writematrix(matrix,spm_file(saveName,'suffix','_physio','ext','.tsv'),'FileType','text','Delimiter','\t');
                        gzip(spm_file(saveName,'suffix','_physio','ext','.tsv'));
                        delete(spm_file(saveName,'suffix','_physio','ext','.tsv'));
                    end
                
            end
            clear data
           end
    end
end

end


