function convert_physio_smr(subIDs)

    %% Path specifications
    [path,vars]     = get_study_specs;
    nSubs           = length(subIDs);
    physioDir       = fullfile(path.baseDir,'physio');
    physio4cbdir    = fullfile(path.baseDir,'physio4cb');
    utilsdir        = fullfile(physio4cbdir,'utils');
    smrRawDir       = fullfile(physioDir,'smr');
    matDir          = fullfile(physioDir,'mat');
    
    addpath(physio4cbdir);
    addpath(utilsdir);
    addpath(fullfile(utilsdir,'smrReader'));

    %% Main Subject Loop
    for sub = 1:nSubs
        sub_name = sprintf('sub-%02d',subIDs(sub));
        fprintf('\n\nRunning physio preprocessing for %s ...',sub_name)
       
        % Delete existing Mats
        fprintf('\nDeleting existing physio mats for %s ...\n',sub_name);
        if exist(fullfile(matDir,sub_name),'dir')
            rmdir(fullfile(matDir,sub_name),'s')
        end
        mkdir(fullfile(matDir,sub_name));
    
        % Session loop
        for ses = 1:vars.nSess
            % SMR raw dir and filename for each session
            ses_file_raw = fullfile(smrRawDir,sprintf('sub-%02d_ses-%02d_physio.smr',subIDs(sub),ses));
    
            % Skip if we do not have physio smr file
            if isempty(regexp(ses_file_raw,'.*\.smr','MATCH'))
                fprintf('\nNo SMR file for sub %02d ses %02d. Skipping ...\n',subIDs(sub),ses);
                continue;
            end

            % Make filenames for raw and preprocessed mat files
            mat_file_raw   = fullfile(matDir,sprintf('sub-%02d',subIDs(sub)),sprintf('ses-%02d',ses),sprintf('sub-%02d_ses-%02d_physio.mat',subIDs(sub),ses));

            % Mat Raw and Preprocessed dirs for each session
            mkdir(fullfile(matDir,sprintf('sub-%02d',subIDs(sub)),sprintf('ses-%02d',ses)))

            import = ImportSMR(ses_file_raw);                               % This is what is imported by readSMR (https://github.com/tjrantal/Spike-smr-reader)
            NChannels  = import(end).hdr.channel;                           % Number of the final channel
            NChannels  = length(import);
    
            fprintf('\n\nConversion SMR to MAT: file %s\n',ses_file_raw);
    
            sub_name = sprintf('sub%02d',subIDs(sub));
    
            m = matfile(mat_file_raw,'writable',true); % instantiate file
            for iChan = 1:NChannels
    
                fprintf('\nchannel %d %s',iChan,sub_name);
    
                tempStruct = struct;
                tempStruct.title = import(iChan).hdr.title;
                tempStruct.comment = import(iChan).hdr.comment;
    
                fprintf('(%s)... ',tempStruct.title);
    
                if strcmp(import(iChan).hdr.channeltype, 'Continuous Waveform') % data channel
                    
                    tempStruct.interval = 1/import(iChan).hdr.adc.SampleInterval(1); % resolution
                    tempStruct.scale  = 0; %
                    tempStruct.offset = 0;
                    tempStruct.units  = import(iChan).hdr.adc.Units;
                    tempStruct.start  = import(iChan).imp.tim(1);
                    tempStruct.values = double(import(iChan).imp.adc)*import(iChan).hdr.adc.Scale;
                    tempStruct.length = numel(tempStruct.values);
                
                elseif strcmp(import(iChan).hdr.channeltype, 'Rising Edge')     % trigger channel
    
                    i64Times = import(iChan).imp.tim;
                    
                    if isempty(i64Times)
                            warning('Channel %d of %s is empty, skipping...',iChan,sub_name)
                        continue;
                    end
                    
                    tempStruct.resolution = import(iChan).hdr.tim.Scale*import(iChan).hdr.tim.Units;
                    dSeconds = double(i64Times) * tempStruct.resolution;
                    tempStruct.resolution = import(iChan).hdr.tim.Scale*import(iChan).hdr.tim.Units;
                    tempStruct.length = numel(dSeconds);
                    tempStruct.times = dSeconds;
                else
                    warning("Could not determine channel type for Channel %d, skipping ...",iChan)
                end
    
                    m.(sprintf('%s_Ch%d',sub_name,iChan)) = tempStruct;
    
            end
        end
    end
end
