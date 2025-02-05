function sbp_convert_physio_smr(subIDs)

[path,vars]     = get_study_specs;
nSubs           = length(subIDs);

physioDir       = fullfile(path.baseDir,'physio');
smrRawDir       = fullfile(physioDir,'smr');
matDir          = fullfile(physioDir,'mat');

% Add CED
if isempty(which('CEDS64LoadLib'))
    error('Cannot find CED toolbox in your path. This toolbox is required for the next steps. \nFor more information see the readme.txt');
end
cedpath = fullfile(physioDir,'CEDMATLAB','CEDS64ML');
addpath(cedpath);
CEDS64LoadLib(cedpath);
loadlibrary ceds64int

for sub = 1:nSubs
    sub_name = sprintf('sub-%02d',subIDs(sub));
    fprintf('\n\n Running physio preprocessing for %s ...',sub_name)
    % Delete existing Mats
    fprintf('\ndeleting existing physio mats for %s ...\n',sub_name);
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

        fhand = CEDS64Open(ses_file_raw,1);
        [ NChannels ] = CEDS64MaxChan( fhand );
        [ dTBaseOut ] = CEDS64TimeBase( fhand ); % this is the time base, that is, the highest resolution time over all channels, on which i64Div is applied in each to get actual resolution

        fprintf('\n\nConversion SMR to MAT: file %s\n',ses_file_raw);

        sub_name = sprintf('sub%02d',subIDs(sub));

        m = matfile(mat_file_raw,'writable',true); % instantiate file
        for iChan = 1:NChannels
            [ iType ] = CEDS64ChanType( fhand, iChan ); % find out if the channel is a wave channel

            if ~iType
                continue;
            end

            fprintf('\nchannel %d %s',iChan,sub_name);

            tempStruct = struct;

            [ ~, sTitleOut ] = CEDS64ChanTitle( fhand, iChan );
            [ ~, sCommentOut ] = CEDS64ChanComment( fhand, iChan);

            i64MaxTime = CEDS64ChanMaxTime( fhand, iChan );
            if i64MaxTime==-1
                if iChan ~= 31                                                  % Chan 31 always gives a warning for some reason, but is being disregarded as it is not used with the settings used by our psychophysics settings.
                    warning('Channel %d is empty, skipping...',iChan)
                end
                continue;
            end

            [ i64Div ] = CEDS64ChanDiv( fhand, iChan );

            tempStruct.title = sTitleOut;
            tempStruct.comment = sCommentOut;

            fprintf('(%s)... ',sTitleOut);

            [ i64MaxTime ] = CEDS64ChanMaxTime( fhand, iChan );

            if iType==1 % data channel

                [ ~, dOffsetOut ]   = CEDS64ChanOffset( fhand, iChan );
                [ ~, sUnitsOut ]    = CEDS64ChanUnits( fhand, iChan );
                [ ~, fVals, i64Time ] = CEDS64ReadWaveF( fhand, iChan, ceil(i64MaxTime/i64Div), 0 );
                fVals = double(fVals);

                tempStruct.interval = dTBaseOut*i64Div; % ms resolution
                tempStruct.scale = 0; % unclear
                tempStruct.offset = dOffsetOut;
                tempStruct.units = sUnitsOut;
                tempStruct.start = i64Time;
                tempStruct.length = numel(fVals);
                tempStruct.values = fVals;

            elseif iType==3

                [ ~, i64Times ] = CEDS64ReadEvents( fhand, iChan, ceil(i64MaxTime/i64Div), 0);
                if isempty(i64Times)
                    if iChan ~= 31
                        warning('Channel %d of %s is empty, skipping...',iChan,sub_name)
                    end
                    continue;
                end
                [ dSeconds ] = CEDS64TicksToSecs( fhand, i64Times );

                tempStruct.resolution = dTBaseOut;
                tempStruct.length = numel(dSeconds);
                tempStruct.times = dSeconds;

            end

            m.(sprintf('%s_Ch%d',sub_name,iChan)) = tempStruct;

        end
        CEDS64Close(fhand); 
    end
end
unloadlibrary ceds64int;
end
