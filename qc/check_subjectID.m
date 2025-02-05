function check_subjectID(subIDs)

% define variables and get paths/files
nSubs        = length(subIDs);
[path,vars]  = get_study_specs;
baseDir      = path.baseDir;
preprocPath  = path.preprocDir;

nRuns        = vars.nRuns;
nSess        = vars.nSess;
TR           = vars.sliceTiming.tr;

BIDS         = spm_BIDS(preprocPath);

table        = readtable(fullfile(preprocPath,'participants.tsv'),'FileType','text','Delimiter','\t');
subVar       = str2double(extract(table.participant_id,digitsPattern));
prismaIDs    = [subVar table.prisma_id]; %table.prisma_id

for sub = 1:nSubs
    prismaID   = prismaIDs(prismaIDs(:,1) == subIDs(sub),2);
    fprintf('Provided PRISMA ID for sub%02d is %d\n',subIDs(sub),prismaID);
    nScans = [];
    for sess = 1:nSess
        if nSess > 1
            fprintf('Analysing session %d now...\n',sess);
        end
        for run = 1:nRuns
            %% extract time stamp from epi files
            epi = spm_BIDS(BIDS,'data','sub',sprintf('%02d',subIDs(sub)),'ses',sprintf('%02d',sess),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            epi = epi(1); % if there are brain and spinal, just take brain
            V   = spm_vol(char(epi));
            nScans = [nScans size(V,1)];
            timeStampMRI{run} = datetime(V(1).descrip,'Format','dd/MM/yyyy HH:mm:ss');

            %% extract time stamp from behavioral files
            subDir = fullfile(baseDir,sprintf('sub-%02d/ses-%02d',subIDs(sub),sess),'beh');
            logfile = spm_select('FPList',subDir,sprintf('sub-%02d_ses-%02d_task-%s_run%02d.*.mat',subIDs(sub),sess,vars.task,run));
%             logfile = spm_select('FPList',subDir,sprintf('sub%03d_session%03d.*.mat',subIDs(sub),1));
            if size(logfile,1) > 1
                logfile = logfile(end,:);
            end
            file = dir(logfile);
            timeStampBeh{run} = datetime(file.date,'Format','dd/MM/yyyy HH:mm:ss');
        end
    end
    
    %% check both time stamps

    [y1,m1,d1] = ymd(timeStampMRI{1});
    [y2,m2,d2] = ymd(timeStampBeh{1});
    if ~isequal([y1,m1,d1],[y2,m2,d2])
        fprintf('Sub%02d: date for MRI images and logfile does not match\n',subIDs(sub));
        fprintf('One reason might be a mismatch between PRISMA and sub ID. Please verify!\n');
        rand1 = 1;
    end

    startMRI = timeofday(timeStampMRI{1}-seconds(nScans(1)*TR));
    stopMRI  = timeofday(timeStampMRI{end});
    midMRI   = (stopMRI-startMRI)/2;
    startBeh = timeofday(timeStampBeh{1}-seconds(nScans(1)*TR));
    stopBeh  = timeofday(timeStampBeh{end});
    midBeh   = (stopBeh-startBeh)/2;
    if ~isbetween(midBeh-minutes(8),midMRI+minutes(8))
        fprintf('Sub%02d: time stamps for MRI images and logfiles do not match\n',subIDs(sub));
        diffTime = minutes(timeofday(timeStampBeh{6})-timeofday(timeStampMRI{6}));
        fprintf('Time difference between last scan and behav logfile is %0.2f minutes\n',diffTime);
        fprintf('Please check wether this is plausible \nOtherwise this might indicate a mismatch between PRISMA and sub ID\n');
        rand2 = 1;
    end
    if ~exist('rand1','var') && ~exist('rand2','var')
        fprintf('Sub%02d: PRISMA ID seems to match for this sub\n',subIDs(sub));
    end

end
end