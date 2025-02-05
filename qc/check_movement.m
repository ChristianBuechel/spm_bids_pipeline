function check_movement(subIDs,movReg,eventName)

[path,vars,~,~,qc]  = get_study_specs;
preprocPath         = path.preprocDir;
nRuns               = vars.nRuns;
nSubs               = length(subIDs);
nSess               = vars.nSess;

BIDS                = spm_BIDS(preprocPath);

for sub = 1:nSubs
    mParams = [];
    pOnsets = [];
    nScans  = [];
    for sess = 1:nSess
        fprintf('Checking sub%02d',subIDs(sub));
        if nSess > 1
            fprintf(' - session%d',sess);
        end
        for run = 1:nRuns
            epi = spm_BIDS(BIDS,'data','sub',sprintf('%02d',subIDs(sub)),'ses',sprintf('%02d',sess),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            epi = epi(1); % if there are brain and spinal, just take brain
            rpFile  = char(spm_file(epi,'prefix','rp_a','ext','txt'));
            matFile = char(spm_file(epi,'prefix','a','ext','mat'));
            runData = load(rpFile);
            load(matFile,'mat');
            maxMov = max(abs(runData(:,1:3)))';
            if any(maxMov > qc.threshMov)
                maxM = max(maxMov(maxMov > qc.threshMov));
                fprintf('Sub%02d Run%d: maximum movement of %.2f mm detected, please check participant\n',subIDs(sub),run,maxM);
            end
            matFirst(:,:,run)   = mat(:,:,2);
            matLast(:,:,run)    = mat(:,:,end);
            withinSessMov       = matLast(:,:,run)/matFirst(:,:,run);
            wsm(:,run)          = abs(withinSessMov(1:3,4));
            if any(wsm > qc.threshMov)
                fprintf('Sub%02d Run%d: movement between first and last image of %.2f mm detected, please check participant\n',subIDs(sub),run,max(wsm));
            end
            if ~isempty(eventName)
                ons = spm_BIDS(BIDS,'data','sub',sprintf('%02d',subIDs(sub)),'ses',sprintf('%02d',sess),'run',sprintf('%02d',run),'task',vars.task,'type','events');
                onsets = readtable(char(ons),'FileType','text','Delimiter','\t');
                pOnsets = [pOnsets; (onsets{(contains(onsets{:,3},eventName)),1}/vars.sliceTiming.tr)+size(mParams,1)];
            end
            mParams = [mParams; runData];
            nScans  = [nScans size(mParams,1)];
        end
        
        %now plot movement and between run movement
        hFig = figure;
        set(hFig,'units','normalized','pos',[0.3,0.25,0.45,0.3]);
        plot(mParams(:,1:3));hold on;
        if ~isempty(eventName)
            xline(pOnsets); %,'HandleVisibility','off'
            legend('x','y','z',eventName);
        else
            legend('x','y','z');
        end
        if nSess > 1
            title(sprintf('Sub%02d - Sess%d',subIDs(sub),sess));
        else
            title(sprintf('Sub%02d',subIDs(sub)));
        end
        xlabel('Number of volumes');
        ylabel('Movement (mm)');
        xlim([1 size(mParams,1)]);
        ylim([min(mParams(:))+min(mParams(:))*0.2 max(mParams(:))+max(mParams(:))*0.2]);
        xticks([1 nScans]);
        set(gcf,'color','w');
        if nRuns > 1
            for run = 1:nRuns-1
                betSessMov(:,:,run) = matLast(:,:,run)/matFirst(:,:,run+1);
                mov2First(:,:,run) = matLast(:,:,run)/matFirst(:,:,1);
                bsm(:,run) = abs(betSessMov(1:3,4,run));
            end
            if any(bsm(:) > qc.threshMov)
                fprintf('Sub%02d: movement between runs of %.2f mm detected, please check participant\n',subIDs(sub),max(bsm(:)));
            end
            iFig = figure;
            set(iFig,'units','normalized','pos',[0.4,0.55,0.25,0.25]);
            if nSess > 1
                sgtitle(sprintf('Sub%02d - Sess%d',subIDs(sub),sess));
            else
                sgtitle(sprintf('Sub%02d',subIDs(sub)));
            end
            subplot(2,1,1)
            plot([zeros(3,1) bsm]','o-');xlim([0.5 nRuns+0.5]);ylim([-3 3]);title('Movement between runs');
            ylabel('Movement (mm)');
            subplot(2,1,2)
            plot(wsm','o-');xlim([0.5 nRuns+0.5]);ylim([-3 3]);title('Movement within runs (last-first volume)');
            xlabel('Run Number');
            ylabel('Movement (mm)');
            set(gcf,'color','w');       
        end
        fprintf('\nPress enter in command window to continue\n');
        input('');
        close(findobj('type','figure','name','hFig'));
        close(findobj('type','figure','name','iFig'));
    end
end

if movReg == 1
    % find spike movement and create noise regressors for first level
    for sub = 1:nSubs
        for sess = 1:nSess
            fprintf('Checking sub%02d',subIDs(sub));
            if nSess > 1
                fprintf(' - session%d',sess);
            end
            fprintf('\n');
            for run = 1:nRuns
                epi = spm_BIDS(BIDS,'data','sub',sprintf('%02d',subIDs(sub)),'ses',sprintf('%02d',sess),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
                epi = epi(1); % if there are brain and spinal, just take brain
                rpFile  = char(spm_file(epi,'prefix','rp_a','ext','txt'));
                mParams = load(rpFile);
                diffParams = abs(diff(mParams(:,1:3)));
                [nSpikes,c] = find(diffParams>qc.threshSpike);
                nSpikes = unique(nSpikes);
                if numel(nSpikes) ~= 0
                    if numel(nSpikes) > floor(size(mParams,1)*qc.percSpike)
                        fprintf('Run%d: spike movement in more than %.1f%% of volumes detected, please check participant\n',run,qc.percSpike*100);
                    end
                    fprintf('Run%d: %d spikes detected\n',run,numel(nSpikes));
                    mov_reg = zeros(size(mParams,1),numel(nSpikes));
                    for m = 1:length(nSpikes)
                        mov_reg(nSpikes(m)+1,m) = 1;
                    end
                    fprintf('Saving noise regressors for respective images\n');
                    outFilename = sprintf(qc.movFilename,subIDs(sub),sess,vars.task,run);
                    save(fullfile(char(spm_file(epi,'path')),outFilename),'mov_reg');
                end
            end
        end
    end
end

