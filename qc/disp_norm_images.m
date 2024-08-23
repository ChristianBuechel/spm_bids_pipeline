function disp_norm_images(subIDs,whichFile)


% resolve paths
[path,vars,~,~,qc]  = get_study_specs;
preprocPath         = path.preprocDir;
nSubs               = length(subIDs);

BIDS                = spm_BIDS(preprocPath);
maxImg              = qc.maxDisImg;

tempFile            = spm_select('FPList',path.templateDir,vars.templateStripID);
allFiles            = [];

for sub = 1:nSubs
    if whichFile == 1
         anatFile   = spm_BIDS(BIDS,'data','sub',sprintf('%02d',subIDs(sub)),'type','T1w');
         anatFile   = anatFile(1); %in case there are more sessions with T1w
         anatDir    = spm_file(anatFile,'path');
         workFile   = char(fullfile(anatDir,spm_file(vars.skullStripID,'prefix','w')));
    elseif whichFile == 2
        epi         = spm_BIDS(BIDS,'data','sub',sprintf('%02d',subIDs(sub)),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
        epi         = epi(1); % if there are brain and spinal, just take brain
        workFile    = char(spm_file(epi,'prefix','wmeana'));
    end
    allFiles   = [allFiles;workFile]; 
end

for i = 1:ceil(nSubs/maxImg)  
    if nSubs < maxImg
        spm_check_registration(char(tempFile,allFiles(maxImg*(i-1)+1:end,:)));
        contInd = nSubs + 1;
    elseif i == ceil(nSubs/maxImg)
        spm_check_registration(char(tempFile,allFiles(maxImg*(i-1)+1:end,:)));
        contInd = ceil(nSubs/maxImg) + 1;
    else
        spm_check_registration(char(tempFile,allFiles(maxImg*(i-1)+1:maxImg*(i-1)+maxImg,:)));
        contInd = maxImg + 1;
    end
    if qc.contour == 1
        spm_orthviews('contour','display',1,2:contInd)
    end
    fprintf('Press enter in command window to continue\n');
    input('');
end

end






