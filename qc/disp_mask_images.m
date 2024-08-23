function disp_mask_images(subIDs,fLevelDir)


% resolve paths
[path,vars,~,~,qc]  = get_study_specs;
fLevelPath          = path.firstlevelDir;
nSubs               = length(subIDs);

maxImg              = qc.maxDisImg;

tempFile            = spm_select('FPList',path.templateDir,vars.templateStripID);
allFiles            = [];
    
for sub = 1:nSubs
    workDir    = fullfile(fLevelPath,sprintf('sub-%02d',subIDs(sub)),fLevelDir);
    subFile    = fullfile(workDir,spm_file('mask.nii','prefix','w'));
    allFiles   = [allFiles;subFile]; 
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






