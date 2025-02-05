function disp_coreg_images(subIDs)


% resolve paths
[path,vars,~,~,qc]  = get_study_specs;
preprocPath         = path.preprocDir;
nSubs               = length(subIDs);

BIDS                = spm_BIDS(preprocPath);

for sub = 1:nSubs
    anatFile  = spm_BIDS(BIDS,'data','sub',sprintf('%02d',subIDs(sub)),'type','T1w');
    anatFile  = char(anatFile(1)); %in case there are more sessions with T1w
    anatDir   = spm_file(anatFile,'path');
    skullFile = fullfile(anatDir,vars.skullStripID);
    
    epi         = spm_BIDS(BIDS,'data','sub',sprintf('%02d',subIDs(sub)),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    epi         = epi(1); % if there are brain and spinal, just take brain
    meanFile    = char(spm_file(epi,'prefix','cmeana'));
    
    dispImages = char(skullFile,meanFile);
    spm_check_registration(dispImages);
    if qc.contour == 1
        spm_orthviews('contour','display',1,2:size(dispImages,1))
    end
    fprintf('Checking Subject %d \nPress enter in command window to continue\n',subIDs(sub));
    input('');
end

end






