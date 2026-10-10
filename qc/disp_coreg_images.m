function disp_coreg_images(subIDs)


% resolve paths
[path,vars,~,~,qc]  = get_study_specs;
preprocPath         = path.preprocDir;
nSubs               = length(subIDs);

BIDS                = spm_BIDS(preprocPath);

for sub = 1:nSubs
    anatFile  = char(get_anat(BIDS,subIDs(sub)));
    anatDir   = spm_file(anatFile,'path');
    skullFile = fullfile(anatDir,vars.skullStripID);
    
    epi         = get_ref_epi(BIDS,subIDs(sub));
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






