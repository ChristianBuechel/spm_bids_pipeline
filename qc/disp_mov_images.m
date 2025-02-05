function disp_mov_images(subIDs)

% resolve paths
[path,vars,~,~,qc]  = get_study_specs;
preprocPath         = path.preprocDir;
nSubs               = length(subIDs);
nRuns               = vars.nRuns;
nSess               = vars.nSess;

BIDS                = spm_BIDS(preprocPath);

for sub = 1:nSubs 
    epi         = spm_BIDS(BIDS,'data','sub',sprintf('%02d',subIDs(sub)),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    epi         = epi(1); % if there are brain and spinal, just take brain
    meanFile    = char(spm_file(epi,'prefix','meana'));
    allFiles = [];
    for sess = 1:nSess
        for run = 1:nRuns
            epi                =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',subIDs(sub)),'ses',sprintf('%02d',sess),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            epi                =  epi(1); % if there are brain and spinal, just take brain
            workFile = spm_select('ExtFPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^ra'),1);
            allFiles = [allFiles;workFile];
        end
    end
    
    dispImages = char(meanFile,allFiles);
    spm_check_registration(dispImages);
    if qc.contour == 1
        spm_orthviews('contour','display',1,2:size(dispImages,1))
    end
    
%     global st
%         % your desired params
%     nLines = 1;
%     styleLines = 'r-' % Line properties z.B. via line
%     % contour to modify
%     ctm = 1;
%     hM = findobj(st.vols{ctm}.ax{1}.cm,'Label','Contour');
%     UD = get(hM,'UserData');
%     UD.nblines = nLines;
%     UD.style = styleLines;
%     set(hM,'UserData',UD);
%     spm_ov_contour('display',ctm,Inf)
    
    fprintf('Checking Subject %d \nPress enter in command window to continue\n',subIDs(sub));
    input('');
end

end






