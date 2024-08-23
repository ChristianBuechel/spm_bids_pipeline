function check_flevel_tmap(subIDs,fLevelDir,tMapNumber)

[path,vars,~,~,qc]  = get_study_specs;
preprocPath         = path.preprocDir;
fLevelPath          = path.firstlevelDir;
nSubs               = length(subIDs);

BIDS                = spm_BIDS(preprocPath);

tempFile            = spm_select('FPList',path.templateDir,vars.templateStripID);
rawtMapID           = sprintf('spmT_0%03d.nii',tMapNumber);

for sub = 1:nSubs
    workDir   = fullfile(fLevelPath,sprintf('sub-%02d',subIDs(sub)),fLevelDir);
    tMapID    = spm_file(rawtMapID,'prefix',['s' num2str(qc.tMapSmoothK) 'w'],'suffix',['_t' num2str(qc.tThresh*10)]);
    tMapFile  = fullfile(workDir,tMapID);
    
    if ~exist(tMapFile,'file')
        epi        = spm_BIDS(BIDS,'data','sub',sprintf('%02d',subIDs(sub)),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
        epi        = epi(1); % if there are brain and spinal, just take brain
        funcDir    = spm_file(epi,'path');
        flowFieldFile =  fullfile(funcDir,'y_epi_2_template.nii');
        
        rawtMapFile    = fullfile(workDir,rawtMapID);
        warptMapFile   = fullfile(workDir,spm_file(rawtMapID,'prefix','w'));
        smoothtMapFile = fullfile(workDir,spm_file(rawtMapID,'prefix',['s' num2str(qc.tMapSmoothK) 'w']));
  
        matlabbatch{1}.spm.util.defs.comp{1}.def = flowFieldFile;
        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = cellstr(rawtMapFile);
        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
        matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
        matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';
        
        matlabbatch{2}.spm.spatial.smooth.data = cellstr(warptMapFile);
        matlabbatch{2}.spm.spatial.smooth.fwhm = repmat(qc.tMapSmoothK,1,3);
        matlabbatch{2}.spm.spatial.smooth.prefix = ['s' num2str(qc.tMapSmoothK)];

        matlabbatch{3}.spm.util.imcalc.input = {smoothtMapFile};
        matlabbatch{3}.spm.util.imcalc.output = tMapID;
        matlabbatch{3}.spm.util.imcalc.outdir = {workDir};
        matlabbatch{3}.spm.util.imcalc.expression = ['i1.*(i1>' num2str(qc.tThresh) ')'];
        matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        
        spm_jobman('run',matlabbatch);
        clear matlabbatch
    end
    
    spm_check_registration(tempFile);
    
    tempFileInd = 1;
    overlayInd = 1;
    spm_orthviews('addcolouredimage',tempFileInd,tMapFile,qc.overlayColor)
    spm_orthviews('redraw',tempFileInd);
    spm_orthviews('addcolourbar',tempFileInd,overlayInd);
    spm_orthviews('redraw',tempFileInd);
    
    fprintf('Checking sub%02d \nPress enter in command window to continue\n',subIDs(sub));
    input('');
    
end
