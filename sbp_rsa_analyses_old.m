function sbp_rsa_analyses(all_sub_ids)

addpath(genpath('c:\Users\buechel\Documents\MATLAB\rsatoolbox_matlab\'));
% get all infos
[path,vars,~,~,rsav] = get_study_specs;
BIDS                 = spm_BIDS(path.preprocDir);

% some variables
xA = spm_atlas('load',rsav.atlas);

% define which ROIs
if isempty(rsav.roi_ind)
    roi_ind = cat(2,xA.labels.index); %all ROIs in atlas
    roi_names = strvcat(xA.labels.name);
else
    roi_ind   = rsav.roi_ind;
    roi_names = [];%??? --> ToDo
end


n_subs        = numel(all_sub_ids);



% load and rearrange
good_rois = zeros(n_subs,numel(roi_ind));
for sub=1:n_subs
    sub_id     = all_sub_ids(sub);
    ana_dir    = fullfile(path.firstlevelDir,sprintf('sub-%2.2d',sub_id),rsav.ana_name);
    if rsav.lss
        a          = load(fullfile(ana_dir,sprintf('RDMi_%s.mat',rsav.atlas)));
    else
        a          = load(fullfile(ana_dir,sprintf('lsa_RDMi_%s.mat',rsav.atlas)));
    end
    for gg=1:numel(roi_ind)
        roi(gg).RDMs(sub) = a.RDMi.dist(gg);
        roi(gg).RDMs(sub).RDM = roi(gg).RDMs(sub).RDM; % use raw
        % roi(gg).RDMs(sub).RDM = roi(gg).RDMs(sub).RDM_raw; % use raw
    
        if sub == 1 
            roi(gg).avgRDM.RDM     = a.RDMi.dist(gg).RDM./n_subs;
            roi(gg).avgRDM.RDM_raw = a.RDMi.dist(gg).RDM_raw./n_subs;
        else
            roi(gg).avgRDM.RDM     = roi(gg).avgRDM.RDM + a.RDMi.dist(gg).RDM./n_subs;
            roi(gg).avgRDM.RDM_raw = roi(gg).avgRDM.RDM_raw + a.RDMi.dist(gg).RDM_raw./n_subs;
        end
    roi(gg).avgRDM.name = a.RDMi.dist(gg).name;
    end
    good_rois(sub,:) = a.RDMi.good_rois;
end

mask = all(good_rois);
% now prune by good_rois
pruned_roi_ind = roi_ind(find(mask));
for gg=1:numel(pruned_roi_ind)
    
    % we do the sign rank test ourselves here
    % first get the similarity between data RDM and candiate RDMs
    for c_i = 1:numel(rsav.models)
        for s_i = 1:nSubjects
            cand2refSims(subI,candI)corr(vectorizeRDMs(meanCandRDMs(:,:,candI))',vectorizeRDMs(refRDM_stack(:,:,subI))','type',userOptions.RDMcorrelationTyp
        end
    end
    
    for candI = 1:numel(candRDMs)
        [ps(candI)] = signrank_onesided(cand2refSims(:,candI));
    end
    
        
     
     
    userOptions.saveFiguresPS = 0;
    userOptions.saveFiguresFig = 0;
    userOptions.saveFiguresPDF = 0;
    userOptions.saveFiguresPS  = 0;

    userOptions.projectName   = rsav.ana_name;
    userOptions.analysisName  = sprintf('%s',roi_names(pruned_roi_ind(gg),:));
    userOptions.rootPath      = fullfile(path.secondlevelDir,rsav.ana_name);
    localOptions.figureNumber = 101;
    %rsa.figureRDMs(roi(pruned_roi_ind(gg)).RDMs,userOptions,localOptions)
    localOptions.figureNumber = 102;
    %rsa.figureRDMs(roi(pruned_roi_ind(gg)).avgRDM.RDM,userOptions,localOptions);    
    
    
    
    % start with MDS
    userOptions.conditionLabels = rsav.stims;
    userOptions.rubberbands     = false;
    localOptions.figureNumber = 103;
    clf(localOptions.figureNumber);
    userOptions.conditionColours = rsav.condition_colors';
    rsa.MDSConditions(roi(pruned_roi_ind(gg)).avgRDM,userOptions, localOptions);

    %userOptions.RDMcorrelationType = 'Kendall_taua';
    %userOptions.RDMcorrelationType = 'Pearson';
    userOptions.RDMcorrelationType = 'Spearman';
    
    userOptions.RDMrelatednessTest = 'subjectRFXsignedRank';
    %userOptions.RDMrelatednessTest = 'randomisation';
    userOptions.RDMrelatednessThreshold = 0.05;
    userOptions.figureIndex = [104 105];
    userOptions.RDMrelatednessMultipleTesting = 'FDR';
    
    userOptions.candRDMdifferencesTest = 'subjectRFXsignedRank';
    userOptions.candRDMdifferencesThreshold = 0.05;
    userOptions.candRDMdifferencesMultipleTesting = 'FDR';
    userOptions.barsOrderedByRDMCorr = 0;
    stats_p_r = rsa.compareRefRDM2candRDMs(roi(pruned_roi_ind(gg)).RDMs, rsav.models, userOptions);
end
