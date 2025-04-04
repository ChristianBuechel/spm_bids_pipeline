function sbp_rsa_analyses(all_sub_ids)

% analyse atlas based RSA data and write statistical maps for the
% comparison with each RDM as specified in rsa.models
% uses non-parametric signed rank tests (like the RSA toolbox)
% the images are stored as -logP (negative log of the p-value)

% get all infos
[path,vars,~,~,rsa] = get_study_specs;
BIDS                 = spm_BIDS(path.preprocDir);

% some variables
xA = spm_atlas('load',rsa.atlas);

% define which ROIs
if isempty(rsa.roi_ind)
    roi_ind = cat(2,xA.labels.index); %all ROIs in atlas
    roi_names = strvcat(xA.labels.name);
else
    roi_ind   = rsa.roi_ind;
    roi_names = [];%??? --> ToDo
end

n_subs        = numel(all_sub_ids);
% load and rearrange
good_rois = zeros(n_subs,numel(roi_ind));
for sub=1:n_subs
    sub_id     = all_sub_ids(sub);
    ana_dir    = fullfile(path.firstlevelDir,sprintf('sub-%2.2d',sub_id),rsa.ana_name);
    if rsa.lss
        a          = load(fullfile(ana_dir,sprintf('lss_RDMi_%s.mat',rsa.atlas)));
    else
        a          = load(fullfile(ana_dir,sprintf('lsa_RDMi_%s.mat',rsa.atlas)));
    end
    for gg=1:numel(roi_ind)
        roi(gg).RDMs(sub) = a.RDMi.dist(gg);
        roi(gg).RDMs(sub).RDM = roi(gg).RDMs(sub).RDM; % use raw
    
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
    n_valid_subs = numel(roi(pruned_roi_ind(gg)).RDMs);
    
    % get Spearman between data and candiadte RDMs
    for c_i = 1:numel(rsa.models)
        if size(rsa.models{c_i}.RDM,1) > 1
            c_model_vec = squareform(rsa.models{c_i}.RDM);
        else
            c_model_vec = rsa.models{c_i}.RDM;
        end
        for s_i = 1:n_valid_subs
            if rsa.mnn == 1
                d_vec = roi(pruned_roi_ind(gg)).RDMs(s_i).RDM;
            else
                d_vec = roi(pruned_roi_ind(gg)).RDMs(s_i).RDM_raw;
            end
            c2r_corr(s_i,c_i) = corr(c_model_vec',d_vec','type','spearman','rows','pairwise');
        end
    end
    
    % we do the sign rank test ourselves here    
    for c_i = 1:numel(rsa.models)
        p = signrank(c2r_corr(:,c_i),[],'method','exact');
        if median(c2r_corr(:,c_i)) > 0
            p_values(gg,c_i) = p/2;
        else
            p_values(gg,c_i) = 1-p/2;
        end
    end    
end    

% write atlas images
for i=1:numel(roi_ind)
    xY(i).name = xA.labels(i).name;
    xY(i).ind  = xA.labels(i).index;
    xY(i).spec = xA.info.files.images;
    xY(i).def  = 'mask';
    xY(i).xyz  = Inf;
    xY(i).Ic   = NaN;
    xY(i).sk   = 0;  
end

V_atlas = spm_vol(xA.info.files.images);

new_image_V       = V_atlas;
new_image_V.dt(1) = spm_type('float32');
if rsa.mnn == 1
    descr = 'mnn_';
else
    descr = 'raw_';
end
if rsa.lss == 1
    descr = [descr 'lss_'];
else
    descr = [descr 'lsa_'];
end
for g = 1:numel(rsa.models)    
    new_image_V.fname = fullfile(path.secondlevelDir,rsa.ana_name,sprintf('%s%s_%s',descr, rsa.models{g}.name,'-logP.nii'));    
    new_data = nan(new_image_V.dim);
    for roi = 1:numel(pruned_roi_ind)
        [~, XYZmm, ~] = spm_ROI_ind(xY(pruned_roi_ind(roi)), V_atlas);
        vox = new_image_V.mat\[XYZmm;ones(1,size(XYZmm,2))];
        new_data(sub2ind(size(new_data),vox(1,:),vox(2,:),vox(3,:))) = -log(p_values(roi,g));
    end
    spm_write_vol(new_image_V,new_data);
end
     