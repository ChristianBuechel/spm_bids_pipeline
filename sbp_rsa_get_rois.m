function sbp_rsa_get_rois(all_sub_ids)

% We define which betas should be averaged in RSA by using
% "F-contrasts". We  then simply multiply the beta matrix with this F-contrast
% to average across some conditions. Careful : lss_betas are a subset of all
% betas! We have also adopted the regexp idea of TDT 


% get all infos
[path,vars,~,~,rsa] = get_study_specs;
BIDS                = spm_BIDS(path.preprocDir);

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

% define file names for output
target_f = fullfile(path.secondlevelDir,rsa.ana_name);
if not(isfolder(target_f))
    mkdir(target_f)
end

% setup up classic SPM xY structure with our xY.ind extension
for i=1:numel(roi_ind)
    xY(i).name = xA.labels(i).name;
    xY(i).ind  = xA.labels(i).index;
    xY(i).spec = xA.info.files.images;
    xY(i).def  = 'mask';
    xY(i).xyz  = Inf;
    xY(i).Ic   = NaN; % maybe set to eoi
    xY(i).sk   = 0;   % smoothing kernel
end

n_subs        = numel(all_sub_ids);
run_parallel  = 1; % to do

% loop over subjects
for sub = 1:n_subs
    fprintf('\tdoing subject %d\n',sub);
    sub_id     = all_sub_ids(sub);
    ana_dir    = fullfile(path.firstlevelDir,sprintf('sub-%2.2d',sub_id),rsa.ana_name);
    c          = load(fullfile(ana_dir,'SPM.mat'));
    % add analysis path to beta file names
    if rsa.lss == 1
        for ii=1:numel(c.SPM.Vbeta_lss)
            c.SPM.Vbeta_lss(ii).fname = fullfile(ana_dir,c.SPM.Vbeta_lss(ii).fname);
        end
    else
        for ii=1:numel(c.SPM.Vbeta)
            c.SPM.Vbeta(ii).fname = fullfile(ana_dir,c.SPM.Vbeta(ii).fname);
        end
    end
    
    if numel(c.SPM.Sess) > 1
        error('for safety reasons we just do it for concatenated designs');
    end
    % get all regressor names
    for ii=1:numel(c.SPM.xX.name)
        all_reg_names(ii) = c.SPM.xX.name(ii);
    end
    lss_mask = isfinite(c.SPM.lss_ind(1,:));
    % create F-type contrast matrices
    F_con = [];
    for ii=1:numel(rsa.stims)
        ind = ~cellfun(@isempty,regexp(all_reg_names,rsa.stims(ii)));
        if rsa.lss == 1
            ind = ind(:,lss_mask); %prune down
        end
        fprintf('grouping:\n');disp(strvcat(all_reg_names(ind)));
        F_con(ii,:) = ind ./ sum(ind);
    end
    % get transformation stuff
    epi            = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');
    epi            = epi(1); % if there are brain and spinal, just take brain
    func_dir       = spm_file(epi,'path');
    Vwm_file       = spm_vol(char(spm_file(epi,'prefix','wmeana')));
    Nepi_2_templ   = nifti(fullfile(func_dir, 'y_epi_2_template.nii'));
    good_rois      = zeros(1,numel(roi_ind));
    tic;
    parfor gg=1:numel(roi_ind)
    %for gg=1:numel(roi_ind)
        [dist(gg), good_rois(gg)] = par_extract(gg,roi_names,roi_ind,xY,c,Vwm_file,Nepi_2_templ,rsa,F_con,sub_id);
    end
    toc
    RDMi.dist      = dist;
    RDMi.good_rois = good_rois;
    RDMi.thresh    = rsa.valid;
    if rsa.lss == 1
        save(fullfile(ana_dir,sprintf('lss_RDMi_%s.mat',rsa.atlas)),'RDMi'); % save single subject
    else
        save(fullfile(ana_dir,sprintf('lsa_RDMi_%s.mat',rsa.atlas)),'RDMi'); % save single subject
    end
end
end


function [dist, good_rois] = par_extract(gg,roi_names,roi_ind,xY,c,Vwm_file,Nepi_2_templ,rsa,F_con,sub_id)

Y    = get_roi_ts(xY(gg),c.SPM.xY.VY, Vwm_file, Nepi_2_templ);
fprintf('%s ROI %d/%d with %d voxels. Thereof ',strtrim(roi_names(gg,:)),gg,numel(roi_ind),size(Y,2));
n_cons = numel(rsa.stims); % or use size(F_con,1)

% now read lss betas from file
if rsa.lss
    betas     = get_roi_ts(xY(gg),c.SPM.Vbeta_lss, Vwm_file, Nepi_2_templ);
else
    betas     = get_roi_ts(xY(gg),c.SPM.Vbeta, Vwm_file, Nepi_2_templ);
end
% invalid data (Y==0) or Beta == NaN
voxel_zero      = (sum(Y)==0);
voxel_nan       = ~isfinite(sum(betas));
bad_voxel_mask  = voxel_zero | voxel_nan;
bad_voxels      = sum(bad_voxel_mask);

% prune
betas  = betas(:,~bad_voxel_mask);
Y      = Y(:,~bad_voxel_mask);
fprintf('Zero=%d NaN=%d\t-->\t%d final voxels.',sum(voxel_zero),sum(voxel_nan),size(Y,2));
good_voxels = (numel(bad_voxel_mask)-bad_voxels)./numel(bad_voxel_mask); % in raw EPIs missing voxels are zero, as all values are > 0 we can use sum()
if good_voxels > rsa.valid
    fprintf('\n');
    Y    = spm_filter(c.SPM.xX.K,c.SPM.xX.W*Y); % filter
    res  = spm_sp('r',c.SPM.xX.xKXs,Y);       % get residuals : OK
    
    % beta = c.SPM.xX.pKX*Y; % ordinary least squares estimate of beta = inv(X'*X)*X'*Y
    % res2 = Y - c.SPM.xX.xKXs.X * beta; %remove all --> same as above
    
    [sigma,shrinkage] = covshrink_lw2(res);
    % average as required using F-con and normalize
    ave_beta      = betas'*F_con';
    % ave_beta_norm = ave_beta'*sigma^(-1/2);
    [E,D]          = eig(sigma); % (TDT) --> much faster
    ave_beta_norm  = ave_beta'*E*diag(1./real(sqrt(diag(D))))*E';
    
    dist.RDM     = pdist(ave_beta_norm,'correlation'); % later we can make this a full matrix using squareform()
    dist.RDM_raw = pdist(ave_beta','correlation'); % non-normalized
    good_rois = 1;
else
    fprintf(' Skipped\n');
    dist.RDM     = nan(1,(n_cons.^2-n_cons)/2); % we need nans here to make est of average easier later
    dist.RDM_raw = nan(1,(n_cons.^2-n_cons)/2);
    good_rois = 0;
end
dist.color = [0 0 1];
dist.name  = sprintf('%s | Subject %d',strtrim(roi_names(gg,:)),sub_id);

end