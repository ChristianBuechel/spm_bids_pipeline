function sbp_rsa_searchlight(all_sub_ids)
% does a full searchlight analysis with parameters specified in
% get_study_specs.m 
% the analysis is performed in native space, the resulting
% images for distance to candidate RDMs are stored in the firstlevel
% directory (SL*.nii) 
% in addition the images will be warped to template space (wSL*.nii) which then allows
% second level models


% get all infos
[path,vars,~,~,rsa] = get_study_specs;
BIDS                = spm_BIDS(path.preprocDir);

n_subs              = numel(all_sub_ids);

% loop over subjects

for sub=1:n_subs
    warp_files  = '';
    mbi         = 1;
    sub_id      = all_sub_ids(sub);
    m_epi       = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');%ses 1, run 1
    m_epi       = m_epi(1); % if there are brain and spinal, just take brain
    mean_dir    = spm_file(m_epi,'path'); %
    
    ana_dir    = fullfile(path.firstlevelDir,sprintf('sub-%2.2d',sub_id),rsa.ana_name);
    % get SPM
    c          = load(fullfile(ana_dir,'SPM.mat'));
    % add analysis path to beta file names and move to c.SPM.xY.VY (this is from where searchlight samples data))
    c.SPM.VM.fname = fullfile(ana_dir,c.SPM.VM.fname);
    
    if rsa.lss == 1
        for ii=1:numel(c.SPM.Vbeta_lss)
            c.SPM.Vbeta_lss(ii).fname = fullfile(ana_dir,c.SPM.Vbeta_lss(ii).fname);
        end
        sl_options{1}.raw = c.SPM.xY.VY;
        c.SPM.xY.VY = c.SPM.Vbeta_lss;
    else
        for ii=1:numel(c.SPM.Vbeta)
            c.SPM.Vbeta(ii).fname = fullfile(ana_dir,c.SPM.Vbeta(ii).fname);
        end
        sl_options{1}.raw = c.SPM.xY.VY;
        c.SPM.xY.VY = c.SPM.Vbeta;
    end
    if rsa.mnn == 1
        sl_options{1}.ind_beta = 1:numel(c.SPM.xY.VY);
        if isfield(rsa,'mnn_lim')
            c.SPM.xY.VY = [c.SPM.xY.VY sl_options{1}.raw(1:rsa.mnn_lim)'];
            sl_options{1}.ind_raw  = numel(sl_options{1}.ind_beta)+1:rsa.mnn_lim+numel(sl_options{1}.ind_beta); % take reduced amount
        else
            c.SPM.xY.VY = [c.SPM.xY.VY sl_options{1}.raw'];
            sl_options{1}.ind_raw  = numel(sl_options{1}.ind_beta)+1:numel(c.SPM.xY.VY); %take all scans for residuals
        end
        sl_options{1}.SPM      = c.SPM;
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
    sl_options{1}.F_con = F_con;
    sl_options{2}       = rsa;
    
    skip = 1;
    if ~skip
        R = spm_searchlight(c.SPM,rsa.searchopt,@rsa_searchlight,sl_options);
    end
    % now write images
    for m = 1:numel(rsa.models)
        VO = c.SPM.xY.VY(1);
        if rsa.mnn == 1
            if isfield(rsa,'mnn_lim')
                mnn_str = num2str(rsa.mnn_lim);
            else
                mnn_str = 'all';
            end
        else
            mnn_str = 'none';
        end
        VO.fname = fullfile(ana_dir,sprintf('SL_%s_%dmm_mnn_%s.nii',rsa.models{m}.name,rsa.searchopt.spec,mnn_str));
        if ~skip
            spm_write_vol(VO,R{m});
        end
        warp_files = strvcat(warp_files,VO.fname);
    end
    % using nlin coreg + DARTEL or vdm for warping
    matlabbatch{mbi,sub}.spm.util.defs.comp{1}.def = fullfile(mean_dir,'y_epi_2_template.nii');
    matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.fnames = cellstr(warp_files);
    matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
    matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.interp = 4;
    matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.prefix = 'w';
    mbi = mbi + 1;
end
run_local(matlabbatch);
end

function dF = rsa_searchlight(Y,XYZ,varargin)
% Runs RSA search light
va      = varargin;
rsa     = va{1}{2};
options = va{1}{1};

% here we check whether Y is stacked (ie [betas; raw_data] for mnn
if rsa.mnn ==  1
    ave_beta  = Y(options.ind_beta,:)'*options.F_con';
    if isfield(rsa,'mnn_lim')
        newK      = options.SPM.xX.K(1); %assume rsa.mnn_lim is within first run
        newK.row  = 1:rsa.mnn_lim;
        newK.X0   = newK.X0(1:rsa.mnn_lim,:);
        raw       = spm_filter(newK,options.SPM.xX.W(1:numel(options.ind_raw),1:numel(options.ind_raw))*Y(options.ind_raw,:)); % filter
        newxKXs   = options.SPM.xX.xKXs;
        newxKXs.X = newxKXs.X(1:rsa.mnn_lim,:);
        newxKXs.u = newxKXs.u(1:rsa.mnn_lim,:);
        res       = spm_sp('r',newxKXs,raw);       % get residuals : OK
    else
        raw       = spm_filter(options.SPM.xX.K,options.SPM.xX.W*Y(options.ind_raw,:)); % filter
        res       = spm_sp('r',options.SPM.xX.xKXs,raw);       % get residuals : OK
    end
    
    [sigma,shrinkage] = covshrink_lw2(res);
    [E,D]          = eig(sigma); % (TDT) --> much faster
    ave_beta_norm  = ave_beta'*E*diag(1./real(sqrt(diag(D))))*E';
    dist           = pdist(ave_beta_norm,'correlation');
else
    ave_beta       = Y'*options.F_con';
    dist           = pdist(ave_beta','correlation');
end
% get Spearman between data and candidate RDMs
for c_i = 1:numel(rsa.models)
    if size(rsa.models{c_i}.RDM,1) > 1
        c_model_vec = squareform(rsa.models{c_i}.RDM);
    else
        c_model_vec = rsa.models{c_i}.RDM;
    end
    dF(c_i) = atanh(corr(c_model_vec',dist','type','spearman','rows','pairwise')); %Spearman --> Fisher-Z transformed
end
end

function run_local(matlabbatch)

if size(matlabbatch,2) > 1 % we have multiple subjects
    matlabbatch = matlabbatch(:);
end
spm_jobman('initcfg');
spm('defaults', 'FMRI');
spm_jobman('run',matlabbatch);
end
