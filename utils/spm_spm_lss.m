function SPM = spm_spm_lss(SPM, lss_ind)
% fork of spm_spm to quickly do a Mumfordian LS-S analysis from a full
% LS-A design matrix
% Single session only --> not true anymore

% pooling of different events into one regressor (i.e. LSS) is difficult in designs
% where regerssors are quite different eg different durations ... --> solution:
% lss_ind is a matrix where the 1st column depicts an index and the 2nd the
% duration

% The only thing needed is to add the field SPM.lss_ind to SPM
% here all columns for which we want an LSS beta image are numbered
% confounds are marked with NaN eg [1 2 3 4 5 NaN NaN NaN NaN NaN NaN]
% for 5 events and 6 motion paraneters
% this will give lss_beta_0001.nii etc images
% Note that conditions will always be included from left to right, the
% number is used to name the beta image
%
% Caveat: assumes conv(onset(a) + onset(b),hrf) == conv(onset(a),hrf) + conv(onset(b),hrf) 
% as it adds sampled columns of desmtx 
% This might suffer from long TRs
%
% no changes are made to SPM.mat
% Non-sphericity is assumed to be from the full model
% Best use 'None'
%
% FORMAT SPM = spm_spm(SPM)
%
% Required fields of SPM:
%
% xY.VY - nScan x 1 struct array of image handles (see spm_vol)
%         Images must have the same orientation, voxel size and data type
%       - Any scaling should have already been applied via the image handle
%         scalefactors.
%
% xX    - Structure containing design matrix information
%       - Required fields are:
%         xX.X      - Design matrix (raw, not temporally smoothed)
%         xX.name   - cellstr of parameter names corresponding to columns
%                     of design matrix
%       - Optional fields are:
%         xX.K      - cell of session-specific structures (see spm_filter)
%                   - Design & data are pre-multiplied by K
%                     (K*Y = K*X*beta + K*e)
%                   - Note that K should not smooth across block boundaries
%                   - defaults to speye(size(xX.X,1))
%         xX.W      - Optional whitening/weighting matrix used to give
%                     weighted least squares estimates (WLS). If not
%                     specified spm_spm will set this to whiten the data
%                     and render the OLS estimates maximum likelihood
%                     i.e. W*W' = inv(xVi.V).
%
% xVi   - Structure describing intrinsic temporal non-sphericity
%       - Required fields are:
%         xVi.Vi    - array of non-sphericity components
%                   - defaults to {speye(size(xX.X,1))} - i.i.d.
%                   - specifying a cell array of constraints (Qi)
%                     These constraints invoke spm_reml to estimate
%                     hyperparameters assuming V is constant over voxels.
%                     that provide a high precise estimate of xX.V
%       - Optional fields are:
%         xX.V      - Optional non-sphericity matrix.  Cov(e) = sigma^2*V
%                     If not specified spm_spm will compute this using
%                     a 1st pass to identify significant voxels over which
%                     to estimate V.  A 2nd pass is then used to re-estimate
%                     the parameters with WLS and save the ML estimates
%                     (unless xX.W is already specified).
%
% xM    - Structure containing masking information, or a simple column vector
%         of thresholds corresponding to the images in VY [default: -Inf]
%       - If a structure, the required fields are:
%         xM.TH - nVar x nScan matrix of analysis thresholds, one per image
%         xM.I  - Implicit masking (0=>none, 1 => implicit zero/NaN mask)
%         xM.VM - struct array of explicit mask image handles
%       - (empty if no explicit masks)
%               - Explicit mask images are >0 for valid voxels to assess.
%               - Mask images can have any orientation, voxel size or data
%                 type. They are interpolated using nearest neighbour
%                 interpolation to the voxel locations of the data Y.
%       - Note that voxels with constant data (i.e. the same value across
%         scans) are also automatically masked out.
%
% swd   - Directory where the output files will be saved [default: pwd]
%         If exists, it becomes the current working directory.
%
% In addition, global SPM "defaults" variable is used (see spm_defaults):
% 
% stats.<modality>.UFp - critical F-threshold for selecting voxels over 
%                        which the non-sphericity is estimated (if 
%                        required) [default: 0.001]
% 
% stats.maxres         - maximum number of residual images for smoothness
%                        estimation
%
% stats.maxmem         - maximum amount of data processed at a time (in bytes)
%
% modality             - SPM modality {'PET','FMRI','EEG'}
%
%__________________________________________________________________________
%
% spm_spm is the heart of the SPM package. Given image files and a
% General Linear Model, it estimates the model parameters, variance
% hyperparameters, and smoothness of standardised residual fields, writing
% these out to disk in the current working directory for later
% interrogation in the results section. (NB: Existing analyses in the
% current working directory are overwritten).  This directory
% now becomes the working directory for this analysis and all saved
% images are relative to this directory.
%
% The model is expressed via the design matrix (xX.X). The basic model
% at each voxel is of the form is Y = X*B + e, for data Y, design
% matrix X, (unknown) parameters B and residual errors e. The errors
% are assumed to have a normal distribution.
%
% Sometimes confounds (e.g. drift terms in fMRI) are necessary. These
% can be specified directly in the design matrix or implicitly, in terms
% of a residual forming matrix K to give a generalised linear model
% K*Y = K*X*B + K*e.  In fact K can be any matrix (e.g. a convolution
% matrix).
%
% In some instances i.i.d. assumptions about errors do not hold. For
% example, with serially correlated (fMRI) data or correlations among the
% levels of a factor in repeated measures designs. This non-sphericity
% can be specified in terms of components (SPM.xVi.Vi{i}). If specified
% these covariance components will then be estimated with ReML (restricted
% maximum likelihood) hyperparameters. This estimation assumes the same
% non-sphericity for voxels that exceed the global F-threshold. The ReML
% estimates can then be used to whiten the data giving maximum likelihood
% (ML) or Gauss-Markov estimators. This entails a second pass of the data
% with an augmented model K*W*Y = K*W*X*B + K*W*e where W*W' = inv(xVi.V).
% xVi.V is the non-sphericity based on the hyperparameter estimates.
% W is stored in xX.W and cov(K*W*e) in xX.V. The covariance of the
% parameter estimates is then xX.Bcov = pinv(K*W*X)*xX.V*pinv(K*W*X)'.
%
% If you do not want ML estimates but want to use ordinary least squares
% (OLS) then simply set SPM.xX.W to the identity matrix. Any non-sphericity
% V will still be estimated but will be used to adjust the degrees of freedom
% of the ensuing statistics using the Satterthwaite approximation (c.f.
% the Greenhouse-Geisser corrections).
%
% If [non-spherical] variance components Vi are not specified xVi.Vi and
% xVi.V default to the identity matrix (i.e. i.i.d). The parameters are
% then estimated by OLS.  In this instance the OLS and ML estimates are
% the same.
%
% Note that only a single voxel-specific hyperparameter (i.e. variance
% component) is estimated, even if V is not i.i.d.  This means spm_spm
% always implements a fixed-effects model.
% Random effects models can be emulated using a multi-stage procedure:
% This entails summarising the data with contrasts such that the fixed
% effects in a second model on the summary data are those effects of
% interest (i.e. the population effects). This means contrasts are
% re-entered into spm_spm to make an inference (SPM) at the next
% level. At this higher hierarchical level the residual variance for the
% model contains the appropriate variance components from lower levels.
%
% Under the additional assumption that the standardised error fields
% are non-stationary standard Gaussian random fields, results from
% Random field theory can be applied to estimate the significance
% statistic images (SPM's) adjusting p values for the multiple tests
% at all voxels in the search volume. The parameters required for
% this random field correction are the volume, and Lambda, the covariance
% matrix of partial derivatives of the standardised error fields, estimated
% by spm_est_smoothness.
%
%                            ----------------
%
% The volume analysed is the intersection of the threshold masks,
% explicit masks and implicit masks. See spm_spm_ui for further details
% on masking options.
%__________________________________________________________________________
%
% The output of spm_spm takes the form of an SPM.mat file of the analysis
% parameters, and 'float' flat-file images of the parameter and variance
% [hyperparameter] estimates. An 8bit zero-one mask image indicating the
% voxels assessed is also written out, with zero indicating voxels outside
% tha analysed volume.
%
%                            ----------------
%
% The following SPM.fields are set by spm_spm (unless specified)
%
%     xVi.V      - estimated non-sphericity trace(V) = rank(V)
%     xVi.h      - hyperparameters  xVi.V = xVi.h(1)*xVi.Vi{1} + ...
%     xVi.Cy     - spatially whitened <Y*Y'> (used by ReML to estimate h)
%
%                            ----------------
%
%     Vbeta     - struct array of beta image handles (relative)
%     VResMS    - file struct of ResMS image handle  (relative)
%     VM        - file struct of Mask  image handle  (relative)
%
%                            ----------------
%
%     xX.W      - if not specified W*W' = inv(x.Vi.V)
%     xX.V      - V matrix (K*W*Vi*W'*K') = correlations after K*W is applied
%     xX.xKXs   - space structure for K*W*X, the 'filtered and whitened'
%                 design matrix
%               - given as spm_sp('Set',xX.K*xX.W*xX.X) - see spm_sp
%     xX.pKX    - pseudoinverse of K*W*X, computed by spm_sp
%     xX.Bcov   - xX.pKX*xX.V*xX.pKX - variance-covariance matrix of
%                 parameter estimates
%                 (when multiplied by the voxel-specific hyperparameter ResMS
%                 of the parameter estimates (ResSS/xX.trRV = ResMS) )
%     xX.trRV   - trace of R*V
%     xX.trRVRV - trace of RVRV
%     xX.erdf   - effective residual degrees of freedom (trRV^2/trRVRV)
%     xX.nKX    - design matrix (xX.xKXs.X) scaled for display
%                 (see spm_DesMtx('sca',... for details)
%
%                            ----------------
%
%     xVol.M    - 4x4 voxel->mm transformation matrix
%     xVol.iM   - 4x4 mm->voxel transformation matrix
%     xVol.DIM  - image dimensions - column vector (in voxels)
%     xVol.XYZ  - 3 x S vector of in-mask voxel coordinates
%     xVol.S    - Lebesgue measure or volume       (in voxels)
%     xVol.R    - vector of resel counts           (in resels)
%     xVol.FWHM - Smoothness of components - FWHM, (in voxels)
%
%                            ----------------
%
%     xCon      - Contrast structure (created by spm_FcUtil.m)
%     xCon.name - Name of contrast
%     xCon.STAT - 'F', 'T' or 'P' - for F/T-contrast ('P' for PPMs)
%     xCon.c    - (F) Contrast weights
%     xCon.X0   - Reduced design matrix (spans design space under Ho)
%                 It is in the form of a matrix (spm99b) or the
%                 coordinates of this matrix in the orthogonal basis
%                 of xX.X defined in spm_sp.
%     xCon.iX0  - Indicates how contrast was specified:
%                 If by columns for reduced design matrix then iX0 contains
%                 the column indices. Otherwise, it's a string containing
%                 the spm_FcUtil 'Set' action: Usually one of {'c','c+','X0'}
%                 (Usually this is the input argument F_iX0.)
%     xCon.X1o  - Remaining design space (orthogonal to X0).
%                 It is in the form of the coordinates of this matrix in
%                 the orthogonal basis of xX.X defined in spm_sp.
%     xCon.eidf - Effective interest degrees of freedom (numerator df)
%     xCon.Vcon - ...for handle of contrast/ESS image (empty at this stage)
%     xCon.Vspm - ...for handle of SPM image (empty at this stage)
%__________________________________________________________________________
%
% The following images are written to disk:
%
% mask.<ext>                                          - analysis mask image
% 8-bit (uint8) image of zero-s & one's indicating which voxels were
% included in the analysis. This mask image is the intersection of the
% explicit, implicit and threshold masks specified in the xM argument.
% The XYZ matrix contains the voxel coordinates of all voxels in the
% analysis mask. The mask image is included for reference, but is not
% explicitly used by the results section.
%
%                            ----------------
%
% beta_????.<ext>                                     - parameter images
% These are 32-bit (float32) images of the parameter estimates. The image
% files are numbered according to the corresponding column of the
% design matrix. Voxels outside the analysis mask (mask.<ext>) are given
% value NaN.
%
%                            ----------------
%
% ResMS.<ext>                           - estimated residual variance image
% This is a 64-bit (float64) image of the residual variance estimate.
% Voxels outside the analysis mask are given value NaN.
%
%                            ----------------
%
% RPV.<ext>                              - estimated resels per voxel image
% This is a 64-bit (float64) image of the RESELs per voxel estimate.
% Voxels outside the analysis mask are given value 0.  These images
% reflect the nonstationary aspects the spatial autocorrelations.
%
%                            ----------------
%
% ResI_????.<ext>                - standardised residual (temporary) images
% These are 64-bit (float64) images of standardised residuals. At most
% maxres images will be saved and used by spm_est_smoothness, after which
% they will be deleted.
%__________________________________________________________________________
%
% References:
%
% Statistical Parametric Maps in Functional Imaging: A General Linear
% Approach. Friston KJ, Holmes AP, Worsley KJ, Poline JB, Frith CD,
% Frackowiak RSJ. (1995) Human Brain Mapping 2:189-210.
%
% Analysis of fMRI Time-Series Revisited - Again. Worsley KJ, Friston KJ.
% (1995) NeuroImage 2:173-181.
%__________________________________________________________________________
% Copyright (C) 1994-2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston & Guillaume Flandin
% $Id: spm_spm.m 7738 2019-12-02 12:45:37Z guillaume $


SVNid = '$Rev: 7738 $';

%-Say hello
%--------------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SVNid);
spm('Pointer','Watch');

%-Get SPM
%--------------------------------------------------------------------------
if ~nargin
    [P, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
    if ~sts, SPM = []; return; end
    swd      = spm_file(P,'fpath');
    load(fullfile(swd,'SPM.mat'));
    SPM.swd  = swd;
end

if nargin>0 && ischar(SPM)
    load(SPM);
end

if nargin > 1 % we expect an lss_ind vector
    SPM.lss_ind = lss_ind;
end
%==========================================================================
%-            C H E C K   F I L E S   A N D   F O L D E R S
%==========================================================================

% if numel(SPM.Sess)>1 %more than 1 session
%     error('SPM.mat should only have a single session');
% end

%-Change directory
%--------------------------------------------------------------------------
try
    cd(SPM.swd);
catch
    SPM.swd = pwd;
end

%-Check input files
%--------------------------------------------------------------------------
try
    VY = SPM.xY.VY;
catch
    error('Data have not been specified.');
end

for i = 1:numel(VY)
    if ~spm_existfile(VY(i).fname)
        error('File not found: %s',VY(i).fname);
    end
    if ~spm_mesh_detect(VY)
        % Backward compatibility: propagate scaling (see spm_fmri_spm_ui.m)
        VY(i).private.dat.scl_slope = VY(i).pinfo(1);
        VY(i).private.dat.scl_inter = VY(i).pinfo(2);
    end
end

spm_check_orientations(VY);

M       = VY(1).mat;
DIM     = VY(1).dim;
YNaNrep = spm_type(VY(1).dt(1),'nanrep');
if spm_mesh_detect(VY)
    file_ext = '.gii';
    g        = VY(1).private;
    metadata = g.private.metadata;
    name     = {metadata.name};
    if any(ismember(name,'SurfaceID'))
        metadata = metadata(ismember(name,'SurfaceID'));
        metadata = {metadata.name, metadata.value};
    elseif isfield(g,'faces') && ~isempty(g.faces)
        metadata = {'SurfaceID', VY(1).fname};
    else
        error('SurfaceID not found in GIfTI''s metadata.');
    end
    if isempty(spm_file(metadata{2},'path'))
        metadata{2} = fullfile(spm_file(VY(1).fname,'path'),metadata{2});
    end
    SPM.xVol.G = metadata{2};
else
    file_ext = spm_file_ext;
    metadata = {};
end

%-Delete files from previous analyses
%--------------------------------------------------------------------------
% if ~isempty(spm_select('List',SPM.swd,'^lss_mask\..{3}$'))
%     
%     str = {'Current directory contains SPM estimation files:',...
%         'pwd = ',SPM.swd,...
%         'Existing results will be overwritten!'};
%     if spm_input(str,1,'bd','stop|continue',[1,0],1)
%         spm('Pointer','Arrow')
%         return
%     else
%         sw = warning('off','backtrace');
%         warning('Overwriting old results\n\t (pwd = %s) ',SPM.swd);
%         warning(sw);
%         try, SPM.xX  = rmfield(SPM.xX, 'W'); end
%         try,
%             if isfield(SPM.xVi,'Vi') && numel(SPM.xVi.Vi)>1
%                 SPM.xVi = rmfield(SPM.xVi, 'V');
%             end
%         end
%     end
% end

files = {'^lss_mask\..{3}$','^lss_beta_.{4}\..{3}$'};

for i = 1:numel(files)
    j = cellstr(spm_select('FPList',SPM.swd,files{i}));
    for k = 1:numel(j)
        spm_unlink(j{k});
    end
end


%==========================================================================
%-             A N A L Y S I S   P R E L I M I N A R I E S
%==========================================================================

%-Get design
%--------------------------------------------------------------------------
xX             = SPM.xX;
lss_ind        = SPM.lss_ind;
[nScan, nBeta] = size(xX.X);
%lss_cond       = lss_ind(isfinite(lss_ind));
lss_cond       = find(isfinite(lss_ind(1,:)));
nLss           = numel(lss_cond);

%-Get masking settings
%--------------------------------------------------------------------------
if isfield(SPM,'xM')
    xM         = SPM.xM;
else
    xM         = -Inf(nScan,1);
end
if ~isstruct(xM)
    xM         = struct(...
                    'T',  [],...
                    'TH', xM,...
                    'I',  0,...
                    'VM', {[]},...
                    'xs', struct('Masking','analysis threshold'));
end

mask           = true(DIM);

%-Check confounds (xX.K)
%--------------------------------------------------------------------------
if ~isfield(xX,'K')
    xX.K       = 1;
end

%-Get non-sphericity (xVi), otherwise assume i.i.d.
%--------------------------------------------------------------------------
if isfield(SPM,'xVi')
    xVi        = SPM.xVi;
else
    xVi        = struct('form', 'i.i.d.',...
                        'V',    speye(nScan,nScan));
end

%-Evoke ReML for hyperparameter estimation
%--------------------------------------------------------------------------
if ~isfield(xVi,'V')
    SPM.xY.VY  = VY;
    SPM.xM     = xM;
    SPM.xX.K   = xX.K;
    [xVi, am]  = spm_est_non_sphericity(SPM);
    mask       = mask & am;
    spm('FnBanner',mfilename,SVNid);
end

%-Get weight/whitening matrix:  W*W' = inv(V)
%--------------------------------------------------------------------------
if isfield(xX,'W')
    W          = xX.W;
else
    W          = spm_sqrtm(spm_inv(xVi.V));
    W          = W.*(abs(W) > 1e-6);
    xX.W       = sparse(W);
end




%==========================================================================
%-             I N I T I A L I S E   O U T P U T   F I L E S
%==========================================================================

%-Initialise mask file
%--------------------------------------------------------------------------
VM = struct(...
    'fname',   ['lss_mask' file_ext],...
    'dim',     DIM,...
    'dt',      [spm_type('uint8') spm_platform('bigend')],...
    'mat',     M,...
    'pinfo',   [1 0 0]',...
    'descrip', 'spm_spm:lss_mask',...
    metadata{:});
VM = spm_data_hdr_write(VM);

%-Initialise beta files
%--------------------------------------------------------------------------
Vbeta(1:nLss) = deal(struct(...
    'fname',   [],...
    'dim',     DIM,...
    'dt',      [spm_type('float32') spm_platform('bigend')],...
    'mat',     M,...
    'pinfo',   [1 0 0]',...
    'descrip', 'spm_spm:lss_beta',...
    metadata{:}));

for i = 1:nLss
    Vbeta(i).fname   = [sprintf('lss_beta_%04d',lss_cond(i)) file_ext];
    Vbeta(i).descrip = sprintf('spm_spm:lss_beta (%04d) - %s',lss_cond(i),xX.name{i});
end
Vbeta = spm_data_hdr_write(Vbeta);


%==========================================================================
%-               G E N E R A L   L I N E A R   M O D E L
%==========================================================================


%-Get explicit mask(s)
%==========================================================================
for i = 1:numel(xM.VM)
    if ~(isfield(SPM,'xVol') && isfield(SPM.xVol,'G'))
        %-Assume it fits entirely in memory
        C = spm_bsplinc(xM.VM(i), [0 0 0 0 0 0]');
        v = true(DIM);
        [x1,x2] = ndgrid(1:DIM(1),1:DIM(2));
        for x3 = 1:DIM(3)
            M2  = inv(M\xM.VM(i).mat);
            y1 = M2(1,1)*x1+M2(1,2)*x2+(M2(1,3)*x3+M2(1,4));
            y2 = M2(2,1)*x1+M2(2,2)*x2+(M2(2,3)*x3+M2(2,4));
            y3 = M2(3,1)*x1+M2(3,2)*x2+(M2(3,3)*x3+M2(3,4));
            v(:,:,x3) = spm_bsplins(C, y1,y2,y3, [0 0 0 0 0 0]') > 0;
        end
        mask = mask & v;
        clear C v x1 x2 x3 M2 y1 y2 y3
    else
        if spm_mesh_detect(xM.VM(i))
            v = full(xM.VM(i).private.cdata) > 0;
        else
            v = spm_mesh_project(gifti(SPM.xVol.G), xM.VM(i)) > 0;
        end
        mask = mask & v(:);
        clear v
    end
end

%-Split data into chunks
%==========================================================================
chunksize = floor(spm_get_defaults('stats.maxmem') / 8 / nScan);
nbchunks  = ceil(prod(DIM) / chunksize);
chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),prod(DIM)+1);

spm_progress_bar('Init',nbchunks,'LS-S estimation','Chunks');

for i=1:nbchunks
    chunk = chunks(i):chunks(i+1)-1;
    
    %-Report progress
    %======================================================================
    if i > 1, fprintf(repmat(sprintf('\b'),1,72)); end                  %-# 
    fprintf('%-40s: %30s', sprintf('Chunk %3d/%-3d',i,nbchunks),...
                           '...processing');                            %-#
    
    %-Get data & construct analysis mask
    %======================================================================
    Y     = zeros(nScan,numel(chunk));
    cmask = mask(chunk);
    for j=1:nScan
        if ~any(cmask), break, end                 %-Break if empty mask
        
        Y(j,cmask) = spm_data_read(VY(j),chunk(cmask));%-Read chunk of data
        
        cmask(cmask) = Y(j,cmask) > xM.TH(j);      %-Threshold (& NaN) mask
        if xM.I && ~YNaNrep && xM.TH(j) < 0        %-Use implicit mask
            cmask(cmask) = abs(Y(j,cmask)) > eps;
        end
    end
    cmask(cmask) = any(diff(Y(:,cmask),1));        %-Mask constant data
    Y            = Y(:,cmask);                     %-Data within mask
    
    
    
    
    %-Whiten/Weight data and remove filter confounds
    %======================================================================
    KWY          = spm_filter(xX.K,W*Y);
    
    confounds    = xX.X(:,~isfinite(SPM.lss_ind(1,:))); % we put NaNs for confounds
       
    for l=1:nLss
        
        lss_mask = zeros(1,size(xX.X,2));
        lss_mask(lss_cond) = 1; % all in
        lss_mask(lss_cond(l)) = 0; % the one of interest out
        [classes,~,~] = unique(SPM.lss_ind(2,isfinite(SPM.lss_ind(2,:)))); %gets us the # of different durations (Nans behave strange in unique ...)
        [~,~,types] = unique(SPM.lss_ind(2,:)); % get us all different indices some are funny due to NaN but we know how many there are (numel(classes))
        types(types>numel(classes)) = 0; % zero those that are funny 
        collapse = nan(size(xX.X,1),max(types));
        for t = 1:max(types)
            collapse(:,t) = (xX.X*(lss_mask'.*(types == t))); % make as many collapsed regressors as there are durations
        end
        xXX_lss = [xX.X(:,l) collapse  confounds]; % 1st column is regressor of interest then as many collapsed regressors as we have durations and then confounds
        
        %-Design space and projector matrix [pseudoinverse] for WLS
        %--------------------------------------------------------------------------
        xX.xKXs        = spm_sp('Set',spm_filter(xX.K,W*xXX_lss));    % KWX
        xX.xKXs.X      = full(xX.xKXs.X);
        xX.pKX         = spm_sp('x-',xX.xKXs);                     % Projector
        %erdf           = spm_SpUtil('trRV',xX.xKXs);               % error df
        
        %-Use non-sphericity xVi.V to compute [effective] degrees of freedom
        %--------------------------------------------------------------------------
%         xX.V           = spm_filter(xX.K,spm_filter(xX.K,W*xVi.V*W')'); % KWVW'K'
%         [trRV, trRVRV] = spm_SpUtil('trRV',xX.xKXs,xX.V);          % trRV (for X)
%         xX.trRV        = trRV;                                     % <R'*y'*y*R>
%         xX.trRVRV      = trRVRV;                                   %-Satterthwaite
%         xX.erdf        = trRV^2/trRVRV;                            % approximation
%         xX.Bcov        = xX.pKX*xX.V*xX.pKX';                      % Cov(beta)
        
        
        %-Weighted Least Squares estimation
        %======================================================================
        beta         = xX.pKX*KWY;                     %-Parameter estimates
        
        %-Write output files
        %======================================================================
        c            = NaN(numel(chunk),1);
        
        %-Write mask file
        %----------------------------------------------------------------------
        if l == 1 %only once
            mask(chunk)  = cmask;
            VM           = spm_data_write(VM, cmask', chunk);
        end
        %-Write beta file
        %----------------------------------------------------------------------
        c(cmask) = beta(1,:); %it's only the first we want
        Vbeta(l) = spm_data_write(Vbeta(l), c, chunk);
    end
    
    %-Report progress
    %======================================================================
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...done');             %-#
    spm_progress_bar('Set',i);
end

fprintf('\n');                                                          %-#
spm_progress_bar('Clear');

if ~any(mask(:))
    error('Please check your data: There are no inmask voxels.');
end


SPM.Vbeta_lss  = SPM.Vbeta; % copy ...  
SPM.Vbeta_lss(isfinite(SPM.lss_ind(1,:))) = Vbeta; % ...and insert lss_betas
% now prune this to only those betas for which we get an lss_beta
SPM.Vbeta_lss = SPM.Vbeta_lss(1:nLss);

%-Save SPM.mat
%--------------------------------------------------------------------------
fprintf('%-40s: %30s','Saving SPM.mat','...writing');                   %-#
fmt = spm_get_defaults('mat.format');
s = whos('SPM');
if s.bytes > 2147483647, fmt = '-v7.3'; end
save('SPM.mat','SPM', fmt);
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#

%==========================================================================
%-                         S A V E   &   E X I T
%==========================================================================

%-Exit
%--------------------------------------------------------------------------
spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#
