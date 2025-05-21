function PPI = spm_peb_gppi(varargin)
% Bold deconvolution to create physio- or psycho-physiologic interactions

% CB edits: If Uu is a cell array (for gPPI) the resulting PPI will also be
% a cell array
% also pruned down the interface stuff

% FORMAT PPI = spm_peb_ppi(SPMname,ppiflag,VOI,Uu,ppiname,showGraphics)
%
% SPM          - Structure containing generic details about the analysis or
%                the fully qualified filename of such a structure.
% VOI          - Structure containing details about a VOI (as produced by
%                spm_regions) or the fully qualified filename of such a
%                structure. If a structure, then VOI should be of size 1x1
%                in the case of simple deconvolution, and psychophysiologic
%                interactions) or 1x2, in the case of physiophysiologic
%                interactions. If a file name it should be 1xN or 2xN.
% Uu           - Matrix of input variables and contrast weights. This is an
%                [n x 3] matrix. The first column indexes SPM.Sess.U(i). The
%                second column indexes the name of the input or cause, see
%                SPM.Sess.U(i).name{j}. The third column is the contrast
%                weight. Unless there are parametric effects the second
%                column will generally be a 1.
% ppiname      - Basename of the PPI file to save. The saved file will be:
%                <PATH_TO_SPM.MAT>/PPI_<ppiname>.mat
%
% PPI.ppi      - (PSY*xn  or xn1*xn2) convolved with the HRF
% PPI.Y        - Original BOLD eigenvariate. Use as covariate of no interest
% PPI.P        - PSY convolved with HRF for psychophysiologic interactions,
%                or in the case of physiophysologic interactions contains
%                the eigenvariate of the second region.
% PPI.name     - Name of PPI
% PPI.xY       - Original VOI information
% PPI.xn       - Deconvolved neural signal(s)
% PPI.psy.u    - Psychological variable or input function (PPIs only)
% PPI.psy.w    - Contrast weights for psychological variable (PPIs only)
% PPI.psy.name - Names of psychological conditions (PPIs only)
%__________________________________________________________________________
%
% This routine is effectively a hemodynamic deconvolution using full priors
% and EM to deconvolve the HRF from a hemodynamic time series to give a
% neuronal time series [that can be found in PPI.xn].  This deconvolution
% conforms to Wiener filtering. The neuronal process is then used to form
% PPIs. See help text within function for more details.
%__________________________________________________________________________
% Copyright (C) 2002-2014 Wellcome Trust Centre for Neuroimaging

% Darren Gitelman
% $Id: spm_peb_ppi.m 6556 2015-09-15 15:42:04Z guillaume $


% SETTING UP A PPI THAT ACCOUNTS FOR THE HRF
% =========================================================================
% PPI's were initially conceived as a means of identifying regions whose
% reponses can be explained in terms of an interaction between activity in
% a specified source (the physiological factor) and some experimental
% effect (the psychological factor). However, a problem in setting up PPI's
% is that in order to derive a proper estimate of the interaction between
% a psychological variable (P) and measured hemodynamic signal (x), one
% cannot simply convolve the psychological variable with the hrf (HRF) and
% multiply by the signal. Thus:
%
%                  conv(P,HRF).* x ~= conv((P.*xn),HRF)
%
% P   = psychological variable
% HRF = hemodynamic response function
% xn  = underlying neural signal which in fMRI is convolved with the hrf to
%       give the signal one measures -- x.
% x   = measured fmri signal
%
% It is actually the right hand side of the equation one wants.
% Thus one has to work backwards, in a sense, and deconvolve the hrf
% from x to get xn. This can then be multiplied by P and the resulting
% vector (or matrix) reconvolved with the hrf.
%
% This algorithm uses a least squares strategy to solve for xn.
%
% The source's hemodynamics are x = HRF*xn;
%
% Using the constraint that xn should have a uniform spectral density
% we can expand x in terms of a discrete cosine set (xb)
%
%      xn  = xb*B
%       B  = parameter estimate
%
% The estimator of x is then
%
%       x  = HRF(k,:)*xn
%       x  = HRF(k,:) * xb * B
%
% This accounts for different time resolutions between our hemodynamic
% signal and the discrete representation of the psychological variable. In
% this case k is a vector representing the time resolution of the scans.
%
% Conditional estimates of B allow for priors that ensure uniform variance
% over frequencies.
%
% PPI STATISTICAL MODEL
% =========================================================================
% Once the PPI.ppi interaction term has been calculated a new GLM must be
% setup to search for the interaction effects across the brain. This is
% done using a standard, first level, fMRI model, which must include 3
% covariates, PPI.ppi (interaction), PPI.Y (main effect: source region bold
% signal) and PPI.P (main effect: "psychological" condition), plus any
% nuisance regressors according to the particular design.
%
% NB: Designs that include only the interaction term without the main
% effects are not proper as inferences on the interaction will include a
% mixture of both main and interaction effects.
%
% Once the model has been setup and run, a contrast of [1 0 0] over the
% PPI.ppi, PPI.Y and PPI.P columns respectively, will show regions with a
% positive relationship to the interaction term, discounting any main
% effects. Negative regressions can be examined with [-1 0 0]. A PPI random
% effects analysis would involve taking the con* files from the [1 0 0]
% t-contrast for each subject and forwarding them to a second level
% analysis.



% Check inputs
%--------------------------------------------------------------------------
SPM = varargin{1};
try, swd = SPM.pwd; catch, swd = pwd; end

SPM.swd = spm_file(swd,'cpath');
cwd     = pwd;
cd(SPM.swd)

%======================================================================
p.xY = varargin{2};
xY(1) = p.xY;
Sess  = SPM.Sess(xY(1).Sess);


% Name of PPI file to be saved
%--------------------------------------------------------------------------

Uu = varargin{3};

for g=1:numel(Uu) %Uu needs to be a cell array
    
    % get 'causes' or inputs U
    %----------------------------------------------------------------------
    U.name = {};
    U.u    = [];
    U.w    = [];
    for i = 1:size(Uu{g},1)
        U.u           = [U.u Sess.U(Uu{g}(i,1)).u(33:end,Uu{g}(i,2))];
        U.name{end+1} = Sess.U(Uu{g}(i,1)).name{Uu{g}(i,2)};
        U.w           = [U.w Uu{g}(i,3)];
    end
    
    % Setup variables
    %--------------------------------------------------------------------------
    RT      = SPM.xY.RT;
    dt      = SPM.xBF.dt;
    NT      = round(RT/dt);
    fMRI_T0 = SPM.xBF.T0;
    %N       = length(xY(1).u);
    N       = length(xY(1).mean);
    k       = 1:NT:N*NT;                       % microtime to scan time indices
    
    
    % Setup other output variables
    %--------------------------------------------------------------------------
    PPI{g}.name = varargin{4};
    %PPI{g}.xY = xY;
    PPI{g}.RT = RT;
    PPI{g}.dt = dt;
    
    
    % Create basis functions and hrf in scan time and microtime
    %--------------------------------------------------------------------------
    hrf = spm_hrf(dt);
    
    
    % Create convolved explanatory {Hxb} variables in scan time
    %--------------------------------------------------------------------------
    xb  = spm_dctmtx(N*NT + 128,N);
    Hxb = zeros(N,N);
    for i = 1:N
        Hx       = conv(xb(:,i),hrf);
        Hxb(:,i) = Hx(k + 128);
    end
    xb = xb(129:end,:);
    
    
    % Get confounds (in scan time) and constant term
    %--------------------------------------------------------------------------
    X0 = xY(1).X0;
    M  = size(X0,2);
    
    
    % Get response variable
    %--------------------------------------------------------------------------
    for i = 1:size(xY,2)
        %Y(:,i) = xY(i).u;
        Y(:,i) = xY(i).mean;
    end
    
    
    % Remove confounds and save Y in ouput structure
    %--------------------------------------------------------------------------
    Yc    = Y - X0*inv(X0'*X0)*X0'*Y;
    PPI{g}.Y = Yc(:,1);
    if size(Y,2) == 2
        PPI{g}.P = Yc(:,2);
    end
    
    deconv = 0;
    if deconv
        
        % Specify covariance components; assume neuronal response is white
        % treating confounds as fixed effects
        %--------------------------------------------------------------------------
        Q = speye(N,N)*N/trace(Hxb'*Hxb);
        Q = blkdiag(Q, speye(M,M)*1e6  );
        
        
        % Get whitening matrix (NB: confounds have already been whitened)
        %--------------------------------------------------------------------------
        W = SPM.xX.W(Sess.row,Sess.row);
        
        
        % Create structure for spm_PEB
        %--------------------------------------------------------------------------
        clear P
        P{1}.X = [W*Hxb X0];        % Design matrix for lowest level
        P{1}.C = speye(N,N)/4;      % i.i.d assumptions
        P{2}.X = sparse(N + M,1);   % Design matrix for parameters (0's)
        P{2}.C = Q;
        
        
        %======================================================================
        
        % COMPUTE PSYCHOPHYSIOLOGIC INTERACTIONS
        % use basis set in microtime
        %----------------------------------------------------------------------
        % get parameter estimates and neural signal; beta (C) is in scan time
        % This clever trick allows us to compute the betas in scan time which
        % is much quicker than with the large microtime vectors. Then the betas
        % are applied to a microtime basis set generating the correct neural
        % activity to convolve with the psychological variable in microtime
        %----------------------------------------------------------------------
        %    deconv = 0;
        %    if deconv
        C  = spm_PEB(Y,P);
        xn = xb*C{2}.E(1:N);
        xn = spm_detrend(xn);
    end
    % Setup psychological variable from inputs and contrast weights
    %----------------------------------------------------------------------
    PSY = zeros(N*NT,1);
    for i = 1:size(U.u,2)
        PSY = PSY + full(U.u(:,i) * U.w(i));
    end
    %PSY = spm_detrend(PSY);
    
    % Multiply psychological variable by neural signal
    %----------------------------------------------------------------------
    if deconv
        PSYxn = PSY.*xn;
        
        % Convolve, convert to scan time, and account for slice timing shift
        %----------------------------------------------------------------------
        ppi = conv(PSYxn,hrf);
        ppi = ppi((k-1) + fMRI_T0);
    end
    % Convolve psych effect, convert to scan time, and account for slice
    % timing shift
    %----------------------------------------------------------------------
    PSYHRF = conv(PSY,hrf);
    PSYHRF = PSYHRF((k-1) + fMRI_T0);
    
    % Save psychological variables
    %----------------------------------------------------------------------
    PPI{g}.psy   = U;
    PPI{g}.P     = PSYHRF;
    if deconv
        PPI{g}.xn    = xn;
        PPI{g}.ppi   = spm_detrend(ppi);
    end
    PPI{g}.o_ppi = PSYHRF.*Y;
    
end

% Save
%--------------------------------------------------------------------------
%str    = ['gPPI_' PPI{1}.name '.mat'];

%save(fullfile(SPM.swd,str),'PPI', spm_get_defaults('mat.format'))

%fprintf('   PPI saved as %s\n',spm_file(fullfile(SPM.swd,str)));

% Clean up
%--------------------------------------------------------------------------
cd(cwd);


