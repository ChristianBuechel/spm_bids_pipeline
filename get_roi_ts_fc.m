function get_roi_ts_fc(SPM, wmean_file, back_file, sk, xY, Uu, prep)

% Versatile routine to extract time courses from sphere, masks etc
% SPM refers to a design, necessary to account for nuisanec variables
% EPI and ROIs are in different spaces (e.g. through nonlin coreg or nonlin coreg + Dartel)
% Therefore one has to supply not only the epi_files (through SPM), but also an epi file
% normalized to the space of the refernce space in which the ROIs are
% defined
% i.e. cmean* for nonlin coreg, wmean* for nonlin coreg + Dartel
% in addition one needs to provide the deformation (back_file) FROM epi to ROI reference space
% i.e. y_mean* for nonlin coreg, y_epi_2_template* for nonlin coreg + Dartel
% the time-series are saved in an xY structure saved to sprintf('VOIs_%d.mat',session)
% sk is a smoothing kernel applied to all images before extraction (can take some time)
% we use an extended xY definition

% classical :
% xY(n).def      - VOI definition [sphere, box, mask, cluster, all]
% xY(n).ind      - if def == mask take only parts of the image with this value (e.g. indexed atlas)
% xY(n).xyz      - centre of VOI {mm} if sphere or box
% xY(n).spec     - VOI definition parameters (e.g. radius in mm for sphere or filename for mask)
% xY(n).str      - description of the VOI (string)

% additional:
% xY(n).Ic       - number of (F) contrast to clean data i.e. EOI F-con
% xY(n).T        - number of T (or F) con to find maximum for 2 pass
% xY(n).rad      - radius to search for max for 2 pass
% xY(n).xyz_o    - original center for 2 pass

% if Uu is defined we also do gPPIs
% if Uu == [] this is skipped
tic;
if nargin < 7
    prep = '';
end
if ~isstruct(SPM)
    swd = spm_file(SPM,'fpath');
    try
        load(fullfile(swd,'SPM.mat'));
        SPM.swd = swd;
    catch
        error(['Cannot read ' fullfile(swd,'SPM.mat')]);
    end
end

try, SPM.swd; catch, SPM.swd = pwd; end
cwd = pwd; cd(SPM.swd);

%-Confounds
%--------------------------------------------------------------------------
X0     = SPM.xX.xKXs.X(:,[SPM.xX.iB SPM.xX.iG]);

V_wmean_scan     = spm_vol(wmean_file);
if isempty(back_file)
    N_back_scan = [];
else
    N_back_scan    = nifti(back_file);
end
V_epi            = SPM.xY.VY;
V_beta           = SPM.Vbeta;
if sk > 0
    V_epi        = spm_smoothto16bit_us(V_epi,sk);
    V_beta       = spm_smoothto16bit_s(V_beta,sk);
end

for i = 1:numel(xY) %loop across all ROIs
    fprintf('Doing ROI: %d\n',i);
    %start with finding the peak in the vicinity (xyz - rad)
    if isfield(xY(i),'T')
        V_T   = SPM.xCon(xY(i).T).Vspm; % T image to search peak
        if sk > 0
            V_T = spm_smoothto16bit_s(V_T,sk);
        end
        xY(i).xyz_o  = xY(i).xyz;
        xY(i).xyz    = get_max_coord(xY(i),V_T, V_wmean_scan, N_back_scan, xY(i).rad);
    else
        xY(i).xyz_o = xY(i).xyz;
    end
    
    [y, xyz]   = get_ts(SPM,xY(i),V_epi, V_wmean_scan, N_back_scan, V_beta);
    y_fin      = (~any(~isfinite(y))) & (~(sum(y)==0));  %if outside mask remove NaNs and columns of zeros!!
%     if sum(~y_fin)>0
%         keyboard; %let's see
%     end
    y        = y(:,y_fin);
            
    for sess = 1:numel(SPM.Sess)
        data(sess).sk    = sk;
        xY(i).sk         = sk;
        if ~isfinite(xY(i).xyz)
            xY(i).xyz = xyz;
        end
        xY(i).Sess       = sess;
        
        
        %-Extract session-specific rows from data and confounds
        %--------------------------------------------------------------------------
        if isfield(SPM,'Sess') && isfield(xY(i),'Sess')
            ii     = SPM.Sess(xY(i).Sess).row;
            ny     = y(ii,:);
            xY(i).X0 = X0(ii,:);
        end
        
        % and add session-specific filter confounds
        %--------------------------------------------------------------------------
        if isfield(SPM,'Sess') && isfield(xY(i),'Sess')
            if numel(SPM.Sess) == 1 && numel(SPM.xX.K) > 1
                xY(i).X0 = [xY(i).X0 blkdiag(SPM.xX.K.X0)]; % concatenated
            else
                xY(i).X0 = [xY(i).X0 SPM.xX.K(xY(i).Sess).X0];
            end
        end
        
        %-Remove null space of X0
        %--------------------------------------------------------------------------
        xY(i).X0         = xY(i).X0(:,any(xY(i).X0));
            
        temp_xY          = xY(i);
        %temp_xY.y        = ny;        
        
        %temp_xY.u        = get_pc(ny);
        temp_xY.mean     = mean(ny,2);
        if ~isempty(Uu)
            temp_xY.PPI = spm_peb_gppi(SPM,temp_xY,Uu,sprintf('%s_%d',xY(i).str,sess));
        end
        data(sess).xY(i) = temp_xY;
    end
end
%now save
mat_f = [SPM.swd filesep prep 'VOIs.mat'];
save(mat_f,'data');

cd(cwd);
toc
end

function xyz_max = get_max_coord(xY,V_T, V_wmean_scan, N_back_scan, rad)
xY.spec = rad; % just for the search of the max
[~, XYZmm, ~] = spm_ROI_ind(xY, V_T); % get the voxels in the image

[oXYZvox, ~]  = transform_back(XYZmm, V_wmean_scan, V_T, N_back_scan, 0);

y         = spm_get_data(V_T,oXYZvox); % get T values and find max
[~,max_T] = max(y);
xyz_max   = XYZmm(:,max_T);
end


function [y, XYZ] = get_ts(SPM,xY,V_epi, V_wmean_scan, N_back_scan, V_beta)
[~, XYZmm, ~] = spm_ROI_ind(xY, V_epi(1)); % get the voxels in the image
XYZ = mean(XYZmm,2);
[oXYZvox, ~]  = transform_back(XYZmm, V_wmean_scan, V_epi(1), N_back_scan, 1);

y         = spm_get_data(V_epi,oXYZvox); % get relevant EPI time-series
y         = y(:,any(y)); % prune out zero time-series
y         = spm_filter(SPM.xX.K,SPM.xX.W*y);
y_fin     = any(isfinite(y)) & sum(y);  %if outside mask remove NaNs and columns of zeros!!
oXYZvox   = oXYZvox(:,y_fin);
y         = y(:,y_fin);


%%simple
%beta      = spm_get_data(V_beta,oXYZvox);
%y         = y - spm_FcUtil('Y0',SPM.xCon(xY.Ic),SPM.xX.xKXs,beta);
%%complex
%-Remove null space of contrast
%--------------------------------------------------------------------------
if xY.Ic ~= 0
    
    %-Parameter estimates: beta = xX.pKX*xX.K*y
    %----------------------------------------------------------------------
    %beta  = spm_data_read(SPM.Vbeta,'xyz',xSPM.XYZ(:,Q));
    beta      = spm_get_data(V_beta,oXYZvox);
    
    %-subtract Y0 = XO*beta,  Y = Yc + Y0 + e
    %----------------------------------------------------------------------
    if ~isnan(xY.Ic)
        y = y - spm_FcUtil('Y0',SPM.xCon(xY.Ic),SPM.xX.xKXs,beta);
    else
        y = y - SPM.xX.xKXs.X * beta;
    end
    
end


%%end complex

end


function Y = get_pc(y)
[m n]   = size(y);
if m > n
    [v s v] = svd(y'*y);
    s       = diag(s);
    v       = v(:,1);
    u       = y*v/sqrt(s(1));
else
    [u s u] = svd(y*y');
    s       = diag(s);
    u       = u(:,1);
    v       = y'*u/sqrt(s(1));
end
d       = sign(sum(v));
u       = u*d;
v       = v*d;
Y       = u*sqrt(s(1)/n);
end

function [oXYZvox, oXYZmm]  = transform_back(XYZmm, Vnormalized, Vorig, Nback, prune)

if nargin < 5
    prune = 1; %prune if not specified
end

if isempty(Nback)
    oXYZmm = XYZmm;
else
    bsplin = [2 2 2]; %order for b-spline interpolation
    XYZvox           = inv(Vnormalized.mat)*[XYZmm; ones(1,size(XYZmm,2))];
    XYZvox           = XYZvox(1:3,:);
    oXYZmm           = zeros(size(XYZvox));
    % get coords in original space
    % let's do bspline interp
    for i=1:3
        c           = spm_bsplinc(squeeze(Nback.dat(:,:,:,1,i)),[bsplin  0 0 0]); % no wrap
        oXYZmm(i,:) = spm_bsplins(c,XYZvox(1,:),XYZvox(2,:),XYZvox(3,:),[bsplin  0 0 0]);
    end
    
end
oXYZvox          = inv(Vorig.mat)*[oXYZmm; ones(1,size(oXYZmm,2))];
oXYZvox          = (oXYZvox(1:3,:));

if (prune == 1)
    oXYZvox          = oXYZvox'; %so that unique 'rows' works
    oXYZmm           = oXYZmm';
    %XYZvox           = XYZvox';
    [oXYZvox, ia, ~] = unique(round(oXYZvox),'stable','rows'); %only get rounded voxel coordinates
    oXYZvox          = oXYZvox'; %and back
    oXYZmm           = oXYZmm(ia,:)';
    %XYZvox           = XYZvox(ia,:)';
end

end



