function cb_vasa(SPM, Ic)

if nargin < 2
    Ic = NaN;
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

M   = SPM.xY.VY(1).mat;
DIM = SPM.xY.VY(1).dim(1:min(numel(SPM.xY.VY(1).dim),3));
[nScan, nBeta] = size(SPM.xX.X);
TR  = SPM.xY.RT;

if spm_mesh_detect(SPM.xY.VY)
    file_ext = '.gii';
else
    file_ext = spm_file_ext;
end

% fft stuff
HighCut    = 0.08;
LowCut     = 0.01;

sFreq    =  1/TR;
sLength  =  nScan;
if (LowCut >= sFreq/2)
    idx_LowCut = sLength/2 + 1;
else
    idx_LowCut = ceil(LowCut * sLength * TR + 1);
end
if (HighCut>=sFreq/2)||(HighCut==0)
    idx_HighCut = sLength/2 + 1;
else
    idx_HighCut = fix(HighCut *sLength *TR + 1);
end


%-Initialise VASA image
%--------------------------------------------------------------------------
V_vasa = deal(struct(...
    'fname',    'vasa_res.nii',...
    'dim',      DIM,...
    'dt',       [spm_type('float64') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  'VasA correction'));
V_vasa = spm_data_hdr_write(V_vasa);


%-Loop over chunks
%--------------------------------------------------------------------------
chunksize = floor(spm_get_defaults('stats.maxmem') / 8 / nScan);
nbchunks  = ceil(prod(DIM) / chunksize);
chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),prod(DIM)+1);

spm_progress_bar('Init',nbchunks,'Estimating VASA','Chunks');

for i=1:nbchunks
    chunk = chunks(i):chunks(i+1)-1;
    
    %-Get mask
    %----------------------------------------------------------------------
    m = spm_data_read(SPM.VM,chunk) > 0;
    m = m(:)';
    
    %-Get raw data, whiten and filter
    %----------------------------------------------------------------------
    y = zeros(nScan,numel(chunk));
    for j=1:nScan
        y(j,:) = spm_data_read(SPM.xY.VY(j),chunk);
    end
    y(:,~m) = [];
    
    X0_detrend = SPM.xX.K;
    for j=1:numel(X0_detrend)
        X0_detrend(j).X0 = X0_detrend(j).X0(:,1); % just leave the lowest frequ ie detrend term
    end
    
    %y = spm_filter(X0_detrend,SPM.xX.W*y);
    y = spm_filter(SPM.xX.K,SPM.xX.W*y);
    
    if Ic ~= 0 % if zero remove nothing
        
        %-Parameter estimates: beta = xX.pKX*xX.K*y
        %------------------------------------------------------------------
        beta = zeros(nBeta,numel(chunk));
        for j=1:nBeta
            beta(j,:) = spm_data_read(SPM.Vbeta(j),chunk);
        end
        beta(:,~m) = [];
        
        %-Subtract Y0 = XO*beta,  Y = Yc + Y0 + e
        %------------------------------------------------------------------
        if ~isnan(Ic)
            y = y - spm_FcUtil('Y0',SPM.xCon(Ic),SPM.xX.xKXs,beta); %remove F-con
        else
            y = y - SPM.xX.xKXs.X * beta; %remove all
        end
    end
    % now do fft and get band power
    y      = 2*abs(fft(y))/sLength;
    y_vasa = mean(sqrt(y(idx_LowCut:idx_HighCut,:))); % see Kazan 2017 methods page 4, 3rd paragraph
    
    % write vasa
    yy = NaN(numel(chunk),1);
    yy(m)   = y_vasa;
    V_vasa  = spm_data_write(V_vasa, yy, chunk);
    
end
cd(cwd);
end
