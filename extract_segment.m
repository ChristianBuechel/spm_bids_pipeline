function extract_segment(epi_files, wmean_file, back_file, c_files, interaction, thresh, ex_var, fname)
% Identifies voxels in segments (e.g. c2 and c3 (WM CSF) partitions) above
% thresh (and eroded once)
% c_files is a cellstr
% epi_files is a char array
% interaction n x 2 specifies if one also wants the joint prob of the segments (ie c2.*c3)
% note that thresh than also needs to have an added entry (cave joint probs of segemnst are rather small when unsmoothed try 0.01 or so)
% EPI and segments are in different spaces (e.g. through nonlin coreg or nonlin coreg + Dartel)
% Therefore one has to supply not only the epi_files, but also an epi file normalzed to the space of the segments
% (--> cmean* for nonlin coreg, wmean* for nonlin coreg + Dartel)
% in addition one needs to provide the deformation (back_file) FROM epi to segments (--> y_mean* for nonlin coreg, y_epi_2_template* for nonlin coreg + Dartel)
% voxels above thresh are then read from the EPIs and PCs (ex_var) extracted
% if ex_var < 1 the number of PCs that explain ex_var of the total variance
% are extracted
% if ex_var >= 1 ex_var PCs are returned
% the time-series are saved to fname in the directrory of the
% EPIs so session specific

%__________________________________________________________________________
%
n_rois = numel(c_files);
[epi_dir,~,~]    = fileparts(epi_files(1,:));
V_epi            = spm_vol(epi_files);
V_wmean_scan     = spm_vol(wmean_file);
N_back_scan      = nifti(back_file);

for i=1:n_rois
    [~,name{i},~] = fileparts(c_files{i});
    V_part     = spm_vol(char(c_files{i}));
    xY(i).name = name{1};
    xY(i).spec = V_part;
    xY(i).def  = 'mask';
    xY(i).xyz  = Inf;
    xY(i).thresh  = thresh(i);
    
    data{i}    = spm_read_vols(V_part); %make data available in RAM for spm_erode
    V_part.dt  = [spm_type('float64') spm_platform('bigend')];
    V_part.dat = spm_erode(data{i}); %erode once
    V_part.pinfo(1) = 1; % scale factor dtype now double so = 1
    V_part.pinfo(3) = 0; % read from .dat
    
    yy         = get_roi_ts(xY(i),V_part, V_epi, V_wmean_scan, N_back_scan);
    segment(i).data = get_pcs(yy,ex_var); % get principal components
    segment(i).name = name{i};
end

if ~isempty(interaction)
    for g=1:size(interaction,1)
        % just reuse V_part
        V_part.dat      = data{interaction(g,1)}.*data{interaction(g,2)};
        xY(i+g).name    = [name{interaction(g,1)} 'x' name{interaction(g,2)}];
        xY(i+g).thresh  = thresh(i+g);
        xY(i+g).def     = 'mask';
        xY(i+g).xyz     = Inf;
        xY(i+g).spec    = V_part; %use last file
        
        yy         = get_roi_ts(xY(i+g),V_part, V_epi, V_wmean_scan, N_back_scan);
        segment(i+g).name = xY(i+g).name;
        if size(yy,2)>=ex_var
            segment(i+g).data = get_pcs(yy,ex_var); % get principal components
        else
            segment(i+g).data = []; % get principal components
            disp('could not get enough voxels at this thresh');
        end
    end
end

%-Save
%==========================================================================
noise_f = spm_file(V_epi(1).fname,'ext','.mat','prefix',fname);
save(noise_f,'segment');


function y = get_roi_ts(xY,V_part, V_epi, V_wmean_scan, N_back_scan)
[~, XYZmm, ~] = spm_ROI(xY, V_wmean_scan); % get the voxels in the image
XYZvoxM       = inv(V_part.mat)*[XYZmm; ones(1,size(XYZmm,2))];
XYZvoxM       = XYZvoxM(1:3,:); %map mm to vox
mask          = spm_get_data(V_part,XYZvoxM); %extract voxel values
prune_ind     = find(mask>xY.thresh); % threshold
XYZmm         = XYZmm(:,prune_ind); % kickout voxels below threshold

% to account for the nonlin coreg we have to sample from the
% deformation field
%[oXYZvox, oXYZmm]  = transform_back([[-28 -57 15]' [31 -53 15]'], V_wmean_scan, V_epi(1), N_back_scan); %just one coord for debugging

[oXYZvox, ~]  = transform_back(XYZmm, V_wmean_scan, V_epi(1), N_back_scan);

y         = spm_get_data(V_epi,oXYZvox); % get relevant EPI time-series
y         = y(:,any(y)); % prune out zero time-series


