function extract_segment(epi_files, wmean_file, back_file, c_files, interaction, thresh, ex_var, fname)
% Identifies voxels in segments (e.g. c2 and c3 (WM CSF) partitions) above
% thresh (and eroded once)
% c_files is a cellstr
% epi_files is a char array
% interaction n x 2 specifies if one also wants the joint prob of the segments (ie c2.*c3)
% note that thresh then also needs to have an added entry (cave joint probs of segemnst are rather small when unsmoothed try 0.01 or so)
% EPI and segments might be in different spaces (e.g. through nonlin coreg)
% Therefore one has to supply not only the epi_files, but also an epi file normalzed to the space of the segments
% (--> cmean* for nonlin coreg or fieldmap (vdm) correction)
% in addition one needs to provide the deformation (back_file) FROM epi to segments (--> y_epi_2_T1* for nonlin coreg and vdm)
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
    
    data{i}    = spm_read_vols(V_part); % make data available in RAM for spm_erode
    V_part.dt  = [spm_type('float64') spm_platform('bigend')]; % that's how the data is in RAM
    V_part.dat = spm_erode(data{i}); % erode once
    V_part.pinfo(1) = 1; % dtype now double so scale factor = 1
    V_part.pinfo(3) = 0; % read from from .dat
    % I think we could add V_part to xY : xY(i).spec = V_part;
    % and then omit passing V_part to get_roi_ts
    yy         = get_roi_ts(xY(i), V_epi, V_wmean_scan, N_back_scan, V_part);
    yy         = yy(:,any(yy)); % prune out zero time-series

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
        
        yy              = get_roi_ts(xY(i+g), V_epi, V_wmean_scan, N_back_scan,V_part);
        segment(i+g).name = xY(i+g).name;
        if size(yy,2)>=ex_var
            yy = yy(:,any(yy)); % prune out zero time-series
            segment(i+g).data = get_pcs(yy,ex_var); % get principal components
        else
            segment(i+g).data = []; 
            disp('could not get enough voxels at this thresh');
        end
    end
end

%-Save
%==========================================================================
noise_f = spm_file(V_epi(1).fname,'ext','.mat','prefix',fname);
save(noise_f,'segment');
