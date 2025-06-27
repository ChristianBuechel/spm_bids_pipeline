function sbp_analyses(all_sub_ids)
% function sbp_analyses(all_sub_ids)
% performs 1st level and 2nd level analyses of all_sub_ids
%
% calls get_study_specs.m to get various options in the structure analysis.XXX 
% 
% analysis.exclude           : struct array with fields "sub" "ses" and "run" e.g. analysis.exclude(1) = struct('sub',3, 'ses',1,'run',[5 6]); 
%                              exclude runs 5 and 6 in session 1 subject 3; NB only run can be a vector!!
% analysis.sess              : vector of sessions to be analyzed order is preserved e.g. [2 1];
% analysis.prune             : use only those scans for which we have events defined in this analysis (plus 8s after last event+duration)  
% analysis.concatenate       : concatenate runs  
% analysis.ana               : 1 - FIR, 2 - hrf, 3 - lsa
% analysis.name_ana          : name for the analysis (part of directory name)
% analysis.lss               : do least squares separate analysis according to Mumford et al. based on lsa analysis (ana = 3)
% analysis.n_base            : number of basis functions (ie number of FIR bins);
% analysis.bin_size          : bin size for FIR anaylsis in TRs (1 if not sopecified)
% analysis.events            : suffix for the onset TSV file e.g. '_lsa' will use file *_lsa.tsv
% analysis.cond_names        : condition names that the analysis will use, need to be exactly as specified in tsv file
% analysis.p_mod             : cell array with an entry for each condition, listing parametric modulators (more than one possible)
%                            : e.g. for three conditions, where the 2nd will use parametric modulators RT and int {{},{'RT','int'},{}}
%                            : name of the parametric modulator needs to be exactly as in the tsv file (columns 4..n)              
% analysis.t_con             : n t-contrasts as a n x conditions(plus parametric modulators) matrix
% analysis.t_con_names       : cell array of n strings
% analysis.f_con             : cell array of f_con matrices (only performed in 2nd level analyses)
% analysis.f_con_names       : cell array of n strings (only performed in 2nd level analyses)
% analysis.skernel           : smoothing kernel to use can be scalar or 1x3 vector (for anisotropic smoothing) 
% analysis.shift             : add a constant to all onsets (in TRs)
%
%
% analysis.do_model          : create first level model
% analysis.bs                : perform brainstem specific analysis (not fully implemented yet)
% analysis.do_est            : estimate first level model
%
% analysis.do_vasa           : estimate vasa maps
% analysis.do_correct_vasa   : correct con and beta images using vasa estimates
% analysis.use_vasa          : use vasa corrected maps for 2nd level
%
% analysis.do_cons           : estimate contrasts
% analysis.do_warp           : perform appropriate warping to get con/beta images to template space for 2nd level analyses;
% analysis.do_smooth         : smooth warped images 
%
% analysis.noise_corr        : a string that can contain the following words to indicate which noise regressors 
%                              to include: "mov6" "mov24" "physio" "other" "wm" "csf" "roi" "movnoise" (see code for specific forma of the *.txt or *.mat files)
%
% analysis.cvi               : non-sphericity either "none" "AR(1)" "FAST" or "wls" (NB "wls" requires the wls toolbox by Joern Diedrichsen)
% analysis.hpf               : High pass filter in s
% analysis.do_fact           : perform 2nd level anova w/o subject constants
% analysis.do_fact_con       : estimate t-contrasts for group 
% analysis.fact_dept         : non-sphericity 2nd level ANOVA (cov)
% analysis.fact_var          : non-sphericity 2nd level ANOVA (var)
%
% analysis.do_one_t          : perform one sample t-tests for all t-contrasts
%
% one can also define groups of subjects 
% analysis.group_weights     : adds a level of contrasts at the group level (e.g. [1 1;1 -1;-1 1] - tests for common and different patterns for 2 groups)
% analysis.group_ana_names   : name of group levels e.g. {'All','Sporty>Lazy','Lazy>Sporty'};
% analysis.all_subs          : subjects to include
% analysis.group_ind         : vector with indices saying who belongs to which group (e.g. sporty = 1, lazy = 2) same size as analysis.group_ind
% if field "group_ind" does not exist in analysis all subjects are included in a single group 

 
    
% 
%% Prepare everything
% now read in stuff form get_study_specs
[path,vars,analysis]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);
del_all      = 0;

% house keeping
con_temp          = 'con_%04.4d.nii';
beta_temp         = 'beta_%04.4d.nii';
con_templ         = {'con_[0-9][0-9][0-9][0-9]','spmF_[0-9][0-9][0-9][0-9]','spmT_[0-9][0-9][0-9][0-9]','ess_[0-9][0-9][0-9][0-9]'};

TR                = vars.sliceTiming.tr;
parallel          = vars.parallel;

% define what needs to be done
ana               = analysis.ana;  % 1=FIR, 2=HRF 3=LS-A
skernel           = analysis.skernel;

if size(skernel,2) == 1
    skernel = repmat(skernel,1,3);
end
sm_str            = num2str(skernel(1));

if ana == 3
    do_lss        = analysis.lss;  % also do an LSS analysis ?
else
    do_lss        = 0; % never for hrf and fir
end
concatenate     = analysis.concatenate;

do_model        = analysis.do_model;
bs              = analysis.bs;
prune           = analysis.prune;
do_est          = analysis.do_est ;
do_vasa         = analysis.do_vasa;
do_cons         = analysis.do_cons;
do_correct_vasa = analysis.do_correct_vasa;
do_warp         = analysis.do_warp;
do_smooth       = analysis.do_smooth;

% TODO change ppi flag logic
if isfield(analysis,'ppi')
    do_ppi          = analysis.do_ppi; 
    if do_ppi
        ppi_specs = analysis.ppi;
    
        if parallel && ppi_specs.skernel > 0
            warning(['You are about to apply smoothing to the ppi seed time series in parallel mode.' ...
                'This operation requires a significant amount of RAM and may cause your system to crash. Proceed with caution.'])
        end
    end
else
    do_ppi          = false;
end

% second level
do_fact         = analysis.do_fact;
do_fact_con     = analysis.do_fact_con;

do_one_t        = analysis.do_one_t;


if isfield(analysis,'sess')
    sess_2_ana = analysis.sess;
else
    sess_2_ana = 1:vars.nSess;
end

if isfield(analysis,'wRuns')
    ana_runs = analysis.wRuns;
else
    ana_runs = 1:vars.nRuns;
end

if ~isfield(analysis,'group_ind')
    analysis.group_weights = 1;
    analysis.group_ana_names   = {'All'};
    analysis.all_subs      = all_sub_ids;
    analysis.group_ind     = ones(1,numel(all_sub_ids));
end

if ~isempty(analysis.t_con)
    warp_a_con        = 1; % additional cons need to be warped
else
    warp_a_con        = 0; % additional cons need to be warped
end

% other things
noise_corr        = analysis.noise_corr; % the whole lot
cvi               = analysis.cvi;
if strcmp(upper(cvi),'WLS')
    wls = 1;
else
    wls = 0;
end

shift             = analysis.shift; % in TRs

spm_path          = fileparts(which('spm')); %get spm path
mat_name          = which(mfilename);
[~,mat_name,~]    = fileparts(mat_name);

if ana == 1 % FIR
    warp_beta         = 0; 
    warp_s_con        = 1; % simple cons
    
end

if ana == 2 % HRF
    warp_beta         = 0;
    warp_s_con        = 1; % simple cons
end

if ana == 3 % HRF LS-A / LSS
    warp_beta         = 1;
    warp_s_con        = 0; % do not warp simple cons
end

if concatenate % simple cons are not required as we can use the betas, but the other cons should be warped
    warp_beta         = 1;
    warp_s_con        = 0; % not simple cons       
end

n_base            = analysis.n_base; % # of basis functions (FIR model)
name_ana          = analysis.name_ana;


% assemble directory name
anadirname  = [noise_corr '_' name_ana '_' num2str(shift) '_' cvi];

if bs
    anadirname = ['BS_' anadirname];
end

if concatenate
    anadirname = ['Conc_' anadirname];
end

%% Main analyis loop
matlabbatch = [];
for sub = 1:n_subs
    mbi        = 1;
    sub_id     = all_sub_ids(sub);
    if numel(sess_2_ana) > 1
        warning('This routine will collapse all sessions into a single one !')
    end
    cnt        = 0;
    % for now all sessions are converted to runs
    
    %prepare exclusions
    epifiles = [];
    tsvfiles = [];
    
    % harvest epifiles and tsvfiles and honor exclusions
    for i_ses = 1:numel(sess_2_ana)
    ses = sess_2_ana(i_ses);
        for run = ana_runs
            if (~isfield(analysis,'exclude')) || (not_excluded(analysis.exclude,sub_id,ses,run))
                cnt = cnt + 1;
                epi = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
                epi = epi(1); % if there are brain and spinal, just take brain
                epifiles{cnt} = epi;
                tsv = strrep(epi,'_bold.nii','_events.tsv'); % change BOLD filename to proper tsv file name
                tsvfiles{cnt} = char(spm_file(tsv,'suffix', analysis.events)); % ... and add suffix
            end
        end
    end
    n_run      = cnt; 
    n_o_run    = n_run; % keep original # of runs bc they will be set to 1 when we concatenate
    
    % get relevant dirs
    struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'type','T1w');
    struc_file = struc_file(1); % if there are more than 1
    
    % now that we can analyse all sessions seperately, we need the 1st ses, 1st run for the mean epi and the deformations
    m_epi = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',1),'run',sprintf('%02d',1),'task',vars.task,'type','bold');%ses 1, run 1
    m_epi = m_epi(1); % if there are brain and spinal, just take brain
    mean_dir = spm_file(m_epi,'path'); % 
    
    template = []; template_wls = []; % everything will first go into a template struct and added to matlabbatch as needed
    
    template.spm.stats.fmri_spec.timing.units   = 'scans';
    template_wls.spm.tools.rwls.fmri_rwls_spec.timing.units = template.spm.stats.fmri_spec.timing.units; % duplicate everything for WLS
    
    template.spm.stats.fmri_spec.timing.RT      = TR;
    template_wls.spm.tools.rwls.fmri_rwls_spec.timing.RT = template.spm.stats.fmri_spec.timing.RT;
    
    template.spm.stats.fmri_spec.timing.fmri_t  = 16;
    template_wls.spm.tools.rwls.fmri_rwls_spec.timing.fmri_t = template.spm.stats.fmri_spec.timing.fmri_t;
    
    template.spm.stats.fmri_spec.timing.fmri_t0 = round(vars.sliceTiming.refslice./max(vars.sliceTiming.so)*16);
    template_wls.spm.tools.rwls.fmri_rwls_spec.timing.fmri_t0 = template.spm.stats.fmri_spec.timing.fmri_t0;
    
    template.spm.stats.fmri_spec.fact           = struct('name', {}, 'levels', {});
    template_wls.spm.tools.rwls.fmri_rwls_spec.fact = template.spm.stats.fmri_spec.fact;
    
    template.spm.stats.fmri_spec.volt           = 1;
    template_wls.spm.tools.rwls.fmri_rwls_spec.volt = template.spm.stats.fmri_spec.volt;
    
    template.spm.stats.fmri_spec.mthresh        = -Inf;
    template_wls.spm.tools.rwls.fmri_rwls_spec.mthresh = template.spm.stats.fmri_spec.mthresh;
    
    template.spm.stats.fmri_spec.global         = 'None';
    template_wls.spm.tools.rwls.fmri_rwls_spec.global = template.spm.stats.fmri_spec.global;
    
    template.spm.stats.fmri_spec.cvi            = cvi;
    template_wls.spm.tools.rwls.fmri_rwls_spec.cvi = 'wls'; % here we differ
    
    if bs == 1
        bs_file    = fullfile(path.templateDir,vars.brainstemID);
        template.spm.stats.fmri_spec.mask           = cellstr(fullfile(mean_dir,spm_file(spm_file(bs_file,'prefix','bins3c'),'filename'))); % specify which voxels to analyse
        template_wls.spm.tools.rwls.fmri_rwls_spec.mask = template.spm.stats.fmri_spec.mask;
    else
       template.spm.stats.fmri_spec.mask           = cellstr(fullfile(mean_dir,'s3rbrain_mask.nii'));
       template_wls.spm.tools.rwls.fmri_rwls_spec.mask = template.spm.stats.fmri_spec.mask;
    end
    
    
    if ana == 1 % FIR
        if isfield(analysis,'bin_size')
            template.spm.stats.fmri_spec.bases.fir.length = analysis.bin_size*n_base;
        else
            template.spm.stats.fmri_spec.bases.fir.length = TR*n_base;
        end
        template.spm.stats.fmri_spec.bases.fir.order  = n_base;
        
    elseif ana == 2 || ana == 3 % HRF (also for LSA LSS)
        template.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    end
    
    template_wls.spm.tools.rwls.fmri_rwls_spec.bases = template.spm.stats.fmri_spec.bases;


    % directory handling
    if do_ppi
        % we use a temporary directory for the PPI analysis to specify and estimate the normal glm
        temp_dir = fullfile(path.firstlevelDir, sprintf('sub-%02d', sub_id), sprintf('%s_tmp', anadirname));
        a_dir    = temp_dir;

        % Store final PPI directory path
        ppi_dir  = fullfile(path.firstlevelDir, sprintf('sub-%02d', sub_id), sprintf('PPI_%s_%s', ppi_specs.str, anadirname));
    else
        a_dir    = fullfile(path.firstlevelDir,sprintf('sub-%02d',sub_id),anadirname);
    end


    lss_ind  = [];
    cond_use = [];
    z        = [];
    all_nuis = [];
    scan_vec = [];
    
    for run = 1:n_run
        % go through all runs get min, max onset to trim scans and find conditions that are never used
        min_onset = Inf;
        max_onset = -Inf;
        x = spm_load(char(tsvfiles{run})); % load onset *.tsv file; analysis.events allows different tsv files
        if isnumeric(x.trial_type) % the rare case that all conditions are numbers...
            x.trial_type = strtrim(cellstr(num2str(x.trial_type)));
        end
        for ii=1:numel(analysis.cond_names)
            t_ind = find(strcmp(x.trial_type,analysis.cond_names(ii)));
            if ~isempty(t_ind) || (concatenate == 1) % allow empty in some runs when concatenate, bc other runs have the condition
                cond_exist(run,ii) = ~isempty(t_ind); % nevertheless track what we have
                c_onset     = (x.onset(t_ind)')./TR + shift; % times in tsv files are in seconds according to BIDS
                c_dur       = (x.duration(t_ind)')./TR;
                min_onset   = min([min_onset c_onset]); % this can then be used to trim scans
                if ana == 1 %FIR
                    max_onset   = max([max_onset c_onset+template.spm.stats.fmri_spec.bases.fir.length./TR]); % FIR time  
                else
                    max_onset   = max([max_onset c_onset+c_dur]); % add duration !!
                end
            end
        end
        
        func_dir = spm_file(epifiles{run},'path');
        if bs == 1
            scans = spm_select('ExtFPlist', func_dir,spm_file(spm_file(epifiles{run},'filename'),'prefix','^bra'),Inf);
        else
            scans = spm_select('ExtFPlist', func_dir,spm_file(spm_file(epifiles{run},'filename'),'prefix','^ra'),Inf);
        end
        if prune == 1
            scan_range{run} = [max(floor(min_onset),1) min(ceil(max_onset + 8./TR),size(scans,1))]; %leave + 8s at the end,
            fprintf('Sub # %d run %d : orig scans: 1-%d will be pruned to %d-%d \n',sub_id,run,size(scans,1),scan_range{run}(1),scan_range{run}(2));
        else
            scan_range{run} = [1 size(scans,1)];
        end
    end
    
    % fix things for non-existing conditions
    any_cond_exist = any(cond_exist,1); % binary mask for conds we don't have at all
    for ww = find(~any_cond_exist)
        warning(sprintf('Sub # %d misses : %s ',sub_id,char(analysis.cond_names(ww))));
    end
    % make a mask for all regressors that we need (ie zeros for
    % non-existing conds) takes p_mod into account
    i_reg = [];n_reg = [];
    for pp=1:numel(analysis.p_mod)
        n_reg(pp) = n_base + n_base * numel(analysis.p_mod{pp});
        if any_cond_exist(pp)
            i_reg     = [i_reg repmat(1,1,n_reg(pp))];
        else
            i_reg     = [i_reg repmat(0,1,n_reg(pp))];
        end
    end
    conds_done{sub} = find(i_reg); % numbers of regressor we have in this subject

    p_mod           = analysis.p_mod(any_cond_exist); % prune pmods

    if ~do_ppi % TODO: implement handling of missing conds for ppi
        if ~isempty(analysis.t_con)
            if any(~i_reg)
                preserve   = find(~any(analysis.t_con(:,~i_reg),2)); % keep only those cons where we do not have weights in missing conds (incl param mods)
            else
                preserve   = 1:size(analysis.t_con,1);
            end
            t_con      = analysis.t_con(preserve,find(i_reg)); % prune t_con
            t_con_names= analysis.t_con_names(preserve); % ... and their names
        else
            t_con      = []; 
            t_con_names= []; 
        end
    else
        t_con       = analysis.t_con;
        t_con_names = analysis.t_con_names;
    end
    
    cond_names = analysis.cond_names(any_cond_exist); % and prune conditions
    n_cond     = numel(cond_names);
    % from here on the 1st level analysis proceeds as if we only are
    % interested in the reduced amount of conditions
    % However, we have stored important vectors so that we know which
    % subject has which condition for the 2nd level
   
    if isfield(analysis,'f_con_names') && isfield(analysis,'f_con')
        f_con             = analysis.f_con;
        f_con_names       = analysis.f_con_names;
    else
        f_con             = [];
        f_con_names       = [];
    end
    
    for run = 1:n_run
        func_dir = spm_file(epifiles{run},'path');
        
        % then get the nuisance vectors (movement, wm, csf etc)
        if bs == 1
            fm        = spm_file(epifiles{run},'prefix','rp_ba','ext','.txt');
        else
            fm        = spm_file(epifiles{run},'prefix','rp_a','ext','.txt');
        end
        movement  = normit(load(char(fm)));
        movement  = movement(scan_range{run}(1):scan_range{run}(2),:); %cut
        
        mov_final     = normit(movement);
        mov_final_d   = diff(mov_final);
        mov_final_d   = normit([mov_final_d(1,:); mov_final_d]);
        mov_final_2   = normit(mov_final.^2);
        mov_final_d_2 = normit(mov_final_d.^2);
        
        all_nuis{run} = [];
        % now assemble confounds according to the noise_corr string
        if strfind(noise_corr,'mov6')
            all_nuis{run} = [all_nuis{run} mov_final];
        elseif strfind(noise_corr,'mov24')
            all_nuis{run} = [all_nuis{run} mov_final mov_final_2 mov_final_d mov_final_d_2];
        end
        
        % other nuisance variables
        
        if strfind(noise_corr,'physio')
            physio_noise_f = char(spm_file(epifiles{run},'prefix','physio_','ext','.txt'));
            if exist(physio_noise_f,'file')
                physio_noise = load(physio_noise_f);
                physio_noise  = physio_noise(scan_range{run}(1):scan_range{run}(2),:); %cut
                all_nuis{run} = [all_nuis{run} normit(physio_noise)];
            end
        end
        
        if strfind(noise_corr,'other')
            other_f = char(spm_file(epifiles{run},'prefix','other_','ext','.mat'));
            if exist(other_f,'file')
                other_cov = load(other_f);
                other_cov  = other_cov(scan_range{run}(1):scan_range{run}(2),:); %cut                
                all_nuis{run} = [all_nuis{run} normit(other_cov.other)];
            end
        end
        
        if strfind(noise_corr,'wm')
            seg_noise_f = char(spm_file(epifiles{run},'prefix','noise_wm_csf_ra','ext','.mat'));
            if exist(seg_noise_f,'file')
                seg_noise = load(seg_noise_f);
                seg_noise_wm  = normit(seg_noise.segment(1).data(scan_range{run}(1):scan_range{run}(2),:)); %cut                
                all_nuis{run} = [all_nuis{run} seg_noise_wm];
            end
        end
        
        if strfind(noise_corr,'csf')
            if bs == 1
                seg_noise_bs_f = char(spm_file(epifiles{run},'prefix','noise_csf_bra','ext','.mat'));
                if exist(seg_noise_bs_f,'file')
                    seg_noise_bs = load(seg_noise_bs_f);
                    seg_noise_bs_csf  = normit(seg_noise_bs.segment(2).data(scan_range{run}(1):scan_range{run}(2),:)); %cut                                    
                    all_nuis{run} = [all_nuis{run} seg_noise_bs_csf];
                end
            else
                seg_noise_f = char(spm_file(epifiles{run},'prefix','noise_wm_csf_ra','ext','.mat'));
                if exist(seg_noise_f,'file')
                    seg_noise = load(seg_noise_f);
                    seg_noise_csf  = normit(seg_noise.segment(2).data(scan_range{run}(1):scan_range{run}(2),:)); %cut                                    
                    all_nuis{run} = [all_nuis{run} seg_noise_csf];
                end
            end
        end
        
        if strfind(noise_corr,'roi')
            roi_noise_f = char(spm_file(epifiles{run},'prefix','noise_roi_ra','ext','.mat'));
            if exist(roi_noise_f,'file')
                roi_noise = load(roi_noise_f);
                %cut
                noise_roi_1  = normit(roi_noise.roi(1).data(scan_range{run}(1):scan_range{run}(2),:)); %cut                                    
                noise_roi_2  = normit(roi_noise.roi(2).data(scan_range{run}(1):scan_range{run}(2),:)); %cut                                    
                all_nuis{run} = [all_nuis{run} noise_roi_1];
                all_nuis{run} = [all_nuis{run} noise_roi_2];
            end
        end
        
        if strfind(noise_corr,'movnoise')
            mov_noise_f = char(spm_file(epifiles{run},'prefix','noise_mov_ra','ext','.mat'));
            if exist(mov_noise_f,'file')
                mov_noise = load(mov_noise_f);
                mov_noise_1  = normit(mov_noise.mov_reg(scan_range{run}(1):scan_range{run}(2),:)); %cut                                    
                all_nuis{run} = [all_nuis{run} mov_noise_1];
            end
        end
        
        n_nuis         = size(all_nuis{run},2);
        z{run}         = zeros(1,n_nuis); % handy for later contrast def
        
        % now select epis for analysis
        if bs == 1
            scans = spm_select('ExtFPlist', func_dir,spm_file(spm_file(epifiles{run},'filename'),'prefix','^bra'),Inf);
        else
            scans = spm_select('ExtFPlist', func_dir,spm_file(spm_file(epifiles{run},'filename'),'prefix','^ra'),Inf);
        end
        
        scans  = scans(scan_range{run}(1):scan_range{run}(2),:); %cut    
        % create vector of number of scans for possible concatenation
        scan_vec(run) = size(scans,1);
        
        template.spm.stats.fmri_spec.sess(run).scans = cellstr(scans);
        template.spm.stats.fmri_spec.sess(run).multi = {''};
        template.spm.stats.fmri_spec.sess(run).multi_reg = {''};
        template.spm.stats.fmri_spec.sess(run).hpf       = analysis.hpf;
        
        x = spm_load(tsvfiles{run}); % load onset *.tsv file; analysis.events allows different tsv files
        if isnumeric(x.trial_type) % the rare case that all conditions are numbers...
            x.trial_type = strtrim(cellstr(num2str(x.trial_type)));
        end
        RES   = [];cnt = 1;cond_use{run} = []; % initialize some variables
        for ii=1:numel(cond_names)
            t_ind = find(strcmp(x.trial_type,cond_names(ii)));
            if ~isempty(t_ind) || (concatenate == 1)% allow empty in some runs when concatenate, bc other runs have the condition
                cond_use{run} = [cond_use{run} ii]; 
                RES(cnt).name      = cond_names{ii};
                RES(cnt).onset     = (x.onset(t_ind)')./TR + shift; % times in tsv files are in seconds according to BIDS
                RES(cnt).onset     = RES(cnt).onset - (scan_range{run}(1) - 1); %correct for cuts                    
                RES(cnt).duration  = (x.duration(t_ind)')./TR;
                if ana == 1 % FIR always has duration of 0
                    RES(cnt).duration  = zeros(size(x.duration(t_ind)'));
                end
                
                % now add parametric regressors
                RES(cnt).orth = 0; % no orthogonalization
                for pm = 1:numel(p_mod{ii})
                    eval(sprintf('t_pm = x.%s;',char(p_mod{ii}(pm)))); % load according to p_mod
                    if ~isnumeric(t_pm) % if a column has no NaNs spm_load returns a numeric vector, else strings
                        t_pm = str2num(strvcat(t_pm)); % convert to numbers
                    end
                    RES(cnt).pmod(pm).name  = char(p_mod{ii}(pm));
                    RES(cnt).pmod(pm).param = t_pm(t_ind); % and add
                    RES(cnt).pmod(pm).poly  = 1;
                end
                cnt = cnt + 1;
            end
        end
        % expand all conditions for LS-A type analysis
        all = [];
        for r = 1:numel(RES)
            new = []; % help contruct to hold all events
            new(:,1) = RES(r).onset';
            new(:,2) = repmat(r,size(new(:,1),1),1);
            new(:,3) = RES(r).duration';
            all = [all ; new];
        end
        [~,i] = sort(all(:,1)); % sort all events by scans/time (or comment these 2 lines)
        all   = all(i,:);
        RES_LSA = [];
        for j=1:size(all,1) % now put in alternative structure
            RES_LSA(j).name       = RES(all(j,2)).name;
            RES_LSA(j).onset      = all(j,1);
            RES_LSA(j).duration   = all(j,3);
        end
        
        % now put RES or RES_LSA in batch structure
        if ana == 1 || ana == 2 % FIR or HRF
            template.spm.stats.fmri_spec.sess(run).cond = RES;
        elseif ana == 3 % HRF LS-A
            template.spm.stats.fmri_spec.sess(run).cond = RES_LSA;
            cond_use{run}  = [1:size(all,1)]; % recalc cond_use
            if run>1
                cond_use{run} = cond_use{run} + max(cond_use{run-1});
            end
        end
        
        %template.spm.stats.fmri_spec.sess(run).regress = []; %so it exists
        % --> this seems to give new problems ???
        
        for nuis = 1:n_nuis % add nuisance variables to batch
            template.spm.stats.fmri_spec.sess(run).regress(nuis) = struct('name', cellstr(num2str(nuis)), 'val', all_nuis{run}(:,nuis));
        end
        if ana == 3 % only necessary if we do LS-A
            if run>1
                lss_ind = [lss_ind [max(lss_ind(1,:))+1:max(lss_ind(1,:))+numel(template.spm.stats.fmri_spec.sess(run).cond) ones(1,numel(z{run})).*NaN;...
                    cat(1,template.spm.stats.fmri_spec.sess(run).cond.duration)' ones(1,numel(z{run})).*NaN]]; %increase numbers
            else
                lss_ind = [lss_ind [1:numel(template.spm.stats.fmri_spec.sess(run).cond) ones(1,numel(z{run})).*NaN;...
                    cat(1,template.spm.stats.fmri_spec.sess(run).cond.duration)' ones(1,numel(z{run})).*NaN]];
            end
        end
    end
    lss_ind = [lss_ind ones(2,n_run).*NaN];  % add constants
    
    template_wls.spm.tools.rwls.fmri_rwls_spec.sess = template.spm.stats.fmri_spec.sess;
    %% do concatenation
    if concatenate
        % now take template struc and create one single run
        c_scan_vec = cumsum(scan_vec);
        new_t_wls  = template_wls; % duplicate template
        new_t      = template;
        
        new_t.spm.stats.fmri_spec.sess(2:end) = []; % get rid of original runs 2..n
        new_t_wls.spm.tools.rwls.fmri_rwls_spec.sess(2:end) = []; % get rid of original runs 2..n
        for run = 2:n_run % go through these runs to put everything in run 1
            % first do scans
            new_t.spm.stats.fmri_spec.sess(1).scans = [new_t.spm.stats.fmri_spec.sess(1).scans;template.spm.stats.fmri_spec.sess(run).scans];
            new_t_wls.spm.tools.rwls.fmri_rwls_spec.sess(1).scans = [new_t_wls.spm.tools.rwls.fmri_rwls_spec.sess(1).scans;template_wls.spm.tools.rwls.fmri_rwls_spec.sess(run).scans];
            
            % now do nuisance variables (e.g. motion etc)
            for re = 1:numel(new_t.spm.stats.fmri_spec.sess(1).regress)
                if numel(new_t.spm.stats.fmri_spec.sess(1).regress) == numel(template.spm.stats.fmri_spec.sess(run).regress)
                    new_t.spm.stats.fmri_spec.sess(1).regress(re).val = [new_t.spm.stats.fmri_spec.sess(1).regress(re).val; template.spm.stats.fmri_spec.sess(run).regress(re).val];
                    new_t_wls.spm.tools.rwls.fmri_rwls_spec.sess(1).regress(re).val = [new_t_wls.spm.tools.rwls.fmri_rwls_spec.sess(1).regress(re).val; template_wls.spm.tools.rwls.fmri_rwls_spec.sess(run).regress(re).val];
                else
                    error('number of nuisance variables needs to be identical'); % in the future this could be relaxed
                end
            end
            % next do conditions
            l_c = numel(template.spm.stats.fmri_spec.sess(run).cond);
            for c=1:l_c
                if ana == 3 %LS A
                    new_t.spm.stats.fmri_spec.sess(1).cond(end+1).name   = template.spm.stats.fmri_spec.sess(run).cond(c).name; %end is now end+1 ... be careful
                    new_t.spm.stats.fmri_spec.sess(1).cond(end).onset    = template.spm.stats.fmri_spec.sess(run).cond(c).onset + c_scan_vec(run-1);
                    new_t.spm.stats.fmri_spec.sess(1).cond(end).duration = template.spm.stats.fmri_spec.sess(run).cond(c).duration;
                    
                    new_t_wls.spm.tools.rwls.fmri_rwls_spec.sess(1).cond(end+1).name   = template_wls.spm.tools.rwls.fmri_rwls_spec.sess(run).cond(c).name; %end is now end+1 ... be careful
                    new_t_wls.spm.tools.rwls.fmri_rwls_spec.sess(1).cond(end).onset    = template_wls.spm.tools.rwls.fmri_rwls_spec.sess(run).cond(c).onset + c_scan_vec(run-1);
                    new_t_wls.spm.tools.rwls.fmri_rwls_spec.sess(1).cond(end).duration = template_wls.spm.tools.rwls.fmri_rwls_spec.sess(run).cond(c).duration;
                    % here we ignore param modulators as they do not make sense
                else
                    new_t.spm.stats.fmri_spec.sess(1).cond(c).onset     = [new_t.spm.stats.fmri_spec.sess(1).cond(c).onset template.spm.stats.fmri_spec.sess(run).cond(c).onset + c_scan_vec(run-1)];
                    new_t_wls.spm.tools.rwls.fmri_rwls_spec.sess(1).cond(c).onset     = [new_t_wls.spm.tools.rwls.fmri_rwls_spec.sess(1).cond(c).onset template_wls.spm.tools.rwls.fmri_rwls_spec.sess(run).cond(c).onset + c_scan_vec(run-1)];
                    
                    new_t.spm.stats.fmri_spec.sess(1).cond(c).duration  = [new_t.spm.stats.fmri_spec.sess(1).cond(c).duration template.spm.stats.fmri_spec.sess(run).cond(c).duration];
                    new_t_wls.spm.tools.rwls.fmri_rwls_spec.sess(1).cond(c).duration  = [new_t_wls.spm.tools.rwls.fmri_rwls_spec.sess(1).cond(c).duration template_wls.spm.tools.rwls.fmri_rwls_spec.sess(run).cond(c).duration];
                    if isfield(new_t.spm.stats.fmri_spec.sess(1).cond(c),'pmod')
                        for pm = 1:numel(new_t.spm.stats.fmri_spec.sess(1).cond(c).pmod)
                            new_t.spm.stats.fmri_spec.sess(1).cond(c).pmod(pm).param  = [new_t.spm.stats.fmri_spec.sess(1).cond(c).pmod(pm).param; template.spm.stats.fmri_spec.sess(run).cond(c).pmod(pm).param];
                            new_t_wls.spm.tools.rwls.fmri_rwls_spec.sess(1).cond(c).pmod(pm).param  = [new_t_wls.spm.tools.rwls.fmri_rwls_spec.sess(1).cond(c).pmod(pm).param; template_wls.spm.tools.rwls.fmri_rwls_spec.sess(run).cond(c).pmod(pm).param];
                        end
                    end
                end
            end
        end
        if ana == 3 % define for which events we want LS-S estimates (and specify their duration as spm_spm_lss will create different regressors according to duration)
            lss_ind  = [1:numel(new_t.spm.stats.fmri_spec.sess(1).cond) ones(1,n_run + numel(z{1})).*NaN;...
                cat(1,new_t.spm.stats.fmri_spec.sess(1).cond.duration)' ones(1,n_run + numel(z{1})).*NaN];
        end
        n_run    = 1; % set # of runs to 1
        template_wls = new_t_wls;
        template     = new_t; % replace template
        z            = z(1); % adjust nuisance counter
    end
    
    if do_model
        if exist(a_dir,'dir') == 7
            str = {[a_dir ' : contains SPM estimation files']};
            if del_all
                ret = 1;
            else
                ret = spm_input(str,1,'bd','abort mission|del dir|del ALL dirs',[0,1,2],1);
            end
            
            if ret == 0
                spm('Pointer','Arrow')
                return
            end
            if ret == 1
                rmdir(a_dir,'s'); % delete all
            end
            if ret == 2
                del_all = 1;
                rmdir(a_dir,'s'); % delete all
            end
        end
        mkdir(a_dir);
        copyfile(which(mfilename),fullfile(a_dir,spm_file(spm_file(which(mfilename),'filename'),'ext','bak'))); % copy this file into analysis directory as *.bak
        copyfile(which('get_study_specs'),fullfile(a_dir,spm_file(spm_file(which('get_study_specs'),'filename'),'ext','bak'))); % and copy this file as it contains everything necessary to repeat the analysis
        for ts = 1:numel(tsvfiles)
            copyfile(tsvfiles{ts},a_dir); % and copy the respective *.tsv files that were used
        end
        template.spm.stats.fmri_spec.dir = {[a_dir]};
        template_wls.spm.tools.rwls.fmri_rwls_spec.dir = template.spm.stats.fmri_spec.dir;
        
        if wls % decide whether we do WLS or OLS
            matlabbatch{mbi,sub} = template_wls;
        else
            matlabbatch{mbi,sub} = template;
        end
        mbi = mbi + 1;
        
        if concatenate
            matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = fullfile(a_dir,'SPM.mat');
            matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.evaluated = scan_vec;
            matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
            matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'spm_fmri_concatenate';
            mbi = mbi + 1;
        end
    end

    if do_est
        if wls
            matlabbatch{mbi,sub}.spm.tools.rwls.fmri_rwls_est.spmmat = {fullfile(a_dir,'SPM.mat')};
            matlabbatch{mbi,sub}.spm.tools.rwls.fmri_rwls_est.method.Classical = 1;
        else
            matlabbatch{mbi,sub}.spm.stats.fmri_est.spmmat           = {fullfile(a_dir,'SPM.mat')};
            matlabbatch{mbi,sub}.spm.stats.fmri_est.method.Classical = 1;
        end
        mbi = mbi + 1;
    end


    if do_lss
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = fullfile(a_dir,'SPM.mat');
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.evaluated = lss_ind;
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'spm_spm_lss';
        mbi = mbi + 1;
    end
    
    if do_vasa
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = fullfile(a_dir,'SPM.mat');
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'cb_vasa';
        mbi = mbi + 1;
    end


    %% PPI

    if do_ppi && ~(do_fact || do_fact_con || do_one_t) %we only do this on the first level


        % first create f contrast for normal glm, we need this to clean the ppi regressors
        % (do_ppi set to false so that contrast is not created for ppi)
        n_ppi_conds      = numel(analysis.ppi.conds); 
        [glm_f_vec, ~]   = create_f_vec(n_run, cond_use, n_base, n_nuis, n_o_run, false, n_ppi_conds, p_mod, ana);

        % TODO: here we can also have users define t-contrasts (based on
        % the original analysis) that can be used to find activtation centers

        matlabbatch{mbi, sub}.spm.stats.con.spmmat = {fullfile(a_dir,'SPM.mat')};
        matlabbatch{mbi, sub}.spm.stats.con.delete = 1;
        matlabbatch{mbi, sub}.spm.stats.con.consess{1}.fcon.name  = 'eff_of_int';
        matlabbatch{mbi, sub}.spm.stats.con.consess{1}.fcon.convec = {glm_f_vec};
        mbi = mbi + 1;
        
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.string    = fullfile(a_dir,'SPM.mat');
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.evaluated = n_cond; % we need the original number of conditions here
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{3}.evaluated = ppi_specs;
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{4}.string    = char(fullfile(mean_dir,'y_epi_2_template.nii')); %TODO: use BIDS instead of function input
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{5}.string    = char(spm_file(m_epi,'prefix','wmeana'));
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'sbp_ppi';
        mbi = mbi + 1;

    end
    
    
    %% now prepare contrasts
    template = [];
    
    % ppi specific overwrites
    if do_ppi
        a_dir           = ppi_dir;

        % set up ppi regressors  
        n_ppi_conds      = numel(analysis.ppi.conds);
        n_ppi_regressors = n_ppi_conds*2 + 1; % 1 psy + 1 ppi regressor for every condition + seed time series (= n_ppi_conds*2+1)
        conds_done{sub}  = 1:n_ppi_regressors; %TODO: account for missing conds
        n_cond = n_ppi_regressors;

        % Get the actual condition names for PPI regressors
        ppi_cond_names = analysis.cond_names(ppi_specs.conds);

        % Create regressor names
        cond_names = [
            strcat('PSY_', ppi_cond_names), ...  % Psychological regressors
            {sprintf('Y_%s', analysis.ppi.name)}, ... % Seed time series 
            strcat('PPI_', ppi_cond_names)  % PPI interaction regressors
        ];

    else
        a_dir = fullfile(path.firstlevelDir, sprintf('sub-%02d', sub_id), anadirname);

    end 
    

    % create f-contrast
    template.spm.stats.con.spmmat = {fullfile(a_dir,'SPM.mat')};
    template.spm.stats.con.delete = 1;
    fco = 1; % counter for f-contrasts
    template.spm.stats.con.consess{fco}.fcon.name  = 'eff_of_int';
    
    if (ana == 3) & (concatenate == 1)
        cond_use{1} = 1:numel(new_t.spm.stats.fmri_spec.sess(1).cond);
    end
    
    [f_vec, n_reg] = create_f_vec(n_run, cond_use, n_base, n_nuis, n_o_run, do_ppi, 0, p_mod, ana);
    template.spm.stats.con.consess{fco}.fcon.convec = {f_vec};

    % construct the effects of int contrast across runs (complicated because some conds might be missing in some runs)
    % get some helpful vectors first
    % i_reg = [];n_reg = [];
    % for pp=1:numel(p_mod)
    %     n_reg(pp) = 1+numel(p_mod{pp});
    %     i_reg     = [i_reg repmat(pp,1,n_reg(pp))];
    % end
    
    % if ana == 3 % f-con for LSA needs adjustment
    %     max_cond = 0;
    %     for run = 1:n_run
    %         max_cond = max(max_cond, max(cond_use{run}));
    %     end
    %     i_reg = 1:max_cond;
    %     n_reg = ones(size(i_reg));
    % end
    

    % now do the simple cons for 2nd level anova
    % simply go through the eoi (f_vec) li-+ne by line and correct by # of occurrence
    ind_t = 1;ind_b = 1;
    %    if warp_s_con || warp_a_con
    simple_con_ind{sub} = [];simple_beta_ind{sub} = [];
    all_t_con_names = [];
    simple_con = zeros(n_cond*n_base,size(f_vec,2));

    for co = 1:n_cond
        for pm = 1:n_reg(co)
            if pm == 1
                label = [cond_names{co} '_main' ]; % the event itself
            else
                label = [cond_names{co} '_' p_mod{co}{pm-1}]; % the parametric regressors
            end
            for i_fir = 1:n_base % take care of FIR bins
                simple_beta_ind{sub} = [simple_beta_ind{sub} ind_b];
                simple_con(ind_b,:) = f_vec(ind_b,:)./sum(f_vec(ind_b,:)); % scale to 1
                ind_b = ind_b + 1;
                
                if warp_s_con
                    all_t_con_names{ind_t}  = [label '_' num2str(i_fir)];
                    simple_con_ind{sub}  = [simple_con_ind{sub} ind_t+fco];
                    template.spm.stats.con.consess{ind_t+fco}.tcon.name    = [label '_' num2str(i_fir)];
                    template.spm.stats.con.consess{ind_t+fco}.tcon.convec  = simple_con(ind_t,:);
                    template.spm.stats.con.consess{ind_t+fco}.tcon.sessrep = 'none';
                    ind_t = ind_t + 1;
                end
            end
        end
    end
    %    end
    
    % additional t-constrasts as specified in get_study_specs, these will also be used as one sample t-tests
    if ~isempty(t_con)
%         if concatenate
%             ind = 1; % reset counter as we have not created those t-cons
%             all_t_con_names = [];
%         end
        add_con_ind = [];
        for co = 1:numel(t_con_names)
            template.spm.stats.con.consess{ind_t+fco}.tcon.name    = [t_con_names{co}];
            all_t_con_names{ind_t}                                 = [t_con_names{co}];
            template.spm.stats.con.consess{ind_t+fco}.tcon.convec  = (simple_con'*t_con(co,:)')'; % nicely takes the weights from the simple t cons into account
            template.spm.stats.con.consess{ind_t+fco}.tcon.sessrep = 'none';
            add_con_ind  = [add_con_ind ind_t+fco];
            ind_t = ind_t + 1;
        end
    end
    
    if do_cons
        % get old con files to delete them *con* *spmT* *spmF* *ess*
        old_con_files = [''];
        for cc = 1:numel(con_templ)
            old_con_files = strvcat(old_con_files,spm_select('FPList',a_dir,con_templ{cc}));
        end
        matlabbatch{mbi,sub}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(old_con_files);
        matlabbatch{mbi,sub}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false; % delete all the stuff we do not need anymore
        mbi = mbi + 1;
        
        matlabbatch{mbi,sub} = template; % now add constrasts to batch
        mbi = mbi + 1;
    end
    
    % prepare_warp
    template    = [];
    warp_files  = '';
    if warp_a_con || warp_s_con        
        for co = 1:size(all_t_con_names,2)
            warp_files(co,:) = [a_dir filesep sprintf(con_temp,co+fco)]; % need to take f cons into account
        end
    end
    
    if warp_beta
        for be = find(sum(f_vec)) % these are all the beta files that are relevant
            warp_files = strvcat(warp_files,[a_dir filesep sprintf(beta_temp,be)]);
        end
    end
    
    warp_files = strvcat(warp_files,fullfile(a_dir,'mask.nii')); % lets also warp the mask
    
    if do_correct_vasa
        vasa_file = fullfile(a_dir,'vasa_res.nii');
        for co = 1:size(warp_files,1)
            c_f = strtrim(warp_files(co,:));
            matlabbatch{mbi,sub}.spm.util.imcalc.input  = {c_f,vasa_file}';
            matlabbatch{mbi,sub}.spm.util.imcalc.output = spm_file(c_f,'prefix','v');
            matlabbatch{mbi,sub}.spm.util.imcalc.outdir = {''};
            matlabbatch{mbi,sub}.spm.util.imcalc.expression = 'i1./i2'; % simply scale by vasa image
            matlabbatch{mbi,sub}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{mbi,sub}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{mbi,sub}.spm.util.imcalc.options.mask = 0;
            matlabbatch{mbi,sub}.spm.util.imcalc.options.interp = 1;
            matlabbatch{mbi,sub}.spm.util.imcalc.options.dtype = 16;
            mbi = mbi + 1;
        end
    end

    

    %% Warping and smoothing

    if do_warp
        % using nlin coreg + DARTEL or vdm corrected EPIs
        matlabbatch{mbi,sub}.spm.util.defs.comp{1}.def = fullfile(mean_dir,'y_epi_2_template.nii');
        if analysis.use_vasa == 1
            matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.fnames = cellstr(spm_file(warp_files, 'prefix','v'));
        else
            matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.fnames = cellstr(spm_file(warp_files, 'prefix',''));
        end
        matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
        matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.interp = 4;
        matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.prefix = 'w';
        mbi = mbi + 1;
    end
    
    if do_smooth
        if analysis.use_vasa == 1
            matlabbatch{mbi,sub}.spm.spatial.smooth.data   = cellstr(spm_file(warp_files, 'prefix','wv'));
        else
            matlabbatch{mbi,sub}.spm.spatial.smooth.data   = cellstr(spm_file(warp_files, 'prefix','w'));
        end
        matlabbatch{mbi,sub}.spm.spatial.smooth.fwhm   = skernel;
        matlabbatch{mbi,sub}.spm.spatial.smooth.prefix = ['s' sm_str];
        mbi = mbi + 1;
    end

    
end

if ~isempty(matlabbatch)
    % run matlabbatch

    n_procs = n_subs; % to not block to many cores on the server
    
    if n_procs > vars.max_procs
        n_procs = vars.max_procs;
    end
    
    if parallel == 1
        run_spm_parallel(matlabbatch, n_procs);
    else
        run_local(matlabbatch);
    end
end


%% second level analyses
% would be interesting to know what we really need from the above
% and whether we can modularize it

if do_ppi
    anadirname = sprintf('PPI_%s_%s', ppi_specs.str, anadirname);
end

n_groups     = max(analysis.group_ind);
for g = 1:n_groups
    [groups{g},~,groups_ind{g}]  = intersect(analysis.all_subs(analysis.group_ind==g),all_sub_ids); % we need the indicies later
end
group_weights = analysis.group_weights;
group_ana_names   = analysis.group_ana_names;

if analysis.use_vasa
    addon = ['V_'];
else
    addon = [];
end

if do_one_t
    if size(groups,2)> 1
        warning('One sample t-test only works for a single group, will only analyse group 1');
    end
    addon   = [addon 'ONE_T'];
    % first let's get rid of negative contrasts --> commented (see below) as it removes
    % too many cons (e.g. [1 0] and [1 -1] are correlated with -1 ...)
    cc    = corrcoef(t_con');
    [~,c] = find(triu(cc)<-.9999999999); % all in c are redundant
    sa    = find(sum(t_con')~=0); % those that do not sum to 0 , e.g. [1 0]
    c     = setdiff(c,sa); % take sa out of c
    pruned_t_con_names = t_con_names;
    pruned_t_con_names(c) = [];
    pruned_add_con_ind = add_con_ind;
    pruned_add_con_ind(c) = [];
    for co = 1:numel(pruned_t_con_names)
        out_dir = fullfile(path.secondlevelDir,[addon '_' anadirname '_' sm_str '_' pruned_t_con_names{co}]);
        out_dir = strrep(out_dir,'>','_bt_'); % these chars are not allowed in directory names ...
        out_dir = strrep(out_dir,'<','_st_');
        mkdir(out_dir);
        
        matlabbatch = [];
        allfiles    = [];
        
        for su = 1:numel(groups{1})
            name        = sprintf('sub-%02d',groups{1}(su));
            a_dir       = fullfile(path.firstlevelDir, name, anadirname);
            s_string    = sprintf('s%s',sm_str);
            if analysis.skernel == 0;s_string = '';end
            if warp_a_con
                swcon_file = fullfile(a_dir, sprintf(con_temp,pruned_add_con_ind(co)));
            end
            
            if analysis.use_vasa == 1
                swcon_file = spm_file(swcon_file,'prefix',sprintf('%swv',s_string));
            else
                swcon_file = spm_file(swcon_file,'prefix',sprintf('%sw',s_string));
            end
            allfiles    = strvcat(allfiles,swcon_file);
        end
        
        mbi = 1;
        matlabbatch{mbi}.spm.stats.factorial_design.dir = {out_dir};
        matlabbatch{mbi}.spm.stats.factorial_design.des.t1.scans = cellstr(allfiles);
        
        matlabbatch{mbi}.spm.stats.factorial_design.des.pt.gmsca  = 0;
        matlabbatch{mbi}.spm.stats.factorial_design.des.pt.ancova = 0;
        
        
        matlabbatch{mbi}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{mbi}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{mbi}.spm.stats.factorial_design.masking.im = 1;
        
        if bs == 1
            matlabbatch{mbi}.spm.stats.factorial_design.masking.em = {fullfile(path.templateDir, 'brainstem_mask.nii')};
        else
            matlabbatch{mbi}.spm.stats.factorial_design.masking.em = {fullfile(path.templateDir, vars.groupMaskID)};
        end
        matlabbatch{mbi}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{mbi}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{mbi}.spm.stats.factorial_design.globalm.glonorm = 1;
        mbi = mbi + 1;
        
        matlabbatch{mbi}.spm.stats.fmri_est.spmmat = {fullfile(out_dir,'SPM.mat')};
        matlabbatch{mbi}.spm.stats.fmri_est.method.Classical = 1;
        mbi = mbi + 1;
        
        matlabbatch{mbi}.spm.stats.con.spmmat = {fullfile(out_dir,'SPM.mat')};
        matlabbatch{mbi}.spm.stats.con.delete = 1;
        
        co = 1;
        matlabbatch{mbi}.spm.stats.con.consess{co}.tcon.name    = 'pos';
        matlabbatch{mbi}.spm.stats.con.consess{co}.tcon.convec  = [1];
        matlabbatch{mbi}.spm.stats.con.consess{co}.tcon.sessrep = 'none';
        co = co + 1; %increment by 1
        
        matlabbatch{mbi}.spm.stats.con.consess{co}.tcon.name    = 'neg';
        matlabbatch{mbi}.spm.stats.con.consess{co}.tcon.convec  = [-1];
        matlabbatch{mbi}.spm.stats.con.consess{co}.tcon.sessrep = 'none';
        co = co + 1; %increment by 1
        
        mbi = mbi + 1;
        

        copyfile(which(mfilename),fullfile(out_dir,spm_file(spm_file(which(mfilename),'filename'),'ext','bak'))); % copy this file into analysis directory as *.bak
        copyfile(which('get_study_specs'),fullfile(out_dir,spm_file(spm_file(which('get_study_specs'),'filename'),'ext','bak'))); % and copy this file as it contains everything necessary to repeat the analysis

        run_local(matlabbatch);
    end
end

if do_fact || do_fact_con
    addon   = [addon 'FACT'];
    
    out_dir = fullfile(path.secondlevelDir,[addon '_' anadirname '_' sm_str]);
    
    matlabbatch = [];
    matlabbatch{1}.spm.stats.factorial_design.dir = {out_dir};
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'GROUP';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'SUBJECT';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'CONDITION';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = analysis.fact_dept;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = analysis.fact_var;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;
    % if ana == 1  % FIR
    %     matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
    %     matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 0; %too many ...
    %     matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 1;
    % end
    
    %GROUP SUBJECT DAY CONDITION
    allfiles  = [];
    imat      = [];
    inc_su = 1;
    for gr = 1:size(groups,2) 
        for su = 1:size(groups{gr},2)
            for co = 1:numel(conds_done{groups_ind{gr}(su)}) % conds_done contains all relevant regressors
                name        = sprintf('sub-%02d',groups{gr}(su));
                a_dir       = fullfile(path.firstlevelDir, name, anadirname);
                s_string    = sprintf('s%s',sm_str);
                if analysis.skernel == 0;s_string = '';end
                if warp_s_con
                    swcon_file = fullfile(a_dir, sprintf(con_temp,simple_con_ind{groups_ind{gr}(su)}(co)));
                end
                if warp_beta % if we have both beta takes precedence
                    swcon_file = fullfile(a_dir, sprintf(beta_temp,simple_beta_ind{groups_ind{gr}(su)}(co)));
                end
                if analysis.use_vasa == 1
                    swcon_file = spm_file(swcon_file,'prefix',sprintf('%swv',s_string));
                else
                    swcon_file = spm_file(swcon_file,'prefix',sprintf('%sw',s_string));
                end
                allfiles    = strvcat(allfiles,swcon_file);
                mat_entry   = [1 gr inc_su conds_done{groups_ind{gr}(su)}(co)];
                imat        = [imat; mat_entry];
            end
            inc_su = inc_su + 1;
        end
    end
    
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans    = cellstr(allfiles);
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.imatrix  = imat;
    
    
    if n_groups>1
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.inter.fnums  = [1 3]; %group by condition
    else
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum  = [3]; % only conditions
    end
    
    matlabbatch{1}.spm.stats.factorial_design.cov                  = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov            = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none   = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im           = 1;
    
    if bs == 1
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {fullfile(path.templateDir, 'brainstem_mask.nii')};
    else
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {fullfile(path.templateDir, vars.groupMaskID)};
    end
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(out_dir,'SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
end

if do_fact
    mkdir(out_dir);

    copyfile(which(mfilename),fullfile(out_dir,spm_file(spm_file(which(mfilename),'filename'),'ext','bak'))); % copy this file into analysis directory as *.bak
    copyfile(which('get_study_specs'),fullfile(out_dir,spm_file(spm_file(which('get_study_specs'),'filename'),'ext','bak'))); % and copy this file as it contains everything necessary to repeat the analysis

    run_local(matlabbatch);
end

if do_fact_con
    matlabbatch = [];
    clear SPM; load(fullfile(out_dir,'SPM.mat')); %should exist by now
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(out_dir,'SPM.mat')};
    matlabbatch{1}.spm.stats.con.delete = 1;
    
    
    co = 1;
    matlabbatch{1}.spm.stats.con.consess{co}.fcon.name   = 'eff_of_int';
    Fc = spm_FcUtil('Set','F_iXO_Test','F','iX0',[],SPM.xX.X);
    matlabbatch{1}.spm.stats.con.consess{co}.fcon.convec = {Fc.c'};
    co = co + 1; %increment
    
    %do F-contrasts
    if ~do_ppi % hott fix bc we f contrasts are not yet implemented for ppi
        for fc=1:numel(f_con_names)
            matlabbatch{1}.spm.stats.con.consess{co}.fcon.name   = analysis.f_con_names{fc}; % we need the original 
            matlabbatch{1}.spm.stats.con.consess{co}.fcon.convec = analysis.f_con{fc};
            co = co + 1; %increment
        end
    end

    if do_ppi % hot fix for ppi bc for t-constrast definition below original contrasts are used. However, these are not the correct contrasts for the ppi
        analysis.t_con_names = t_con_names;
        analysis.t_con       = t_con;
    end

    for gw = 1:numel(group_ana_names)
        for tc = 1:numel(analysis.t_con_names)
            matlabbatch{1}.spm.stats.con.consess{co}.tcon.name    = [group_ana_names{gw} '_' analysis.t_con_names{tc}]; % we need the original 
            matlabbatch{1}.spm.stats.con.consess{co}.tcon.convec  = [kron(group_weights(gw,:),analysis.t_con(tc,:))]; % we need the original 
            matlabbatch{1}.spm.stats.con.consess{co}.tcon.sessrep = 'none';
            co = co + 1; %increment
        end
    end
    
    run_local(matlabbatch);
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

function Y = normit(X)
%normalise X to zero mean, std unity
Y = X - ones(size(X,1),1)*mean(X); %zero mean
Y = Y ./ (ones(size(Y,1),1)*std(Y)); %std 1
end

function out = not_excluded(ana_ex,sub_id,ses,run)
% check whether this run should be in the analysis
% not very elegant
out = 1;
for i_ae = 1:numel(ana_ex)
    if ana_ex(i_ae).sub == sub_id
        if ana_ex(i_ae).ses == ses
            for i_run = 1:numel(ana_ex(i_ae).run)
                if ana_ex(i_ae).run(i_run) == run
                    out = 0;
                end
            end
        end
    end
end
end


function [f_vec, n_reg] = create_f_vec(n_run, cond_use, n_base, n_nuis, n_o_run, is_ppi, n_ppi_conds, p_mod, ana)
    
    % construct the effects of int contrast across runs (complicated because some conds might be missing in some runs)
    % get some helpful vectors first
    i_reg = [];n_reg = [];
    for pp=1:numel(p_mod)
        n_reg(pp) = 1+numel(p_mod{pp});
        i_reg     = [i_reg repmat(pp,1,n_reg(pp))];
    end
    
    % Adjust for LSA if needed
    if ana == 3 % f-con for LSA needs adjustment
        max_cond = 0;
        for run = 1:n_run
            max_cond = max(max_cond, max(cond_use{run}));
        end
        i_reg = 1:max_cond;
        n_reg = ones(size(i_reg));
    end
    
    % actual contrast creation
    if ~is_ppi

        % Standard GLM F-contrast
        f_vec = [];
        for run = 1:n_run
            cvec{run} = zeros(sum(n_reg));
            for co = 1:numel(cond_use{run})
                cvec{run} = cvec{run} + diag(i_reg==cond_use{run}(co));
            end
            cvec{run} = kron(cvec{run}, eye(n_base));
            cvec{run}(:,find(~sum(cvec{run}))) = [];
            f_vec = [f_vec cvec{run} zeros(size(cvec{run},1), n_nuis)];
        end
        f_vec = [f_vec zeros(size(cvec{run},1), n_o_run)];
    else
        % PPI F-contrast
        n_ppi_regressors = n_ppi_conds*2 + 1;
        f_vec = [eye(n_ppi_regressors) zeros(n_ppi_regressors, n_nuis+n_o_run)];
    end

end