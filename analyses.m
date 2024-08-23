function analyses(all_sub_ids)
% does 1st level and second level analyses 

%% Prepare everything
% now read in stuff form get_study_specs
[path,vars,analysis]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

n_groups     = max(analysis.group_ind);

for g = 1:n_groups
    groups{g}       = intersect(analysis.all_subs(analysis.group_ind==g),all_sub_ids);
end

group_weights = analysis.group_weights;
group_names   = analysis.group_names;

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
    do_lss            = analysis.lss;  % also do an LSS analysis ?
else
    do_lss            = 0; % never for hrf and fir
end
concatenate     = analysis.concatenate;

do_model        = analysis.do_model;
bs              = analysis.bs;
do_est          = analysis.do_est ;
do_vasa         = analysis.do_vasa;
do_cons         = analysis.do_cons;
do_correct_vasa = analysis.do_correct_vasa;
do_warp         = analysis.do_warp;
do_smooth       = analysis.do_smooth;

% second level
do_fact         = analysis.do_fact;
do_fact_con     = analysis.do_fact_con;

do_one_t        = analysis.do_one_t;

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
    warp_con          = 1;
end

if ana == 2 % HRF
    warp_beta         = 0;
    warp_con          = 1;
end

if ana == 3 % HRF LS-A / LSS
    warp_beta         = 1;
    warp_con          = 0;
end

n_base            = analysis.n_base; % # of basis functions (FIR model)

name_ana          = analysis.name_ana;

% cond_names need to be EXACTLY the same as they appear in the .tsv files
cond_names        = analysis.cond_names;
n_cond            = numel(cond_names);
t_con             = analysis.t_con;
t_con_names       = analysis.t_con_names;

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
    if vars.nSess > 1
        warning('This routine will collapse all sessions into a single one !')
    end
    cnt        = 1;
    % for now all sessions are converted to runs
    
    %prepare exclusions
    
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            if (~isfield(analysis,'exclude')) || (not_excluded(analysis.exclude,sub_id,ses,run)) 
                epi = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
                epi = epi(1); % if there are brain and spinal, just take brain
                epifiles{cnt} = epi;
                tsv = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','events');
                tsvfiles{cnt} = tsv;
                cnt           = cnt + 1;
            end
        end
    end
    %n_run      = vars.nRuns*vars.nSess;
    n_run      = cnt - 1; % already incremented in loop
    n_o_run    = n_run; % keep original as # of runs get set to 1 for concatenate
    
    % get relevant dirs
    struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'type','T1w');
    struc_file = struc_file(1); % if there are more than 1
    u_rc1_file = spm_file(struc_file,'prefix','u_rc1');
    
    mean_dir = spm_file(epifiles{1},'path');
    
    template = []; template_wls = []; % everything will first go into a template struc and added to matlabbatch as needed
    
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
        template.spm.stats.fmri_spec.bases.fir.length = TR*n_base;
        template.spm.stats.fmri_spec.bases.fir.order  = n_base;
        
    elseif ana == 2 || ana == 3 % HRF (also for LSA LSS)
        template.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    end
    
    template_wls.spm.tools.rwls.fmri_rwls_spec.bases = template.spm.stats.fmri_spec.bases;
    
    a_dir    = fullfile(path.firstlevelDir,sprintf('sub-%02d',sub_id),anadirname);
    lss_ind = [];
    
    for run = 1:n_run
        func_dir = spm_file(epifiles{run},'path');
        
        % first get min and max onset to trim scans, not yet implemented
        min_onset = Inf;
        max_onset = -Inf;
        x = spm_load(char(spm_file(tsvfiles{run},'suffix', analysis.events))); % load onset *.tsv file; analysis.events allows different tsv files
        if isnumeric(x.trial_type) % the rare case that all conditions are numbers...
            x.trial_type = strtrim(cellstr(num2str(x.trial_type)));
        end
        for ii=1:numel(cond_names)
            t_ind = find(strcmp(x.trial_type,cond_names(ii)));
            if ~isempty(t_ind) || (concatenate == 1) % allow empty in some runs when concatenate, bc other runs have the condition (maybe check in the end whether any run has it)
                c_onset     = (x.onset(t_ind)')./TR + shift; % times in tsv files are in seconds according to BIDS
                min_onset   = min([min_onset c_onset]); % this can then bee used to trim scans
                max_onset   = max([max_onset c_onset]); % think about leaving about 8s worth of scans after last event
            end
        end
             
        % then get the nuisance vectors (movement, wm, csf etc)
        if bs == 1
            fm        = spm_file(epifiles{run},'prefix','rp_ba','ext','.txt');
        else
            fm        = spm_file(epifiles{run},'prefix','rp_a','ext','.txt');
        end
        movement  = normit(load(char(fm)));

       
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
            physio_noise_f = char(spm_file(epifiles{run},'prefix','physio_','ext','.mat'));
            if exist(physio_noise_f,'file')
                physio_noise = load(physio_noise_f);                
                all_nuis{run} = [all_nuis{run} normit(physio_noise.physio)];
            end
        end
        
        if strfind(noise_corr,'other')
            other_f = char(spm_file(epifiles{run},'prefix','other_','ext','.mat'));
            if exist(other_f,'file')
                other_cov = load(other_f);
                all_nuis{run} = [all_nuis{run} normit(other_cov.other)];
            end
        end
         
        
        if strfind(noise_corr,'wm')
            seg_noise_f = char(spm_file(epifiles{run},'prefix','noise_wm_csf_ra','ext','.mat'));
            if exist(seg_noise_f,'file')
                seg_noise = load(seg_noise_f);
                all_nuis{run} = [all_nuis{run} normit(seg_noise.segment(1).data)];
            end
        end
        
        
        if strfind(noise_corr,'csf')
            if bs == 1
                seg_noise_bs_f = char(spm_file(epifiles{run},'prefix','noise_csf_bra','ext','.mat'));
                if exist(seg_noise_bs_f,'file')
                    seg_noise_bs = load(seg_noise_bs_f);
                    all_nuis{run} = [all_nuis{run} normit(seg_noise_bs.segment(1).data)];
                end
            else
                seg_noise_f = char(spm_file(epifiles{run},'prefix','noise_wm_csf_ra','ext','.mat'));
                if exist(seg_noise_f,'file')
                    seg_noise = load(seg_noise_f);
                    all_nuis{run} = [all_nuis{run} normit(seg_noise.segment(2).data)];
                end
            end
        end
        
        if strfind(noise_corr,'roi')
            roi_noise_f = char(spm_file(epifiles{run},'prefix','noise_roi_ra','ext','.mat'));
            if exist(roi_noise_f,'file')
                roi_noise = load(roi_noise_f);                
                all_nuis{run} = [all_nuis{run} normit(roi_noise.roi(1).data)];
                all_nuis{run} = [all_nuis{run} normit(roi_noise.roi(2).data)];
            end
        end
        
        if strfind(noise_corr,'movnoise')
            mov_noise_f = char(spm_file(epifiles{run},'prefix','noise_mov_ra','ext','.mat'));
            if exist(mov_noise_f,'file')
                mov_noise = load(mov_noise_f);
                all_nuis{run} = [all_nuis{run} mov_noise.mov_reg];
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
        % create vector of number of scans for possible concatenation
        scan_vec(run) = size(scans,1);
        
        template.spm.stats.fmri_spec.sess(run).scans = cellstr(scans);
        template.spm.stats.fmri_spec.sess(run).multi = {''};
        template.spm.stats.fmri_spec.sess(run).multi_reg = {''};
        template.spm.stats.fmri_spec.sess(run).hpf       = analysis.hpf;
        
        x = spm_load(char(spm_file(tsvfiles{run},'suffix', analysis.events))); % load onset *.tsv file; analysis.events allows different tsv files
        if isnumeric(x.trial_type) % the rare case that all conditions are numbers...
            x.trial_type = strtrim(cellstr(num2str(x.trial_type)));
        end
        RES   = [];cnt = 1;cond_use{run} = []; %initialize some variables
        for ii=1:numel(cond_names)
            t_ind = find(strcmp(x.trial_type,cond_names(ii)));
            if ~isempty(t_ind) || (concatenate == 1) % allow empty in some runs when concatenate, bc other runs have the condition (maybe check in the end whether any run has it)
                cond_use{run} = [cond_use{run} ii];
                RES(cnt).name      = cond_names{ii};
                RES(cnt).onset     = (x.onset(t_ind)')./TR + shift; % times in tsv files are in seconds according to BIDS
                RES(cnt).duration  = (x.duration(t_ind)')./TR;
                if ana == 1 % FIR always has duration of 0
                    RES(cnt).duration  = zeros(size(x.duration(t_ind)'));
                end
                
                % now add parametric regressors
                RES(cnt).orth = 0; % no orthogonalization
                for pm = 1:numel(analysis.p_mod{ii})
                    eval(sprintf('t_pm = x.%s;',char(analysis.p_mod{ii}(pm)))); %load according to analysis.p_mod
                    if ~isnumeric(t_pm) % if a column has no NaNs spm_load returns a numeric vector, else strings
                        t_pm = str2num(strvcat(t_pm)); % convert to numbers
                    end
                    RES(cnt).pmod(pm).name  = char(analysis.p_mod{ii}(pm));
                    RES(cnt).pmod(pm).param = t_pm(t_ind); % and add
                    RES(cnt).pmod(pm).poly  = 1;
                end
                cnt = cnt + 1;
            end
        end
        % expand all conditions for LS-A type analysis
        ind  = 1;all = [];
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
        if ana == 3 % define for which events we want LS-S estimates (and specify their duration as spm_spm_lss will create different regressors accorsing to duration)
            lss_ind  = [1:numel(new_t.spm.stats.fmri_spec.sess(1).cond) ones(1,n_run + numel(z{1})).*NaN;...
                cat(1,new_t.spm.stats.fmri_spec.sess(1).cond.duration)' ones(1,n_run + numel(z{1})).*NaN];
        end
        n_run    = 1; % set # of runs to 1
        template_wls = new_t_wls;
        template     = new_t; % replace template
        z            = z(1); % adjust nuisance counter
    end
    
    if do_model
        mkdir(a_dir);
        copyfile(which(mfilename),a_dir); % copy this file into analysis directory
        copyfile(which('get_study_specs'),a_dir); % and copy this file as it contains everything necessary to repeat the analysis
        for ts = 1:numel(tsvfiles)
            copyfile(char(spm_file(tsvfiles{ts},'suffix', analysis.events)),a_dir); % and copy the respective *.tsv files that were used
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
    
    
    %% now prepare contrasts
    template = [];
    template.spm.stats.con.spmmat = {fullfile(a_dir,'SPM.mat')};
    template.spm.stats.con.delete = 1;
    fco = 1; % counter for f-contrasts
    template.spm.stats.con.consess{fco}.fcon.name  = 'eff_of_int';
    
    if (ana == 3) & (concatenate == 1)
        cond_use{1} = 1:numel(new_t.spm.stats.fmri_spec.sess(1).cond);
    end
    % construct the effects of int contrast across runs (complicated because some conds might be missing in some runs)
    % get some helpful vectors first
    i_reg = [];
    for pp=1:numel(analysis.p_mod)
        n_reg(pp) = 1+numel(analysis.p_mod{pp});
        i_reg     = [i_reg repmat(pp,1,n_reg(pp))];
    end
    
    if ana == 3 % f-con for LSA needs adjustment
        max_cond = 0;
        for run = 1:n_run
            max_cond = max(max_cond, max(cond_use{run}));
        end
        i_reg = 1:max_cond;
        n_reg = ones(size(i_reg));
    end
    
    
    f_vec = [];
    for run = 1:n_run
        cvec{run} = zeros(sum(n_reg));
        for co = 1:numel(cond_use{run})
            cvec{run} = cvec{run} + diag(i_reg==cond_use{run}(co)); %add in the eye(s)
        end
        cvec{run} = kron(cvec{run},eye(n_base)); % expand basis functions
        cvec{run}(:,find(~sum(cvec{run})))=[]; % remove conds that do not exist in this run
        f_vec = [f_vec cvec{run} zeros(size(cvec{run},1),n_nuis)]; % add zeros for nuisance vars
    end
    f_vec = [f_vec zeros(size(cvec{run},1),n_o_run)]; % add zeros for run constants
    template.spm.stats.con.consess{fco}.fcon.convec = {f_vec}; % complete eoi contrast matrix
    
    % now do the simple cons for 2nd level anova
    % simply go through the eoi (f_vec) line by line and correct by # of occurrence
    if ana == 1 || ana == 2
        simple_con_ind = [];simple_beta_ind = [];
        simple_con = zeros(n_cond*n_base,size(f_vec,2));
        ind = 1;
        for co = 1:n_cond
            for pm = 1:n_reg(co)
                if pm == 1
                    label = [cond_names{co} '_main' ]; % the event itself
                else
                    label = [cond_names{co} '_' analysis.p_mod{co}{pm-1}]; % the parametric regressors
                end
                for i_fir = 1:n_base % take care of FIR bins
                    simple_con(ind,:) = f_vec(ind,:)./sum(f_vec(ind,:)); % scale to 1
                    template.spm.stats.con.consess{ind+fco}.tcon.name    = [label '_' num2str(i_fir)];
                    all_t_con_names{ind}                                 = [label '_' num2str(i_fir)];
                    template.spm.stats.con.consess{ind+fco}.tcon.convec  = simple_con(ind,:);
                    template.spm.stats.con.consess{ind+fco}.tcon.sessrep = 'none';
                    simple_con_ind  = [simple_con_ind ind+fco];
                    simple_beta_ind = [simple_beta_ind ind];
                    ind = ind + 1;
                end
            end
        end
    end
    
    % additional t-constrasts as specified in get_study_specs, these will also be as one sample t-tests
    if ~isempty(t_con)
        add_con_ind = [];
        %for co = 1:size(t_con_names,2)
        for co = 1:numel(t_con_names)
            template.spm.stats.con.consess{ind+fco}.tcon.name    = [t_con_names{co}];
            all_t_con_names{ind}                                 = [t_con_names{co}];
            template.spm.stats.con.consess{ind+fco}.tcon.convec  = (simple_con'*t_con(co,:)')'; % nicely takes the weights from the simple t cons into account
            template.spm.stats.con.consess{ind+fco}.tcon.sessrep = 'none';
            add_con_ind  = [add_con_ind ind+fco];
            ind = ind + 1;
        end
    end
    
    if do_cons
        % get old con files to delete them *con* *spmT* *spmF* *ess* 
        old_con_files = [''];
        for cc = i:numel(con_templ)
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
    if warp_con
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
    
    if do_warp
        %using nlin coreg + DARTEL
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
    
    if n_procs > analysis.max_procs
        n_procs = analysis.max_procs;
    end
    
    if parallel == 1
        run_spm_parallel(matlabbatch, n_procs);
    else
        run_local(matlabbatch);
    end
end

%% second level analyses
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
            if warp_con
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
        
        copyfile(which(mfilename),out_dir);
        copyfile(which('get_study_specs'),out_dir);
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
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;
    if ana == 1  % FIR
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 0; %too many ...
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 1;
    end
    
    %GROUP SUBJECT DAY CONDITION
    allfiles  = [];
    imat      = [];
    inc_su = 1;
    for gr = 1:size(groups,2) % better numel(groups)
        for su = 1:size(groups{gr},2)
            for co = 1:size(simple_con_ind,2)
                name        = sprintf('sub-%02d',groups{gr}(su));
                a_dir       = fullfile(path.firstlevelDir, name, anadirname);
                s_string    = sprintf('s%s',sm_str);
                if analysis.skernel == 0;s_string = '';end
                if warp_beta
                    swcon_file = fullfile(a_dir, sprintf(beta_temp,simple_beta_ind(co)));
                end
                if warp_con
                    swcon_file = fullfile(a_dir, sprintf(con_temp,simple_con_ind(co)));
                end
                if analysis.use_vasa == 1
                    swcon_file = spm_file(swcon_file,'prefix',sprintf('%swv',s_string));
                else
                    swcon_file = spm_file(swcon_file,'prefix',sprintf('%sw',s_string));
                end
                allfiles    = strvcat(allfiles,swcon_file);
                mat_entry   = [1 gr inc_su co];
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
    copyfile(which(mfilename),out_dir);
    copyfile(which('get_study_specs'),out_dir);
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
    
    for gw = 1:numel(group_names)
        for tc = 1:numel(t_con_names)
            matlabbatch{1}.spm.stats.con.consess{co}.tcon.name    = [group_names{gw} '_' t_con_names{tc}];
            matlabbatch{1}.spm.stats.con.consess{co}.tcon.convec  = [kron(group_weights(gw,:),t_con(tc,:))];
            matlabbatch{1}.spm.stats.con.consess{co}.tcon.sessrep = 'none';
            co = co + 1; %increment
        end
    end
    
    run_local(matlabbatch);
end


function run_local(matlabbatch)

if size(matlabbatch,2) > 1 % we have multiple subjects
    matlabbatch = matlabbatch(:);
end
spm_jobman('initcfg');
spm('defaults', 'FMRI');
spm_jobman('run',matlabbatch);

function Y = normit(X)
%normalise X to zero mean, std unity
Y = X - ones(size(X,1),1)*mean(X); %zero mean
Y = Y ./ (ones(size(Y,1),1)*std(Y)); %std 1


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
