function analyses(all_sub_ids)
% example for a 1st löevel and secopnd löevel anaylsis function

% resolve paths
[path,vars,analysis]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

n_groups     = max(analysis.group_ind);

for g = 1:n_groups
    groups{g}       = intersect(analysis.all_subs(analysis.group_ind==g),all_sub_ids);
end

group_weights = analysis.group_weights;
group_names   = analysis.group_names;


%% house keeping
con_temp          = 'con_%04.4d.nii';
beta_temp         = 'beta_%04.4d.nii';

TR                = vars.sliceTiming.tr;


parallel          = analysis.parallel;

%% define what needs to be done
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
concatenate       = analysis.concatenate;
% TO DO move the following the get_study_specs


do_model        = analysis.do_model;
wls             = analysis.wls;
bs              = analysis.bs;
do_est          = analysis.do_est ;
do_vasa         = analysis.do_vasa;
do_cons         = analysis.do_cons;
do_correct_vasa = analysis.do_correct_vasa;
do_warp         = analysis.do_warp;
do_smooth       = analysis.do_smooth;

%second level
do_fact         = analysis.do_fact;
do_fact_con     = analysis.do_fact_con;

do_one_t        = analysis.do_one_t;
% TO DO unil here

% other things
noise_corr        = analysis.noise_corr; %the whole lot
cvi               = analysis.cvi;
shift             = analysis.shift; %in TRs


spm_path          = fileparts(which('spm')); %get spm path
mat_name          = which(mfilename);
[~,mat_name,~]    = fileparts(mat_name);


if ana == 1 % Offset FIR
    warp_beta         = 0;
    warp_con          = 1;
end

if ana == 2 % Offset HRF
    warp_beta         = 0;
    warp_con          = 1;
end

if ana == 3 % Offset HRF LS-A / LSS
    warp_beta         = 1;
    warp_con          = 0;
end

%% now read in stuff form get_study_specs
n_base            = analysis.n_base; %

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

if wls
    anadirname = ['WLS_' anadirname];
end

if concatenate
    anadirname = ['Conc_' anadirname];
end


%% main loop
matlabbatch = [];
for sub = 1:n_subs
    mbi        = 0;
    sub_id     = all_sub_ids(sub);
    if vars.nSess > 1
        warning('This routine will collapse all sessions into a single one !')
    end
    cnt        = 1;
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            path.epifiles{cnt} =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            cnt                =  cnt + 1;
        end
    end
    n_run      = vars.nRuns*vars.nSess;
    n_o_run    = n_run; % keep original #runs for concatenate
    % get relevant dirs
    
    struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'type','T1w');
    u_rc1_file = spm_file(struc_file,'prefix','u_rc1');
    
    func_dir = spm_file(path.epifiles{run},'path');
    
    
    
    template = []; template_wls = []; %everything will go in a template struc and added to matlabbatch as needed
    
    template.spm.stats.fmri_spec.timing.units   = 'scans';
    template_wls.spm.tools.rwls.fmri_rwls_spec.timing.units = template.spm.stats.fmri_spec.timing.units;
    
    template.spm.stats.fmri_spec.timing.RT      = TR;
    template_wls.spm.tools.rwls.fmri_rwls_spec.timing.RT = template.spm.stats.fmri_spec.timing.RT;
    
    template.spm.stats.fmri_spec.timing.fmri_t  = 16;
    template_wls.spm.tools.rwls.fmri_rwls_spec.timing.fmri_t = template.spm.stats.fmri_spec.timing.fmri_t;
    
    template.spm.stats.fmri_spec.timing.fmri_t0 = 8;
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
        template.spm.stats.fmri_spec.mask           = cellstr(fullfile(func_dir,'bins3cbrainstem_mask.nii'));
        template_wls.spm.tools.rwls.fmri_rwls_spec.mask = template.spm.stats.fmri_spec.mask;
    else
        template.spm.stats.fmri_spec.mask           = cellstr(fullfile(func_dir,'s3rbrain_mask.nii'));
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
    
    for run = 1:n_run
        % first get the nuisance vectors (movement, wm, csf etc)
        if bs == 1
            fm        = spm_file(path.epifiles{run},'prefix','rp_ba','ext','.txt');
        else
            fm        = spm_file(path.epifiles{run},'prefix','rp_a','ext','.txt');
        end
        movement  = normit(load(char(fm)));
        physio_noise_f = char(spm_file(path.epifiles{run},'prefix','physio_','ext','.mat'));
        if exist(physio_noise_f,'file')
            physio_noise = load(physio_noise_f);
        end
        seg_noise_bs_f = char(spm_file(path.epifiles{run},'prefix','noise_csf_bra','ext','.mat'));
        if exist(seg_noise_bs_f,'file')
            seg_noise_bs = load(seg_noise_bs_f);
        end
        seg_noise_f = char(spm_file(path.epifiles{run},'prefix','noise_wm_csf_ra','ext','.mat'));
        if exist(seg_noise_f,'file')
            seg_noise = load(seg_noise_f);
        end
        roi_noise_f = char(spm_file(path.epifiles{run},'prefix','noise_roi_ra','ext','.mat'));
        if exist(roi_noise_f,'file')
            roi_noise = load(roi_noise_f);
        end
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
        if strfind(noise_corr,'physio')
            all_nuis{run} = [all_nuis{run} normit(physio_noise.physio)];
        end
        if strfind(noise_corr,'wm')
            all_nuis{run} = [all_nuis{run} normit(seg_noise.segment(1).data)];
        end
        if strfind(noise_corr,'csf')
            if bs == 1
                all_nuis{run} = [all_nuis{run} normit(seg_noise_bs.segment(1).data)];
            else
                all_nuis{run} = [all_nuis{run} normit(seg_noise.segment(2).data)];
            end
        end
        if strfind(noise_corr,'roi')
            all_nuis{run} = [all_nuis{run} normit(roi_noise.roi(1).data)];
            all_nuis{run} = [all_nuis{run} normit(roi_noise.roi(2).data)];
        end
        
        n_nuis         = size(all_nuis{run},2);
        z{run}         = zeros(1,n_nuis); % handy for later contrast def
        
        if bs == 1
            scans = spm_select('ExtFPlist', func_dir,spm_file(spm_file(path.epifiles{run},'filename'),'prefix','^bra'),Inf);
        else
            scans = spm_select('ExtFPlist', func_dir,spm_file(spm_file(path.epifiles{run},'filename'),'prefix','^ra'),Inf);
        end
        % create vector of scans for concatenation
        scan_vec(run) = size(scans,1);
        
        template.spm.stats.fmri_spec.sess(run).scans = cellstr(scans);
        template.spm.stats.fmri_spec.sess(run).multi = {''};
        template.spm.stats.fmri_spec.sess(run).multi_reg = {''};
        template.spm.stats.fmri_spec.sess(run).hpf       = 180;
        
        % we use tsv files
        tsv_file   = spm_file(path.epifiles{run},'ext','.tsv');
        x = spm_load(char(tsv_file));
        % items = unique(x.trial_type,'sorted'); % with this we could get all the conditions from TSV but it is better to specify the events to use
        % a priori in cond_names this also allows to specifiy the order of the conditions
        
        items = cond_names;
        RES   = [];cnt = 1;
        cond_use{run} = [];
        for ii=1:numel(items)
            t_ind = find(strcmp(x.trial_type,items(ii)));
            if ~isempty(t_ind)
                cond_use{run} = [cond_use{run} ii];
                RES(cnt).name      = items{ii};
                RES(cnt).onset     = x.onset(t_ind)' + shift;
                RES(cnt).duration  = x.duration(t_ind)';
                if ana == 1 % FIR always has duration of 0
                    RES(cnt).duration  = zeros(size(x.duration(t_ind)'));
                end
                cnt = cnt + 1;
            end
        end
        % expand all conditions for LS-A type analysis
        ind  = 1;all = [];
        for r = 1:numel(RES)
            new = [];
            new(:,1) = RES(r).onset';
            new(:,2) = repmat(r,size(new(:,1),1),1);
            new(:,3) = RES(r).duration';
            all = [all ; new];
        end
        [~,i] = sort(all(:,1)); % sort all events by scans (or comment these 2 lines)
        all   = all(i,:);
        for j=1:size(all,1)
            RES_LSA(j).name       = RES(all(j,2)).name;
            RES_LSA(j).onset      = all(j,1);
            RES_LSA(j).duration   = all(j,3);
        end
        
        % now put all in batch structure
        if ana == 1 || ana == 2 % FIR or HRF
            template.spm.stats.fmri_spec.sess(run).cond = RES;
        elseif ana == 3 % HRF LS-A
            template.spm.stats.fmri_spec.sess(run).cond = RES_LSA;
            cond_use{run}  = [1:size(all,1)]; % recalc cond_use
            if run>1
                cond_use{run} = cond_use{run} + max(cond_use{run-1})
            end
        end
        
        for nuis = 1:n_nuis
            template.spm.stats.fmri_spec.sess(run).regress(nuis) = struct('name', cellstr(num2str(nuis)), 'val', all_nuis{run}(:,nuis));
        end
    end
    
    template_wls.spm.tools.rwls.fmri_rwls_spec.sess = template.spm.stats.fmri_spec.sess;
    
    if concatenate
        % now take template struc and create a single run ...
        c_scan_vec = cumsum(scan_vec);
        if wls
            new_t = template_wls;
        else
            new_t = template;
        end
        
        new_t.spm.stats.fmri_spec.sess(2:end) = []; % get rid of original runs 2..n
        for run = 2:n_run
            % first do nuisance variables (eg motion)
            new_t.spm.stats.fmri_spec.sess(1).scans = [new_t.spm.stats.fmri_spec.sess(1).scans;template.spm.stats.fmri_spec.sess(run).scans];
            for re = 1:numel(new_t.spm.stats.fmri_spec.sess(1).regress)
                if numel(new_t.spm.stats.fmri_spec.sess(1).regress) == numel(template.spm.stats.fmri_spec.sess(run).regress)
                    new_t.spm.stats.fmri_spec.sess(1).regress(re).val = [new_t.spm.stats.fmri_spec.sess(1).regress(re).val; template.spm.stats.fmri_spec.sess(run).regress(re).val];
                else
                    error('number of nuisance variables needs to be identical');
                end
            end
            % next do conditions
            l_c = numel(template.spm.stats.fmri_spec.sess(run).cond);
            for c=1:l_c
                if ana == 3 %LS A
                    new_t.spm.stats.fmri_spec.sess(1).cond(end+1).name   = template.spm.stats.fmri_spec.sess(run).cond(c).name; %end is now end+1 ... be careful
                    new_t.spm.stats.fmri_spec.sess(1).cond(end).onset    = template.spm.stats.fmri_spec.sess(run).cond(c).onset + c_scan_vec(run-1);
                    new_t.spm.stats.fmri_spec.sess(1).cond(end).duration = template.spm.stats.fmri_spec.sess(run).cond(c).duration;
                else
                    if strcmp(template.spm.stats.fmri_spec.sess(1).cond(c).name,[template.spm.stats.fmri_spec.sess(run).cond(c).name]) % OK same name
                        new_t.spm.stats.fmri_spec.sess(1).cond(c).onset     = [new_t.spm.stats.fmri_spec.sess(1).cond(c).onset template.spm.stats.fmri_spec.sess(run).cond(c).onset + c_scan_vec(run-1)];
                        new_t.spm.stats.fmri_spec.sess(1).cond(c).duration  = [new_t.spm.stats.fmri_spec.sess(1).cond(c).duration template.spm.stats.fmri_spec.sess(run).cond(c).duration];
                    else
                        error('Conditions or order do not match betweeen runs')
                    end
                end
            end
        end
        lss_ind  = [1:numel(new_t.spm.stats.fmri_spec.sess(1).cond) ones(1,n_run + numel(z{1})).*NaN];
        %cond_use =
        %[1:n_base*numel(new_t.spm.stats.fmri_spec.sess(1).cond)]; --> not needed anymore
        n_run    = 1; % set # of runs to 1
        if wls
            template_wls = new_t;
        else
            template = new_t;
        end
        z        = z(1); % adjust nuisance counter
    end
    
    if do_model
        
        mbi = mbi + 1;
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = sprintf('--------------\ndoing Subject %d\n--------------\n',sub_id);
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'fprintf';
        mbi = mbi + 1;
        
        mkdir(a_dir);
        copyfile(which(mfilename),a_dir);
        copyfile(which('get_study_specs'),a_dir);
        
        template.spm.stats.fmri_spec.dir = {[a_dir]};
        template_wls.spm.tools.rwls.fmri_rwls_spec.dir =  template.spm.stats.fmri_spec.dir;
        
        if wls
            matlabbatch{mbi,sub} = template_wls;
        else
            matlabbatch{mbi,sub} = template;
        end
        %here we decide whether we do WLS or OLS
        
        
        if concatenate
            mbi = mbi + 1;
            matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = fullfile(a_dir,'SPM.mat');
            matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.evaluated = scan_vec;
            matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
            matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'spm_fmri_concatenate';
        end
    end
    
    if do_est
        mbi = mbi + 1;
        if wls
            matlabbatch{mbi,sub}.spm.tools.rwls.fmri_rwls_est.spmmat = {fullfile(a_dir,'SPM.mat')};
            matlabbatch{mbi,sub}.spm.tools.rwls.fmri_rwls_est.method.Classical = 1;
        else
            matlabbatch{mbi,sub}.spm.stats.fmri_est.spmmat           = {fullfile(a_dir,'SPM.mat')};
            matlabbatch{mbi,sub}.spm.stats.fmri_est.method.Classical = 1;
        end
    end
    
    if do_lss
        mbi = mbi + 1;
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = fullfile(a_dir,'SPM.mat');
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{2}.evaluated = lss_ind;
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'spm_spm_lss';
    end
    
    if do_vasa
        mbi = mbi + 1;
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = fullfile(a_dir,'SPM.mat');
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.outputs = {};
        matlabbatch{mbi,sub}.cfg_basicio.run_ops.call_matlab.fun = 'cb_vasa';
    end
    
    
    %%template for contrasts
    template = [];
    template.spm.stats.con.spmmat = {fullfile(a_dir,'SPM.mat')};
    template.spm.stats.con.delete = 1;
    fco = 0;
    fco = fco + 1; %counter for f-contrasts
    template.spm.stats.con.consess{fco}.fcon.name   = 'eff_of_int';
    
    
    
    if ana == 1 || ana == 2
        total_cond = items;
    end
    if ana == 3
        if concatenate == 1
            cond_use{1} = 1:numel(new_t.spm.stats.fmri_spec.sess(1).cond);
        end
        total_cond = 1:max(cond_use{n_run});
    end
    f_vec = []; % here we construct the effects of int contrast (complicated because not all conds could be present in all runs)
    for run = 1:n_run
        cvec{run} = [];
        for co = 1:numel(total_cond)
            d = cond_use{run} == co;
            cvec{run} = [cvec{run};kron(d,eye(n_base))];
        end
        f_vec = [f_vec cvec{run} zeros(size(cvec{run},1),n_nuis)];%add zeros for nuisance vars
    end
    f_vec = [f_vec zeros(size(cvec{run},1),n_o_run)]; % add zeros for run constants
    template.spm.stats.con.consess{fco}.fcon.convec = {f_vec};
    
    % now do the simple cons for 2nd level anova
    % simply go through the eoi (f_vec) line by line and correct by # of occurrence
    if ana == 1 || ana == 2
        co_i  = 0;
        simple_con_ind = [];simple_beta_ind = [];
        simple_con = zeros(n_cond*n_base,size(f_vec,2));
        for co = 1:n_cond
            for i_fir = 1:n_base
                ind = (co-1)*n_base + i_fir;
                simple_con(ind,:) = f_vec(ind,:)./sum(f_vec(ind,:)); % scale to 1
                co_i = co_i + 1;
                template.spm.stats.con.consess{co_i+fco}.tcon.name    = [cond_names{co} '_' num2str(i_fir)];
                all_t_con_names{co_i}                                 = [cond_names{co} '_' num2str(i_fir)];
                template.spm.stats.con.consess{co_i+fco}.tcon.convec  = simple_con(ind,:);
                template.spm.stats.con.consess{co_i+fco}.tcon.sessrep = 'none';
                simple_con_ind  = [simple_con_ind co_i+fco];
                simple_beta_ind = [simple_beta_ind co_i];
            end
        end
    end
    
    % additional t-constrasts e.g. for one sample t test
    if ~isempty(t_con) % additional t-constrasts e.g. for one sample t test
        add_con_ind = [];
        for co = 1:size(t_con_names,2)
            co_i = co_i + 1;
            template.spm.stats.con.consess{co_i+fco}.tcon.name    = [t_con_names{co}];
            all_t_con_names{co_i}                                 = [t_con_names{co}];
            template.spm.stats.con.consess{co_i+fco}.tcon.convec  = (simple_con'*t_con(co,:)')'; % nicely takes the weights from the simple t cons into account
            template.spm.stats.con.consess{co_i+fco}.tcon.sessrep = 'none';
            add_con_ind  = [add_con_ind co_i+fco];
        end
    end
    
    if do_cons
        mbi = mbi + 1;
        matlabbatch{mbi,sub} = template; % now add constrasts to batch
    end
    
    % prepare_warp
    template    = [];
    warp_files  = '';
    if warp_con
        for co = 1:size(all_t_con_names,2)
            warp_files(co,:) = [a_dir filesep sprintf(con_temp,co+fco)]; % need take f cons into account
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
            mbi = mbi + 1;
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
        end
    end
    
    if do_warp
        %using nlin coreg + DARTEL
        mbi = mbi + 1;
        matlabbatch{mbi,sub}.spm.util.defs.comp{1}.def = fullfile(func_dir,'y_epi_2_template.nii');
        matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.fnames = cellstr(spm_file(warp_files, 'prefix','v'));
        %matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.fnames = cellstr(spm_file(warp_files, 'prefix',''));
        matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
        matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.interp = 4;
        matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{mbi,sub}.spm.util.defs.out{1}.pull.prefix = 'w';
    end
    
    if do_smooth
        mbi = mbi + 1;
        matlabbatch{mbi,sub}.spm.spatial.smooth.data   = cellstr(spm_file(warp_files, 'prefix','wv'));
        %matlabbatch{mbi,sub}.spm.spatial.smooth.data   = cellstr(spm_file(warp_files, 'prefix','w'));
        matlabbatch{mbi,sub}.spm.spatial.smooth.fwhm   = skernel;
        matlabbatch{mbi,sub}.spm.spatial.smooth.prefix = ['s' sm_str];
    end
    
    % save('temp.mat','matlabbatch');
    % this can be used for debugging, ie simply open temp.mat in batch editor GUI
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

%% second level analysis


if do_one_t
    if size(groups,2)> 1
        warning('One sample t-test only works for a single group, will only analyse group 1');
    end
    addon   = 'ONE_T';
    % first let's get rid of negative contrasts
    cc    = corrcoef(t_con');
    [~,c] = find(triu(cc)<-.9999999999); % all in c are redundant
    pruned_t_con_names = t_con_names;
    pruned_t_con_names(c) = [];
    pruned_add_con_ind = add_con_ind;
    pruned_add_con_ind(c) = [];
    for co = 1:size(pruned_t_con_names,2)
        out_dir = fullfile(path.secondlevelDir,[addon '_' anadirname '_' sm_str '_' pruned_t_con_names{co}]);
        out_dir = strrep(out_dir,'>','_bt_');
        out_dir = strrep(out_dir,'<','_st_');
        mkdir(out_dir);
        
        matlabbatch = [];
        allfiles    = [];
        
        %% --------------------- MODEL SPECIFICATION --------------------- %%
        
        for su = 1:numel(groups{1})
            name        = sprintf('sub-%02d',groups{1}(su));
            a_dir       = fullfile(path.firstlevelDir, name, anadirname);
            s_string    = sprintf('s%s',sm_str);
            if analysis.skernel == 0;s_string = '';end
            if warp_con
                swcon_file = fullfile(a_dir, sprintf(con_temp,pruned_add_con_ind(co)));
            end
            swcon_file = spm_file(swcon_file,'prefix',sprintf('%swv',s_string));
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
        %% --------------------- MODEL ESTIMATION --------------------- %%
        
        matlabbatch{mbi}.spm.stats.fmri_est.spmmat = {[out_dir '\SPM.mat']};
        matlabbatch{mbi}.spm.stats.fmri_est.method.Classical = 1;
        mbi = mbi + 1;
        % --------------------- CONTRASTS --------------------- %%
        matlabbatch{mbi}.spm.stats.con.spmmat = {[out_dir '\SPM.mat']};
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
    addon   = 'FACT';
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
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
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
                swcon_file = spm_file(swcon_file,'prefix',sprintf('%swv',s_string));
                %swcon_file = spm_file(swcon_file,'prefix',sprintf('%sw',s_string));
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
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum  = [3]; % condition
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
    
    
    %% --------------------- MODEL ESTIMATION --------------------- %%
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {[out_dir '\SPM.mat']};
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
    clear SPM; load([out_dir '\SPM.mat']); %should exist by now
    matlabbatch{1}.spm.stats.con.spmmat = {[out_dir '\SPM.mat']};
    matlabbatch{1}.spm.stats.con.delete = 1;
    
    
    co = 1;
    matlabbatch{1}.spm.stats.con.consess{co}.fcon.name   = 'eff_of_int';
    Fc = spm_FcUtil('Set','F_iXO_Test','F','iX0',[],SPM.xX.X);
    matlabbatch{1}.spm.stats.con.consess{co}.fcon.convec = {Fc.c'};
    co = co + 1; %increment by 1
    
    for gw = 1:numel(group_names)
        for tc = 1:numel(t_con_names)
            matlabbatch{1}.spm.stats.con.consess{co}.tcon.name    = [group_names{gw} '_' t_con_names{tc}];
            matlabbatch{1}.spm.stats.con.consess{co}.tcon.convec  = [kron(group_weights(gw,:),t_con(tc,:))];
            matlabbatch{1}.spm.stats.con.consess{co}.tcon.sessrep = 'none';
            co = co + 1; %increment by 1
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
%save('test.mat','matlabbatch');
spm_jobman('run',matlabbatch);



function Y = normit(X)
%normalise X to zero mean, std unity
Y = X - ones(size(X,1),1)*mean(X); %zero mean
Y = Y ./ (ones(size(Y,1),1)*std(Y)); %std 1