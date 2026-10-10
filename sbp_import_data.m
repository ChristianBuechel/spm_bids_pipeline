function sbp_import_data(all_sub_ids,dry_run)

% function sbp_import_data(all_sub_ids,dry_run)
% Retrieve images from DICOM server using dicq and SFTP
% it also DICOM converts them and makes 4D EPI
% everything is arranged according to light weight BIDS
% who and what to import is specified in get_study_specs or in
% participants.tsv (preferred)
%
% Flow (several subjects are done in parallel, as in the other modules):
%  0) query dicq and PLAN which series to import (plan_series, applies
%     import.data: seq, cond, task) -> skipped series are listed with the reason
%  1) BATCH 1 (parallel, at most import.max_network_procs subjects at a time,
%     default 4, to be nice to the server): SFTP download and DICOM conversion
%  2) serial step: read DICOM headers (mag/phase), assign BIDS run numbers
%     (assign_runs; mag+phase of one acquisition share a run number and get
%     part-mag / part-phase, runs are counted per modality and task), split brain
%     and spinal EPIs (vars.acq_type, classify_geometry), build names (bids_name,
%     acq_tag), write the json files
%  3) BATCH 2 (parallel, vars.max_procs): 3D to 4D EPIs, anat and fmap files moved
%     to their BIDS name, gzip copies in path.baseDir, temp files deleted
% The parallel runs need the Parallel Computing Toolbox (the second step needs the
% first one to be finished), without it everything is done one subject after the
% other.
% dry_run = true stops after the plan table (no download of the images, only the
% first DICOM(s) of each functional series are fetched to read their headers).
%
% optional fields in import.data(n): task (task label, default vars.task),
% e.g. for short reference runs that should not be preprocessed as fMRI
% vars.acq_type: {'brain'} (default), {'spinal'} or {'spinal','brain'}
% import.max_network_procs: how many subjects are downloaded at the same time (4)
% NOTE: the password is part of the download jobs (it is also saved in the batch
% file if the batch is run via run_spm_multiple). On the server none is needed.
%
% HOW TO EXCLUDE SERIES (e.g. aborted runs, repeated T1s)
% A series is only imported if it matches a seq AND its cond (e.g. n > 300) in
% import.data. Series that match a seq but fail cond are listed as SKIPPED.
% To leave out a series that would otherwise be imported, list its series
% number in () directly behind the PRISMA number in the scan_id column of
% participants.tsv (NO space before the bracket, several series separated by
% commas, several sessions by spaces):
%   participant_id  sex  age  scan_id
%   sub-09          F    20   25879 25878(13)       --> series 13 of session 2 is omitted
%   sub-12          M    31   28000(7,12) 28001     --> series 7 and 12 of session 1 omitted
% Run sbp_import_data(ids,true) first to see the series numbers and what would
% be imported. Excluded series are listed in the protocol file (see below).
% Run numbers are counted AFTER the exclusion, i.e. excluding the first T1 makes
% the second one run-01.
%
% The transcript of the SFTP download is saved in <scanner>_<number>_ftp.log (not shown).
% Besides the SFTP batch file (<scanner>_<number>_ftp.txt) a short protocol
% <scanner>_<number>_protocol.txt is written into the subject folder: which
% series became which BIDS files, which series were not imported (and why), and
% how long the two batches took. A session that failed is marked ERROR there and
% in the summary at the end (its temp folders are kept for inspection).

if nargin<1
    error('You need to specify which subjects you want to import, e.g. import_data_ses([2 4 5])');
end
if nargin<2
    dry_run = false;
end
specs  = get_study_specs;
path   = specs.path;
vars   = specs.vars;
import = specs.import;

if isfield(import,'scanner')
    scanner = import.scanner; %if this is TRIO we download TRIO images
else
    scanner = 'PRISMA';
end

if isfield(vars,'acq_type')
    acq_type = vars.acq_type;
else
    acq_type = {'brain'};
end
acq_type = cellstr(acq_type);

check_task_labels(vars,import,path); % stop early with a clear message

import = read_import_specs(import);

machina = get_machine;
pw = [];
if any(strcmp(machina,{'win_loc','lin_loc'})) % only for those we need a pw
    pw = passcode; % assumes that on the server you are already logged in
end
conn = struct('user',import.user,'server',import.server);

if isfield(import,'max_network_procs')
    max_net = import.max_network_procs;
else
    max_net = 4;
end

if ~dry_run
    mkdir(path.preprocDir); % create coarse directory structure
    spm_save(fullfile(path.preprocDir,'.bidsignore.'),{sprintf('**/%s_*_ftp.txt',scanner),sprintf('**/%s_*_ftp.log',scanner),sprintf('**/%s_*_protocol.txt',scanner)}); % creates .bidsignore file
    descr = struct('Name',vars.task,'BIDSVersion','1.0.2');
    spm_jsonwrite(fullfile(path.preprocDir,'dataset_description.json'),descr, struct('indent','  ')); % create description file
    spm_jobman('initcfg');
end

%% Stage 0: plan all sessions (serial)
S = struct('col',{},'sub',{},'ses',{},'sub_str',{},'ses_str',{},'exam_id',{},'planned',{},'skipped',{},'ftp_file',{},'ok',{},'msg',{},'out_paths',{},'out_txt',{});
n_col = 0; % one column of the batch per subject
for s = 1:numel(all_sub_ids) % across volunteers
    c_sub = all_sub_ids(s);
    ind_c = find(import.prisma_no==c_sub);
    if isempty(ind_c)
        error('Subject %d is not listed. Abort\n',c_sub);
    end
    sub_str  = sprintf('sub-%2.2d',c_sub);
    this_col = 0;
    for ses = 1:numel(import.prisma{ind_c})
        ses_str = sprintf('ses-%2.2d',ses);
        exam_id = sprintf('%s_%d',scanner,import.prisma{ind_c}{ses});
        fprintf('Doing sub #%d session %d (%s) \n',c_sub,ses,exam_id);

        [status, result] = remote_cmd(machina,conn,pw,sprintf('netapp dicq -f --series --exam=%s',exam_id));
        if status ~= 0
            error('cannot retrieve output from dicq');
        end
        series = parse_dicq(result);
        [planned,skipped] = plan_series(series,import.data,import.prisma_skip{ind_c}{ses});
        if isempty(planned)
            warning('Nothing to import for sub %d session %d',c_sub,ses);
            continue
        end
        for k = 1:numel(planned)
            planned(k).dir = fullfile(path.preprocDir,sub_str,ses_str,planned(k).mod,sprintf('%s_%d_s%d',planned(k).type,ses,planned(k).ser)); % temp dir for DICOM files
        end

        if dry_run
            % download only the first few DICOM files of each func series, to read mag/phase
            % and the brain/spinal geometry from the headers (nothing else is downloaded)
            tmp_root = tempname;
            is_func  = find(strcmp({planned.mod},'func'));
            dest     = cell(size(is_func));
            for i = 1:numel(is_func)
                dest{i} = fullfile(tmp_root,sprintf('s%d',planned(is_func(i)).ser));
                mkdir(dest{i});
            end
            if numel(acq_type) == 1, n_peek = 1; else, n_peek = 6; end % combined acquisition: both geometries needed
            fetch_first_files(machina,conn,pw,{planned(is_func).path},dest,n_peek);
            for i = 1:numel(is_func)
                k = is_func(i);
                all_dicom = spm_select('FPlist',dest{i},'^MR.*');
                try
                    planned(k).subtype = read_series_headers(all_dicom,acq_type);
                catch err
                    warning('series %d: could not read headers (%s)',planned(k).ser,err.message);
                end
            end
            if isfolder(tmp_root), rmdir(tmp_root,'s'); end
            [~,o]   = sort([planned.ser]);
            planned = planned(o);
            if any(strcmp({planned.mod},'func') & cellfun(@isempty,{planned.subtype}))
                warning('DICOM headers are missing for some functional series: the run numbers and parts below are NOT reliable');
            end
            planned = assign_runs(planned);
            print_import_plan(planned,skipped,vars.task,acq_type);
            continue
        end

        % SFTP batch file: one temp dir per series
        ftp_cmd = cell(1,0);
        for k = 1:numel(planned)
            if isfolder(planned(k).dir) % left over from an aborted import: start from scratch
                rmdir(planned(k).dir,'s');
            end
            mkdir(planned(k).dir);
            ftp_cmd{end+1} = ['lcd ' planned(k).dir]; %#ok<AGROW>
            ftp_cmd{end+1} = ['cd ' planned(k).path]; %#ok<AGROW>
            ftp_cmd{end+1} = 'mget MR.*'; %#ok<AGROW> % get all dicom files in that dir
        end
        ftp_file = fullfile(path.preprocDir,sub_str,sprintf('%s_ftp.txt',exam_id));
        spm_save(ftp_file,ftp_cmd);
        fprintf('  %d series to download, %d not imported\n',numel(planned),numel(skipped));

        if this_col == 0
            n_col    = n_col + 1;
            this_col = n_col;
        end
        S(end+1) = struct('col',this_col,'sub',c_sub,'ses',ses,'sub_str',sub_str,'ses_str',ses_str,'exam_id',exam_id, ...
            'planned',{planned},'skipped',{skipped},'ftp_file',ftp_file,'ok',true,'msg','','out_paths',{{}},'out_txt',{{}}); %#ok<AGROW>
    end
end
if dry_run || isempty(S)
    return
end

use_par  = true;
why_seq  = '';
if ~(isfield(vars,'parallel') && vars.parallel == 1)
    use_par = false; why_seq = 'vars.parallel is not 1 in get_study_specs';
elseif ~license('test','Distrib_Computing_Toolbox')
    use_par = false; why_seq = 'no Parallel Computing Toolbox licence (the second step needs the first to be finished)';
elseif n_col == 1
    use_par = false; why_seq = 'only one subject';
end
n_net = min(n_col,max_net);
n_cpu = min(n_col,vars.max_procs);
if use_par
    fprintf('\n%d subject(s): batch 1 (download, conversion) with %d parallel worker(s) [import.max_network_procs = %d], batch 2 with %d\n',n_col,n_net,max_net,n_cpu);
else
    fprintf('\n%d subject(s): downloading and importing one subject after the other (%s)\n',n_col,why_seq);
end

%% Batch 1: download and DICOM conversion
b1 = cell(1,n_col);
g1 = zeros(1,n_col);
for i = 1:numel(S)
    c = S(i).col;
    dl.machina  = machina;
    dl.conn     = conn;
    dl.pw       = pw;
    dl.ftp_file = S(i).ftp_file;
    dl.dirs     = {S(i).planned.dir};
    g1(c) = g1(c) + 1;
    b1{g1(c),c}.cfg_basicio.run_ops.call_matlab.inputs{1}.evaluated = dl;
    b1{g1(c),c}.cfg_basicio.run_ops.call_matlab.outputs = {};
    b1{g1(c),c}.cfg_basicio.run_ops.call_matlab.fun = 'import_download';
    for k = 1:numel(S(i).planned)
        g1(c) = g1(c) + 1;
        b1{g1(c),c}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = S(i).planned(k).dir;
        b1{g1(c),c}.cfg_basicio.run_ops.call_matlab.outputs = {};
        b1{g1(c),c}.cfg_basicio.run_ops.call_matlab.fun = 'import_convert';
    end
end
t1 = tic;
fprintf('Batch 1: download and DICOM conversion ...\n');
run_batches(b1,n_net,use_par);
timings.batch1 = toc(t1);

%% Serial step: headers, run numbers, names --> batch 2
ctx.path     = path;
ctx.scanner  = scanner;
ctx.acq_type = acq_type;
ctx.dummies  = import.dummies;
ctx.task     = vars.task;
ctx.tr       = vars.sliceTiming.tr;
b2 = cell(1,n_col);
g2 = zeros(1,n_col);
json_done = {};
for i = 1:numel(S)
    try
        fprintf('\nsub #%d session %d (%s)\n',S(i).sub,S(i).ses,S(i).exam_id);
        [S(i).planned,jobs,S(i).out_paths,S(i).out_txt,json_done] = build_session(S(i),ctx,json_done);
    catch err
        S(i).ok  = false;
        S(i).msg = err.message;
        warning('sub #%d session %d FAILED: %s (temp folders are kept)',S(i).sub,S(i).ses,err.message);
        continue
    end
    c = S(i).col;
    for j = 1:numel(jobs)
        g2(c) = g2(c) + 1;
        b2{g2(c),c} = jobs{j};
    end
end

%% Batch 2: 4D files, moves, gzip copies, clean up
t2 = tic;
if any(g2 > 0)
    cols_used = find(g2 > 0);
    fprintf('\nBatch 2: 4D files, moving, gzip copies and clean up ...\n');
    run_batches(b2(:,cols_used),min(numel(cols_used),vars.max_procs),use_par);
end
timings.batch2 = toc(t2);

%% protocol and summary
fprintf('\nImport summary:\n');
for i = 1:numel(S)
    txt = S(i).out_txt;
    if S(i).ok
        for k = 1:numel(S(i).out_paths)
            for j = 1:numel(S(i).out_paths{k})
                if ~isfile(S(i).out_paths{k}{j})
                    txt{k}{j} = [txt{k}{j} '  [MISSING]'];
                    S(i).ok   = false;
                    S(i).msg  = 'some files are missing after the second batch, see the protocol';
                end
            end
        end
    end
    write_protocol(fullfile(path.preprocDir,S(i).sub_str,sprintf('%s_protocol.txt',S(i).exam_id)),S(i).sub_str,S(i).ses_str,S(i).exam_id,S(i).planned,S(i).skipped,txt,timings,S(i).msg);
    if S(i).ok
        fprintf('  %s %s (%s): OK\n',S(i).sub_str,S(i).ses_str,S(i).exam_id);
    else
        fprintf('  %s %s (%s): FAILED - %s\n',S(i).sub_str,S(i).ses_str,S(i).exam_id,S(i).msg);
    end
end
fprintf('Time: download and conversion %.0f s, 4D files and copies %.0f s\n',timings.batch1,timings.batch2);
end

%% ------------------------------------------------------------------------
function run_batches(b,n_procs,use_par)
% matlabbatch{gi,sub} --> one column per subject
% A job that fails makes spm_jobman throw an error: this must not stop the import of
% the other subjects. Missing results are detected afterwards (converted images before
% batch 2, BIDS files after batch 2) and reported per session.
try
    if use_par && n_procs > 1
        run_spm_parallel(b,n_procs);
    else
        run_spm_sequential(b);
    end
catch err
    warning('A job of the batch failed (%s). Checking which sessions are affected ...',err.message);
end
end

function [planned,jobs,out_paths,out_txt,json_done] = build_session(s,ctx,json_done)
% after batch 1: check, read DICOM headers, assign run numbers, build names and the
% jobs of batch 2 (cell row). The json files are written here.
planned  = s.planned;
c_sub    = s.sub;
ses      = s.ses;
out_base = fullfile(ctx.path.preprocDir,s.sub_str,s.ses_str);
zip_base = fullfile(ctx.path.baseDir,s.sub_str,s.ses_str);

for k = 1:numel(planned)
    if isempty(dir(fullfile(planned(k).dir,'*.nii')))
        error('series %d (%s): no converted images, the download or the conversion failed',planned(k).ser,strtrim(planned(k).name));
    end
end

ref_hdr = cell(numel(planned),1); % one DICOM header per label (brain / spinal) for the task json
for k = 1:numel(planned)
    if ~strcmp(planned(k).mod,'func'), continue; end
    all_dicom = spm_select('FPlist',planned(k).dir,'^MR.*');
    [planned(k).subtype,ref_hdr{k}] = read_series_headers(all_dicom,ctx.acq_type);
end
[~,o]   = sort([planned.ser]); % acquisition order (series number); keep ref_hdr in step
planned = planned(o);
ref_hdr = ref_hdr(o);
planned = assign_runs(planned);
print_import_plan(planned,s.skipped,ctx.task,ctx.acq_type);

jobs       = {};
out_paths  = repmat({{}},numel(planned),1); % BIDS files written for each series (for the protocol)
out_txt    = repmat({{}},numel(planned),1);
echo_times = [];
for k = 1:numel(planned)
    p       = planned(k);
    out_dir = fullfile(out_base,p.mod);
    zip_dir = fullfile(zip_base,p.mod);
    if ~isfolder(out_dir), mkdir(out_dir); end
    if ~isfolder(zip_dir), mkdir(zip_dir); end
    switch p.mod
        case 'func'
            task = p.task;
            if isempty(task), task = ctx.task; end
            ff    = dir(fullfile(p.dir,'*.nii'));
            files = fullfile(p.dir,{ff.name})';
            if numel(ctx.acq_type) == 2
                labels = classify_geometry(nifti_signature(files),ctx.acq_type);
            else
                labels = classify_geometry(zeros(numel(files),3),ctx.acq_type);
            end
            for lab = unique(labels,'stable')'
                sel = files(strcmp(labels,lab{1}));
                sel(1:ctx.dummies) = []; % remove dummies
                b   = spm_jsonread(spm_file(sel{end},'ext','json'));
                dat_tim = [datestr(b.acqpar.AcquisitionDate) ' ' char(duration(0,0,b.acqpar.AcquisitionTime))]; % date time str of last EPI
                tag    = acq_tag(lab{1},ctx.acq_type);
                target = fullfile(out_dir,bids_name(c_sub,ses,p.type,'task',task,'run',p.run,'acq',tag,'part',p.part));
                jobs   = [jobs make_4d_jobs(sel,target,dat_tim,ctx.tr,zip_dir)]; %#ok<AGROW>
                out_paths{k}{end+1} = target;
                out_txt{k}{end+1}   = sprintf('%s  (%d vol)',spm_file(target,'filename'),numel(sel));
                json_key = [task '|' tag '|' p.type];
                if ~ismember(json_key,json_done)
                    if isstruct(ref_hdr{k}) && isfield(ref_hdr{k},lab{1})
                        write_task_json(ctx.path.preprocDir,task,tag,p.type,ref_hdr{k}.(lab{1}));
                        json_done{end+1} = json_key; %#ok<AGROW>
                    else
                        warning('No DICOM header of a %s volume among the first ones of series %d, task json not written',lab{1},p.ser);
                    end
                end
            end
        case 'fmap'
            source = spm_select('FPlist',p.dir,sprintf('^s%s.*\\.nii$',ctx.scanner));
            switch p.type
                case {'epi_std','epi_inv'}
                    one_image(source,p);
                    if strcmp(p.type,'epi_std'), dir_lab = 'std'; pe = 'j+'; else, dir_lab = 'inv'; pe = 'j-'; end
                    target = fullfile(out_dir,bids_name(c_sub,ses,'epi','dir',dir_lab,'run',p.run));
                    descr  = struct('PhaseEncodingDirection',pe,'TotalReadoutTime',NaN,'B0FieldIdentifier','pepolar_fmap0');
                    write_json_twice(target,zip_dir,descr);
                    [jobs,out_paths,out_txt] = add_image(jobs,out_paths,out_txt,k,source,target,out_dir,zip_dir);
                case 'magnitude'
                    echo_times = [];
                    for si = 1:size(source,1)
                        target = fullfile(out_dir,bids_name(c_sub,ses,sprintf('%s%d',p.type,si),'run',p.run));
                        b = spm_jsonread(spm_file(source(si,:),'ext','json'));
                        echo_times(si) = b.acqpar.EchoTime./1000; %#ok<AGROW> % echo time in s
                        [jobs,out_paths,out_txt] = add_image(jobs,out_paths,out_txt,k,source(si,:),target,out_dir,zip_dir);
                    end
                case 'phasediff'
                    one_image(source,p);
                    if numel(echo_times) < 2
                        error('no echo times available (series %d): the magnitude series has to be acquired BEFORE the phasediff series',p.ser);
                    end
                    target = fullfile(out_dir,bids_name(c_sub,ses,p.type,'run',p.run));
                    descr  = struct('EchoTime1',echo_times(1),'EchoTime2',echo_times(2),'B0FieldIdentifier','phasediff_fmap0');
                    write_json_twice(target,zip_dir,descr);
                    [jobs,out_paths,out_txt] = add_image(jobs,out_paths,out_txt,k,source,target,out_dir,zip_dir);
                otherwise
                    warning('fmap type %s not implemented, series %d skipped',p.type,p.ser);
            end
        case 'anat'
            source = spm_select('FPlist',p.dir,sprintf('^s%s.*\\.nii$',ctx.scanner));
            one_image(source,p);
            target = fullfile(out_dir,bids_name(c_sub,ses,p.type,'run',p.run));
            [jobs,out_paths,out_txt] = add_image(jobs,out_paths,out_txt,k,source,target,out_dir,zip_dir);
        otherwise
            warning('modality %s not implemented, series %d skipped',p.mod,p.ser);
    end
    % remove the temp dir of this series (after the jobs that use it)
    jobs{end+1}.cfg_basicio.file_dir.dir_ops.dir_move.dir = {p.dir}; %#ok<AGROW>
    jobs{end}.cfg_basicio.file_dir.dir_ops.dir_move.action.delete = true;
end
end

function one_image(source,p)
if size(source,1) ~= 1
    error('series %d (%s): %d converted images found, one expected',p.ser,strtrim(p.name),size(source,1));
end
end

function write_json_twice(target,zip_dir,descr)
% json next to the image in preprocDir and next to the gzip copy in baseDir
js = spm_file(target,'ext','json');
spm_jsonwrite(js,descr,struct('indent','  '));
spm_jsonwrite(fullfile(zip_dir,spm_file(js,'filename')),descr,struct('indent','  '));
end

function [jobs,out_paths,out_txt] = add_image(jobs,out_paths,out_txt,k,source,target,out_dir,zip_dir)
% move a converted image to its BIDS name and gzip a copy into baseDir
jobs{end+1}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(source);
jobs{end}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = {out_dir};
jobs{end}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(1).pattern = '^.*$'; % the whole name ...
jobs{end}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep(1).repl = spm_file(target,'basename'); % ... is replaced
jobs{end}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
jobs{end+1}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.files  = cellstr(target);
jobs{end}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.outdir = {zip_dir};
jobs{end}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.keep   = true;
out_paths{k}{end+1} = target;
out_txt{k}{end+1}   = spm_file(target,'filename');
end

function jobs = make_4d_jobs(files,target,dat_tim,tr,zip_dir)
% three jobs: concatenate to a 4D nifti, add date/time, gzip a copy to zip_dir
jobs = cell(1,3);
jobs{1}.spm.util.cat.vols  = cellstr(files);
jobs{1}.spm.util.cat.name  = target;
jobs{1}.spm.util.cat.dtype = 4;
jobs{1}.spm.util.cat.RT    = tr;
jobs{2}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = target;
jobs{2}.cfg_basicio.run_ops.call_matlab.inputs{2}.string = dat_tim;
jobs{2}.cfg_basicio.run_ops.call_matlab.outputs = {};
jobs{2}.cfg_basicio.run_ops.call_matlab.fun = 'add_time_epi';
jobs{3}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.files  = cellstr(target);
jobs{3}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.outdir = cellstr(zip_dir);
jobs{3}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.keep   = true;
end

function check_task_labels(vars,import,path)
% BIDS labels may only contain letters and digits. vars.task is used as the task
% label in all file names (and as Name in dataset_description.json): with an
% underscore (or hyphen) spm_BIDS can no longer parse the file names and ALL
% calls to spm_BIDS fail. The same holds for the optional import.data(n).task.
labels = {vars.task};
if isfield(import,'data') && isfield(import.data,'task')
    labels = [labels {import.data.task}];
end
dd_file = fullfile(path.preprocDir,'dataset_description.json'); % written by an earlier import
if isfile(dd_file)
    dd = spm_jsonread(dd_file);
    if isfield(dd,'Name'), labels{end+1} = dd.Name; end
end
labels = labels(~cellfun(@isempty,labels));
bad    = labels(~cellfun(@isempty,regexp(labels,'[^a-zA-Z0-9]','once')));
if ~isempty(bad)
    error(['The name "%s" contains characters other than letters and digits (e.g. an underscore).\n' ...
        'vars.task (and import.data(n).task) is used as task label in all BIDS file names and as Name in dataset_description.json,\n' ...
        'and spm_BIDS cannot read such data. Please change it in get_study_specs.m (e.g. to "%s")\n' ...
        'and, if data were already imported, delete or correct the Name in %s.'],...
        bad{1},regexprep(bad{1},'[^a-zA-Z0-9]',''),dd_file);
end
end

function import = read_import_specs(import)
% import.prisma is either given directly or a participants.tsv that
% contains the columns participant_id and scan_id (PRISMA numbers).
% If there are more sessions separate the PRISMA numbers by spaces NOT tabs.
% Series that should NOT be imported are added in () with NO space, e.g.
% participant_id	sex	age	scan_id
% sub-09	F	20	25879 25878(13)
% imports sub-09, with sessions 25879 and 25878. Series 13 in 25878 is omitted
if ischar(import.prisma) && isfile(import.prisma)
    x = spm_load(import.prisma);
    if isnumeric(x.scan_id)
        x.scan_id = cellstr(num2str(x.scan_id));
    end
    import.prisma = [];
    import.prisma_no = nan(1,numel(x.participant_id));
    for ii = 1:numel(x.participant_id)
        import.prisma_no(ii) = sscanf(upper(x.participant_id{ii}),'SUB-%d');
        [import.prisma{ii},import.prisma_skip{ii}] = parse_sessions(x.scan_id{ii});
    end
end
if ~isfield(import,'prisma_skip') % nothing to skip
    import.prisma_skip = cellfun(@(c) repmat({{}},1,numel(c)),import.prisma,'UniformOutput',false);
end
end

function machina = get_machine
if ispc
    machina = 'win_loc';
elseif isunix
    [status, ~] = system('netapp dicq'); % 1 if netapp dicq is found, i.e. we are on a server
    if status == 1
        machina = 'lin_ser';
    else
        machina = 'lin_loc';
    end
end
end

function fetch_first_files(machina,conn,pw,remote_paths,dest_dirs,n)
% download only the first n DICOM files (MR.*) of several series directories,
% with ONE sftp connection for listing and ONE for downloading. The listing of
% each directory is found by the echoed command ("sftp> ls <path>/MR.*") in the
% output; directories where this does not work are listed one by one.
remote_paths = cellstr(remote_paths);
dest_dirs    = cellstr(dest_dirs);
names        = cell(size(remote_paths));

ls_file = [tempname '_ls.txt'];
spm_save(ls_file,strcat('ls',{' '},remote_paths,'/MR.*')); % no options: psftp (PuTTY) does not know ls -1
out = sftp_cmd(machina,conn,pw,ls_file,true);
delete(ls_file);
cur = 0;
for ln_c = strsplit(out,{'\r\n','\n'})
    ln = ln_c{1};
    echo = regexp(ln,'^\s*p?sftp>\s*ls\s+(\S+)','tokens','once');
    if ~isempty(echo)
        cur = find(strcmp(remote_paths,regexprep(echo{1},'/MR\.\*$','')),1);
    elseif cur > 0
        names{cur} = [names{cur} regexp(ln,'MR\.[^\s]+','match')];
    end
end
for i = 1:numel(remote_paths)
    if isempty(names{i}) && numel(remote_paths) > 1 % could not be attributed, try this one alone
        ls_file = [tempname '_ls.txt'];
        spm_save(ls_file,{['ls ' remote_paths{i} '/MR.*']});
        out1 = sftp_cmd(machina,conn,pw,ls_file,true);
        delete(ls_file);
        names{i} = regexp(out1,'MR\.[^\s]+','match');
        if isempty(names{i}), out = out1; end
    end
end

cmds = {};
for i = 1:numel(remote_paths)
    nm = sort(unique(names{i}));
    nm = nm(1:min(n,numel(nm)));
    if isempty(nm)
        warning('No DICOM files (MR.*) found in %s\nsftp said:\n%s',remote_paths{i},out(1:min(end,600)));
        continue
    end
    cmds = [cmds, {['lcd ' dest_dirs{i}],['cd ' remote_paths{i}]}, strcat('get',{' '},nm)]; %#ok<AGROW>
end
if ~isempty(cmds)
    get_file = [tempname '_get.txt'];
    spm_save(get_file,cmds);
    sftp_cmd(machina,conn,pw,get_file,true);
    delete(get_file);
end
end

function [sess,excl] = parse_sessions(str)
% parse multiple session and possible exclusions
% str = '24170(15,25) 24192 24532(3)' becomes
% sess = {[24170]}    {[24192]}    {[24170]}
% excl = {[15 25]}    {0x0 cell}    {[3]}

str = strtrim(str);
err = strfind(str,' (');
if ~isempty(err)
    error('No space between PRISMA no and (exception): %s',str);
end
sessions = strsplit(str);
for i=1:numel(sessions)
    if isempty(strfind(sessions{i},'('))
        sess{i} = str2double(sessions{i});
        excl{i} = {};
    else
        sess{i} = str2double(sessions{i}(1:strfind(sessions{i},'(')-1));
        newStr  = extractBetween(char(sessions(i)),'(',')');
        excl{i} = str2double(strsplit(char(newStr),','));
    end
end
end
