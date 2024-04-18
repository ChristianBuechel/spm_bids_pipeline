function import_data_ses
%
% Finds scans using dicq, use SFTP to get them, DICOM convert and make 4D
% and arrange a la BIDS
% who an what to import is specified in get_study_specs
[path,vars,~,import] = get_study_specs;
if isfield(import,'scanner')
    scanner = import.scanner; %if this is TRIO we download TRIO images
else
    scanner = 'PRISMA';
end

do_once = 1;
if ~isunix
    pw = passcode; %assumes that under unix you are already logged on
end
mkdir(path.preprocDir); % create coarse directory structure
spm_save(fullfile(path.preprocDir,'.bidsignore.'),sprintf('**/%s_*_ftp.txt',scanner)); % creates .bidsignore file
descr = struct('Name',vars.task,'BIDSVersion','1.0.2');
spm_jsonwrite(fullfile(path.preprocDir,'dataset_description.json'),descr, struct('indent','  ')); % create description file

for s = 1:numel(import.prisma) % across volunteers
    mkdir(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)));%
    
    for ses = 1:numel(import.prisma{s})
        fi = 1;ftp_cmd = [];
        Volunteer(s).sess(ses).ID = sprintf('%s_%d',scanner,import.prisma{s}{ses});
        ftp_file = fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('%s_ftp.txt',Volunteer(s).sess(ses).ID)); %create SFTP batch file
        fprintf('Doing #%d session %d (%s) \n',import.prisma_no(s),ses,Volunteer(s).sess(ses).ID);
        if isunix
            [status, result] = system(sprintf('ssh %s@%s netapp dicq -f --series --exam=%s',import.user,import.server,Volunteer(s).sess(ses).ID),'-echo'); % query DICOM database
        else
            [status, result] = system(sprintf('plink -ssh %s@%s -pw %s netapp dicq -f --series --exam=%s',import.user,import.server,pw,Volunteer(s).sess(ses).ID),'-echo');
        end
        if status ~= 0
            error('cannot retrieve output from dicq');
        end
        list_of_series = strsplit(result,'\n');
        %% create SFTP batch file
        for tt = 1:numel(import.data)
            run = 1;
            ind = find(contains(list_of_series,import.data(tt).seq)); % find series that match
            for g = ind
                mkdir(fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses)),import.data(tt).dir); % create BIDS dir
                [n,p,ser] = parse_series(list_of_series{g});
                if eval(import.data(tt).cond)
                    mkdir(fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir),sprintf('%s_%d_%d',import.data(tt).type,ses,run)); % create temp dir for dicom files
                    Volunteer(s).sess(ses).data{tt}(run).path   = p;
                    Volunteer(s).sess(ses).data{tt}(run).nscans = n;
                    Volunteer(s).sess(ses).data{tt}(run).series = ser;
                    Volunteer(s).sess(ses).data{tt}(run).mod    = import.data(tt).dir;
                    ftp_cmd{fi} = ['lcd ' fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('%s_%d_%d',import.data(tt).type,ses,run))]; fi = fi + 1;
                    ftp_cmd{fi} = ['cd ' p]; fi = fi + 1;
                    ftp_cmd{fi} = ['mget MR.*']; fi = fi + 1; % create entry into SFTP batch to get all dicom files in that dir
                    run = run + 1;
                end
            end
        end
        spm_save(ftp_file,ftp_cmd);
        %% get the data
        if isunix
            [status, result] = system(sprintf('sftp -b %s %s@%s',ftp_file,import.user,import.server),'-echo'); %simply run the SFTP batch
        else
            [status, result] = system(sprintf('psftp %s@%s -pw %s -b %s',import.user,import.server,pw,ftp_file),'-echo');
        end
        %% now do DICOM conversion within temp dir
        final_del_dir = []; % keep track of temp dirs to later delete them
        matlabbatch  = [];
        g = 1;
        for tt = 1:numel(Volunteer(s).sess(ses).data)
            for run = 1:numel(Volunteer(s).sess(ses).data{tt})
                all_dicom =  spm_select('FPlist', fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('%s_%d_%d',import.data(tt).type,ses,run)),'^MR.*');
                if do_once && strcmp(import.data(tt).dir,'func') % very first epi --> write task-XXX_bold.json file
                    dc   = spm_dicom_headers(all_dicom(1,:));
                    slt  = dc{1}.Private_0019_1029./1000;
                    rt   = dc{1}.RepetitionTime./1000;
                    bold_json = struct('RepetitionTime',rt,'TaskName',vars.task,'SliceTiming',slt);
                    spm_jsonwrite(fullfile(path.preprocDir,sprintf('task-%s_bold.json',vars.task)),bold_json, struct('indent','  '));
                    do_once = 0;
                end
                matlabbatch{g}.spm.util.import.dicom.data = cellstr(all_dicom);
                matlabbatch{g}.spm.util.import.dicom.root = 'flat';
                final_del_dir = strvcat(final_del_dir,fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('%s_%d_%d',import.data(tt).type,ses,run)));
                matlabbatch{g}.spm.util.import.dicom.outdir = {fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('%s_%d_%d',import.data(tt).type,ses,run))};
                matlabbatch{g}.spm.util.import.dicom.protfilter = '.*';
                matlabbatch{g}.spm.util.import.dicom.convopts.format = 'nii';
                matlabbatch{g}.spm.util.import.dicom.convopts.meta = 1;
                matlabbatch{g}.spm.util.import.dicom.convopts.icedims = 0;
                g = g + 1;
            end
        end
        
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);
        
        %% now convert EPIs to 4d nifti (and move them) and simply move other images
        matlabbatch  = [];
        g = 1;
        for tt = 1:numel(Volunteer(s).sess(ses).data)
            for run = 1:numel(Volunteer(s).sess(ses).data{tt})
                if strcmp(Volunteer(s).sess(ses).data{tt}(run).mod,'func')
                    ff = dir(fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('%s_%d_%d',import.data(tt).type,ses,run),'*.nii'));
                    if size(unique(cat(1,ff.bytes)),1) == 2 % indicates we have 2 series eg spinal + brain, where size(spinal)<size(brain)
                        large = cat(1,ff.bytes) > mean(cat(1,ff.bytes));
                        large_f = spm_file(strvcat(ff(large).name),'path',fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('%s_%d_%d',import.data(tt).type,ses,run)));                        
                        large_f(1:import.dummies,:) = []; %remove dummies
                        b = spm_jsonread(spm_file(large_f(end,:),'ext','json'));
                        dat_tim = [datestr(b.acqpar.AcquisitionDate) ' ' char(duration(0,0,b.acqpar.AcquisitionTime))]; % create date time str of last EPI
                        matlabbatch{g}.spm.util.cat.vols = cellstr(large_f);
                        matlabbatch{g}.spm.util.cat.name = fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('sub-%2.2d_ses-%2.2d_task-%s_run-%2.2d_acq-%s_bold.nii',import.prisma_no(s),ses,vars.task,run,'brain'));
                        matlabbatch{g}.spm.util.cat.dtype = 4;
                        matlabbatch{g}.spm.util.cat.RT = vars.sliceTiming.tr;
                        g = g + 1;                        
                        matlabbatch{g}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('sub-%2.2d_ses-%2.2d_task-%s_run-%2.2d_acq-%s_bold.nii',import.prisma_no(s),ses,vars.task,run,'brain'));
                        matlabbatch{g}.cfg_basicio.run_ops.call_matlab.inputs{2}.string = dat_tim;
                        matlabbatch{g}.cfg_basicio.run_ops.call_matlab.outputs = {};
                        matlabbatch{g}.cfg_basicio.run_ops.call_matlab.fun = 'add_time_epi';
                        g = g + 1;                                                
                        matlabbatch{g}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.files = cellstr(fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('sub-%2.2d_ses-%2.2d_task-%s_run-%2.2d_acq-%s_bold.nii',import.prisma_no(s),ses,vars.task,run,'brain')));
                        matlabbatch{g}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.outdir = cellstr(fullfile(path.baseDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir));
                        matlabbatch{g}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.keep = true;
                        g = g + 1;
                        small_f = spm_file(strvcat(ff(~large).name),'path',fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('%s_%d_%d',import.data(tt).type,ses,run)));
                        small_f(1:import.dummies,:) = []; %remove dummies
                        matlabbatch{g}.spm.util.cat.vols = cellstr(small_f);
                        matlabbatch{g}.spm.util.cat.name = fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('sub-%2.2d_ses-%2.2d_task-%s_run-%2.2d_acq-%s_bold.nii',import.prisma_no(s),ses,vars.task,run,'spinal'));
                        matlabbatch{g}.spm.util.cat.dtype = 4;
                        matlabbatch{g}.spm.util.cat.RT = vars.sliceTiming.tr;
                        g = g + 1;                        
                        matlabbatch{g}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('sub-%2.2d_ses-%2.2d_task-%s_run-%2.2d_acq-%s_bold.nii',import.prisma_no(s),ses,vars.task,run,'spinal'));
                        matlabbatch{g}.cfg_basicio.run_ops.call_matlab.inputs{2}.string = dat_tim;
                        matlabbatch{g}.cfg_basicio.run_ops.call_matlab.outputs = {};
                        matlabbatch{g}.cfg_basicio.run_ops.call_matlab.fun = 'add_time_epi';
                        g = g + 1;                                                                       
                        matlabbatch{g}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.files = cellstr(fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('sub-%2.2d_ses-%2.2d_task-%s_run-%2.2d_acq-%s_bold.nii',import.prisma_no(s),ses,vars.task,run,'spinal')));
                        matlabbatch{g}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.outdir = cellstr(fullfile(path.baseDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir));
                        matlabbatch{g}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.keep = true;
                        g = g + 1;
                        
                    else % just one EPI e.g. brain
                        all_f = spm_file(strvcat(ff.name),'path',fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('%s_%d_%d',import.data(tt).type,ses,run)));
                        all_f(1:import.dummies,:) = []; %remove dummies                        
                        b = spm_jsonread(spm_file(all_f(end,:),'ext','json'));
                        dat_tim = [datestr(b.acqpar.AcquisitionDate) ' ' char(duration(0,0,b.acqpar.AcquisitionTime))]; % create date time str of last EPI                        
                        matlabbatch{g}.spm.util.cat.vols = cellstr(all_f);
                        matlabbatch{g}.spm.util.cat.name = fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('sub-%2.2d_ses-%2.2d_task-%s_run-%2.2d_bold.nii',import.prisma_no(s),ses,vars.task,run));
                        matlabbatch{g}.spm.util.cat.dtype = 4;
                        matlabbatch{g}.spm.util.cat.RT = vars.sliceTiming.tr;
                        g = g + 1;                        
                        matlabbatch{g}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('sub-%2.2d_ses-%2.2d_task-%s_run-%2.2d_bold.nii',import.prisma_no(s),ses,vars.task,run));
                        matlabbatch{g}.cfg_basicio.run_ops.call_matlab.inputs{2}.string = dat_tim;
                        matlabbatch{g}.cfg_basicio.run_ops.call_matlab.outputs = {};
                        matlabbatch{g}.cfg_basicio.run_ops.call_matlab.fun = 'add_time_epi';
                        g = g + 1;                        
                        matlabbatch{g}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.files = cellstr(fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('sub-%2.2d_ses-%2.2d_task-%s_run-%2.2d_bold.nii',import.prisma_no(s),ses,vars.task,run)));
                        matlabbatch{g}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.outdir = cellstr(fullfile(path.baseDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir));
                        matlabbatch{g}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.keep = true;
                        g = g + 1;
                    end
                elseif (strcmp(Volunteer(s).sess(ses).data{tt}(run).mod,'anat') || strcmp(Volunteer(s).sess(ses).data{tt}(run).mod,'blip') || strcmp(Volunteer(s).sess(ses).data{tt}(run).mod,'fmap'))
                    source = spm_select('FPlist',fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('%s_%d_%d',import.data(tt).type,ses,run)),sprintf('^s%s.*\\.nii$',scanner));
                    if size(source,1) > 1 % e.g. 2 magnitude images for fieldmap
                        for si=1:size(source,1)
                            target = fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('sub-%2.2d_ses-%2.2d_%s%d.nii',import.prisma_no(s),ses,import.data(tt).type,si));
                            movefile(source(si,:),target);
                            mkdir(fullfile(path.baseDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir));
                            spm_copy(target,fullfile(path.baseDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir),'gzip',true);
                        end
                    else
                        target = fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir,sprintf('sub-%2.2d_ses-%2.2d_%s.nii',import.prisma_no(s),ses,import.data(tt).type));
                        movefile(source,target);
                        mkdir(fullfile(path.baseDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir));
                        spm_copy(target,fullfile(path.baseDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('ses-%2.2d',ses),import.data(tt).dir),'gzip',true);
                    end
                end
            end
        end
        %% remove temp dirs and files inside them
        for dd = 1: size(final_del_dir,1)
            matlabbatch{g}.cfg_basicio.file_dir.dir_ops.dir_move.dir = cellstr(final_del_dir(dd,:));
            matlabbatch{g}.cfg_basicio.file_dir.dir_ops.dir_move.action.delete = true;
            g = g + 1;
        end
        spm_jobman('run',matlabbatch);
    end
end
end

function [n, path, ser] = parse_series(in)  %quite crude but works as dicq is very well formatted ;-)
eos  = findstr(in,[']']);
n    = str2num(in(eos-6:eos-1));
sp   = findstr(in,['/']);
path = in(sp:end);
se   = findstr(in,['{']);
ser  = str2num(in(se-5:se-1));
end
