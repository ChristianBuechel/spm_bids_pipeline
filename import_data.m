function import_data
%
% Finds scans using dicq, uses SFTP to get them DICOM convert make 4D nifti
% and arrange a la BIDS
[path,vars,analysis,import] = get_study_specs;
min_epi_scans     = import.func.min_vol;
n_mprage_slices   = import.anat.slices; %
if ~isunix
    pw = passcode;
end
% create coarse directory structure
mkdir(path.preprocDir); %
for s = 1:numel(import.prisma)
    fi = 1;ftp_cmd = [];
    mkdir(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s))); %
    Volunteer(s).ID = sprintf('PRISMA_%d',import.prisma(s));
    ftp_file = fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),sprintf('%s_ftp.txt',Volunteer(s).ID));
    fprintf('Doing #%d (%s) \n',import.prisma_no(s),Volunteer(s).ID);
    if isunix
        [status, result] = system(sprintf('netapp dicq -lf --series --exam=%s',Volunteer(s).ID),'-echo');
    else
        [status, result] = system(sprintf('plink -ssh %s@%s -pw %s netapp dicq -lf --series --exam=%s',import.user,import.server,pw,Volunteer(s).ID),'-echo');
    end
    if status == 0
        list_of_series = strsplit(result,'\n');
        
        EPI_ind    =  find(contains(list_of_series,import.func.seq));
        T1_ind     =  find(contains(list_of_series,import.anat.seq));
        
        % now get func
        % ------------
        epi_s = 1;
        for g = EPI_ind
            mkdir(fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s))),'func'); %
            [n,p,ser] = parse_series(list_of_series{g});
            if n > min_epi_scans %session with less were probably aborted
                mkdir(fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),'func'),sprintf('ses%d',epi_s)); % temp_dir for epis
                Volunteer(s).epi(epi_s).path   = p;
                Volunteer(s).epi(epi_s).nscans = n;
                Volunteer(s).epi(epi_s).series = ser;
                ftp_cmd{fi} = ['lcd ' fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),'func',sprintf('ses%d',epi_s))]; fi = fi + 1;
                ftp_cmd{fi} = ['cd ' p]; fi = fi + 1;
                ftp_cmd{fi} = ['mget MR.*']; fi = fi + 1;
                epi_s = epi_s+1;
            end
        end
        
        % now get T1w
        % -----------
        T1_s = 1;
        for g = T1_ind
            mkdir(fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s))),'anat'); %
            [n,p,ser] = parse_series(list_of_series{g});
            if n == n_mprage_slices
                Volunteer(s).T1(T1_s).path   = p;
                Volunteer(s).T1(T1_s).nscans = n;
                Volunteer(s).T1(T1_s).series = ser;
                ftp_cmd{fi} = ['lcd ' fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),'anat')]; fi = fi + 1;
                ftp_cmd{fi} = ['cd ' p]; fi = fi + 1;
                ftp_cmd{fi} = ['mget MR.*']; fi = fi + 1;
                T1_s = T1_s+1;
            end
        end
        spm_save(ftp_file,ftp_cmd);
        if isunix
            [status, result] = system(sprintf('sftp -b %s %s',ftp_file,import.server),'-echo');
        else
            [status, result] = system(sprintf('psftp %s@%s -pw %s -b %s',import.user,import.server,pw,ftp_file),'-echo');
        end
        % now do DICOM conversion
        final_delete = [];
        matlabbatch  = [];
        for g=1:numel(Volunteer(s).epi)
            all_dicom =  spm_select('FPlist', fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),'func',sprintf('ses%d',g)),'^MR.*');
            final_delete = strvcat(final_delete,all_dicom);
            matlabbatch{g}.spm.util.import.dicom.data = cellstr(all_dicom);
            matlabbatch{g}.spm.util.import.dicom.root = 'flat';
            matlabbatch{g}.spm.util.import.dicom.outdir = {fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),'func',sprintf('ses%d',g))};
            matlabbatch{g}.spm.util.import.dicom.protfilter = '.*';
            matlabbatch{g}.spm.util.import.dicom.convopts.format = 'nii';
            matlabbatch{g}.spm.util.import.dicom.convopts.meta = 0;
            matlabbatch{g}.spm.util.import.dicom.convopts.icedims = 0;
        end
        g = g + 1;
        all_dicom =  spm_select('FPlist', fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),'anat'),'^MR.*');
        final_delete = strvcat(final_delete,all_dicom); %gets renamed
        matlabbatch{g}.spm.util.import.dicom.data = cellstr(all_dicom);
        matlabbatch{g}.spm.util.import.dicom.root = 'flat';
        matlabbatch{g}.spm.util.import.dicom.outdir = {fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),'anat')};
        matlabbatch{g}.spm.util.import.dicom.protfilter = '.*';
        matlabbatch{g}.spm.util.import.dicom.convopts.format = 'nii';
        matlabbatch{g}.spm.util.import.dicom.convopts.meta = 0;
        matlabbatch{g}.spm.util.import.dicom.convopts.icedims = 0;
        spm_jobman('initcfg');
        save('test','matlabbatch');
        spm_jobman('run',matlabbatch);
        % now make EPIs 4d nifti and rename T1
        matlabbatch = [];
        for g=1:numel(Volunteer(s).epi)
            all_nifti = spm_select('FPlist', fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),'func',sprintf('ses%d',g)),'^fPRISMA.*\.nii$');
            final_delete = strvcat(final_delete,all_nifti);
            all_nifti(1:import.dummies,:) = []; %remove dummies
            matlabbatch{g}.spm.util.cat.vols = cellstr(all_nifti);
            matlabbatch{g}.spm.util.cat.name = fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),'func',sprintf('sub-%2.2d_ses-01_task-%s_run-%2.2d_bold.nii',import.prisma_no(s),vars.task,g));
            matlabbatch{g}.spm.util.cat.dtype = 4;
            matlabbatch{g}.spm.util.cat.RT = NaN;
        end
        g = g + 1;
        % move and rename T1
        movefile(spm_select('FPlist', fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),'anat'),'^sPRISMA.*\.nii$'),fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),'anat',sprintf('sub-%2.2d_T1w.nii',import.prisma_no(s))));
        % remove files        
        matlabbatch{g}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(final_delete);
        matlabbatch{g}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false; %  
        g = g + 1;
        % remove dirs        
        for dd = 1: numel(Volunteer(s).epi)
            matlabbatch{g}.cfg_basicio.file_dir.dir_ops.dir_move.dir = cellstr(fullfile(path.preprocDir,sprintf('sub-%2.2d',import.prisma_no(s)),'func',sprintf('ses%d',dd)));
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
