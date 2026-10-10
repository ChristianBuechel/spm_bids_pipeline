function write_task_json(preproc_dir,task,tag,type,hdr)

% function write_task_json(preproc_dir,task,tag,type,hdr)
% writes task-<task>[_acq-<tag>]_<type>.json (slice timing etc.) from a DICOM header

st_ind = find(strcmp({hdr.CSAImageHeaderInfo.name},'MosaicRefAcqTimes'),1);
if ~isempty(st_ind) && ~isempty(hdr.CSAImageHeaderInfo(st_ind).item) % not for 3D EPI sequences
    slt = str2num(char(hdr.CSAImageHeaderInfo(st_ind).item.val))./1000; %#ok<ST2NM>
    rt  = hdr.RepetitionTime./1000;
else
    slt = 0;
    rt  = 0;
end
trt  = 1./hdr.Private_0019_1028; % 1./BandwidthPerPixelPhaseEncode (Hz)
js   = struct('RepetitionTime',rt,'TaskName',task,'SliceTiming',slt,'TotalReadoutTime',trt);
if isempty(tag), acq = ''; else, acq = ['_acq-' tag]; end
spm_jsonwrite(fullfile(preproc_dir,sprintf('task-%s%s_%s.json',task,acq,type)),js,struct('indent','  '));
end
