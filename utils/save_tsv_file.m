function save_tsv_file(filename,crop,samp_int)
%% Create json and tsv logs by runs
% write json header
bold_json = struct('SamplingFrequency',1/samp_int,'StartTime',0,'Columns',["cardiac","respiratory","trigger"]);
spm_jsonwrite(spm_file(filename,'suffix','_physio','ext','.json'),bold_json, struct('indent','  '));

matrix = [crop.cardiac,crop.resp,crop.trigger];
writematrix(matrix,spm_file(filename,'suffix','_physio','ext','.tsv'),'FileType','text','Delimiter','\t');
gzip(spm_file(filename,'suffix','_physio','ext','.tsv'));
%delete(spm_file(filename,'suffix','_physio','ext','.tsv'));

end