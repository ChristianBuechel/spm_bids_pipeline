function copy_project(source_dir,dest_dir,doit)
% copies BIDS relevant files from source to dest
% using spm_BIDS
% note that the current version of spm_BIDS misses the magnitude images for
% fieldmaps
if nargin<3
    doit = 0;
end

if source_dir(end) ~= filesep
    source_dir(end+1) = filesep;
end

if dest_dir(end) ~= filesep
    dest_dir(end+1) = filesep;
end

BIDS         = spm_BIDS(source_dir);
source_files = spm_BIDS(BIDS,'data');
dest_files   = strrep(source_files,source_dir,dest_dir);

if doit
    for j=1:numel(source_files)
        [p,~,~] = fileparts(dest_files{j});
        if ~exist(p,'dir')
            mkdir(p);
        end
        copyfile(source_files{j},dest_files{j})
    end
else
    for j=1:numel(source_files)
        fprintf('%s --> %s\n',source_files{j},dest_files{j});
    end
end


end

