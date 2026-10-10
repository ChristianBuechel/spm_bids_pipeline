function copy_sidecar(src_images,dst_images)

% function copy_sidecar(src_images,dst_images)
% copies the realignment sidecar (<image>.mat, the header of every volume of a 4D image
% after realignment/coregistration) of the images in src_images to the images in
% dst_images (cellstr, same length). Used for phase images, which are acquired at the same
% time as their magnitude images and cannot be realigned themselves (phase wraps).

src_images = cellstr(src_images);
dst_images = cellstr(dst_images);
if numel(src_images) ~= numel(dst_images)
    error('copy_sidecar: %d source and %d destination images',numel(src_images),numel(dst_images));
end
for i = 1:numel(src_images)
    src = spm_file(src_images{i},'ext','.mat');
    dst = spm_file(dst_images{i},'ext','.mat');
    if ~isfile(src)
        error('copy_sidecar: %s not found (was the magnitude image realigned?)',src);
    end
    copyfile(src,dst,'f');
end
end
