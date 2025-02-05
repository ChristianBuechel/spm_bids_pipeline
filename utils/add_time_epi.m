function add_time_epi(fname,dat_tim)
% function add_time_epi
% puts the dat_tim string into the description field of the 4D EPI
% usually we use the acq time of the last EPI volume to allow a comparison
% with the psychtoolbox output file
a         = nifti(fname);
a.descrip = dat_tim;
create(a);
