function sbp_shift_images(all_sub_ids)

% function shift_images(all_sub_ids)
% allows to move images (bold and T1w) by specifying a 4x4 matrix 
% this can be helpful when for any reason images are too far away from MNI
% space e.g. in combined spinal-brain data

% resolve paths
specs    = get_study_specs;
path     = specs.path;
vars     = specs.vars;
import   = specs.import;

BIDS                  = spm_BIDS(path.preprocDir);
n_subs                = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    
    %do T1w file
    struc_file = get_anat(BIDS,sub_id);
    
    check = nifti(struc_file);
    if sum(sum((check.mat-check.mat0).^2)) > 1e-3 %it seems the combined niftis have a diff
        fprintf('nifti.mat0 and nifti.mat are different, %s seems to have already been moved - SKIP\n',char(struc_file));
    else
        MM = spm_get_space(char(struc_file));
        spm_get_space(char(struc_file),import.shift*MM);
    end
    
    %do bold data
    for ses = 1:vars.nSess
        for run = get_run_numbers(BIDS,sub_id,ses)
            epi                =  get_func_files(BIDS,sub_id,ses,run);
            check =  nifti(epi);
            if isfield(check.extras,'mat')
                fprintf('nifti.extras has .mat entry, %s seems to have already been moved - SKIP\n',char(epi));
            elseif sum(sum((check.mat-check.mat0).^2)) > 1e-3 %it seems the combined niftis have a diff
                fprintf('nifti.mat0 and nifti.mat are different, %s seems to have already been moved - SKIP\n',char(epi));
            else
                 check.mat = import.shift*check.mat; % shift 
                 create(check); % and write ... doing it like this avoids a sidecar *.mat file ...
            end
        end
    end
    
end

end






