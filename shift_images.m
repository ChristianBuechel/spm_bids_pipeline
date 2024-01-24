function shift_images(all_sub_ids)

% function shift_images(all_sub_ids)
% allows to move images (bold and T1w) by specifying a 4x4 matrix 
% this can be helpful when for any reason images are too far away from MNI
% space e.g. in combined spinal-brain data

% resolve paths
[path,vars,~,import]  = get_study_specs;
BIDS                  = spm_BIDS(path.preprocDir);
n_subs                = length(all_sub_ids);

for sub = 1:n_subs
    sub_id     = all_sub_ids(sub);
    
    %do T1w file
    struc_file = spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'type','T1w');
    struc_file = struc_file(1); % if there are more than 1
    
    check = nifti(struc_file);
    if sum(sum((check.mat-check.mat0).^2)) > 1e-3 %it seems the combined niftis have a diff
        fprintf('nifti.mat0 and nifti.mat are different, %s seems to have already been moved - SKIP\n',char(struc_file));
    else
        MM = spm_get_space(struc_file);
        spm_get_space(struc_file,import.shift*MM);
    end
    
    %do bold data
    for ses = 1:vars.nSess
        for run = 1:vars.nRuns
            epi                =  spm_BIDS(BIDS,'data','sub',sprintf('%02d',sub_id),'ses',sprintf('%02d',ses),'run',sprintf('%02d',run),'task',vars.task,'type','bold');
            epi                =  epi(1); % if there are brain and spinal, just take brain
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






