function check_realign(all_sub_ids)

% function check_realign(all_sub_ids);
% displays 1st images from all sessions in coreg to see whether everything is OK


% resolve paths
[path,vars]  = get_study_specs;
BIDS         = spm_BIDS(path.preprocDir);
n_subs       = length(all_sub_ids);

for sub = 1:n_subs
    sub_id      = all_sub_ids(sub);
        
    alla_niftis = [];
    cnt = 1;
    for ses = 1:vars.nSess
        for run = get_run_numbers(BIDS,sub_id,ses)
            epi                =  get_func_files(BIDS,sub_id,ses,run);
            run_niftis         =  spm_select('ExtFPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^ra'),1);
            runa_niftis        =  spm_select('ExtFPlist', spm_file(epi,'path'), spm_file(spm_file(epi,'filename'),'prefix','^a'),Inf);
            all_niftis{cnt,1}  =  run_niftis;
            alla_niftis        =  strvcat(alla_niftis, char(runa_niftis));
            cnt = cnt + 1;
        end
    end

    matlabbatch{1}.spm.util.checkreg.data = all_niftis;

    plot_realign(alla_niftis);
    spm_jobman('run',matlabbatch);
    spm_orthviews('Reposition',[0 0 0]);
    spm_orthviews('zoom',Inf);
    spm_orthviews('redraw');

    fprintf('Checking Subject %d \nPress enter in command window to continue\n',sub_id);
    input('');
end

end






