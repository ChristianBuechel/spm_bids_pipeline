function do_em_all_bs_physio
% test the new batch system 

% define subjects that should be included
all_sub_ids  = [5:44]; 
all_sub_ids  = [5]; 

% create a brainstem mask and then do a brainstem speciic realignemtn
% takes a brainstemmask in template space and warps it back to epi space
% create_bs_mask(all_sub_ids);


% create noise regressors from CSF around brainstem
% create_noise_bs(all_sub_ids);

% remove noise surces from EPIs before registering (note that the original
% files will be resliced!)
% CSF regressors from brainstem will be removed, if present physio as well

% clean_physio(all_sub_ids);

% does realignment with a weight of 1 over the braistem and zeros elsewhere 
 realign_bs_physio(all_sub_ids);

end


