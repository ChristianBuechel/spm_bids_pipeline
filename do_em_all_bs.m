function do_em_all_bs
% test the new batch system 

% define subjects that should be included
all_sub_ids  = [5:44]; 
%all_sub_ids     = [1:4 6:9]; %maybe exclude 8


% create a brainstem mask and then do a brainstem specific realignment
% takes a brainstem mask in template space and warps it back to epi space
create_bs_mask(all_sub_ids);

% does realignment with a weight of 1 over the braistem and zeros elsewhere 
realign_bs(all_sub_ids);

end


