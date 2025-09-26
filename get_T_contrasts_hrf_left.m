function [con_vec, con_label] = get_T_contrasts_hrf_left(conditions)

% define vectors for experimental factors
% remove the last 5 conditions as they are "nuisance" variables
n_nuis = 5;
test   = conditions(1:end-n_nuis);

% parse the conditions for exp factors
stim = contains(test,'STM').*2-1;
move = (contains(test, '-Move') | contains(test, '_Move')).*2-1;
pain = contains(test, '_Pain').*2-1;

% main effects (NO right factor anymore)
all_con = [pain; move; stim];
labels  = {{'Pain','Nopain'}, {'Move','Nomove'}, {'Stim','Cue'}}; % same order

% masks (NO right factor)
all_mask = [(stim<0); (stim>0); (pain>0); (pain<0); (move>0); (move<0)];
m_labels = {'CUE_', 'STM_', 'PAIN_', 'NOPAIN_', 'MOV_', 'NOMOV_'};

%% main effects selection
% (adapted indices: only 3 factors left)
me      =  [ 1     1     2     2     2     2      1      1];
me_mask = {[1],  [2],  [1],  [2], [2 3], [2 4], [1 6], [2 6]}; 

% interactions
inter      = {[1 2], [1 2]};
inter_mask = {[1],   [2]};

%% build contrasts
con_vec = [];
con_label = [];

for i_me = 1:numel(me)
    % positive
    con_vec   = [con_vec; prod([(all_con(me(i_me),:)>0); prod(all_mask(me_mask{i_me},:),1)])];
    con_label = strvcat(con_label,[strcat(m_labels{me_mask{i_me}}) labels{me(i_me)}{1}]);
    % negative
    con_vec   = [con_vec; prod([(all_con(me(i_me),:)<0); prod(all_mask(me_mask{i_me},:),1)])];
    con_label = strvcat(con_label,[strcat(m_labels{me_mask{i_me}}) labels{me(i_me)}{2}]);
    % diff
    con_vec   = [con_vec; prod([all_con(me(i_me),:); prod(all_mask(me_mask{i_me},:),1)])];
    con_label = strvcat(con_label,[strcat(m_labels{me_mask{i_me}}) labels{me(i_me)}{1} '>' labels{me(i_me)}{2}]);
    con_vec   = [con_vec; prod([all_con(me(i_me),:).*-1; prod(all_mask(me_mask{i_me},:),1)])];
    con_label = strvcat(con_label,[strcat(m_labels{me_mask{i_me}}) labels{me(i_me)}{1} '<' labels{me(i_me)}{2}]);
end

% interactions
for i_int = 1:numel(inter)
    con_vec   = [con_vec; prod([prod(all_con(inter{i_int},:)); prod(all_mask(inter_mask{i_int},:),1)])];
    con_vec   = [con_vec; prod([-prod(all_con(inter{i_int},:)); prod(all_mask(inter_mask{i_int},:),1)])];
    
    t_cname = strcat(m_labels{inter_mask{i_int}});
    for i_comp = 1:numel(inter{i_int})
        t_cname = [t_cname labels{inter{i_int}(i_comp)}{1} '_' labels{inter{i_int}(i_comp)}{2} '_x_'];
    end
    con_label = strvcat(con_label,t_cname(1:end-3));
    con_label = strvcat(con_label,[t_cname(1:end-3) '_neg']);
end

% add zeros for nuisance
con_vec   = [con_vec zeros(size(con_vec,1),n_nuis)];
% make cellstr
con_label = cellstr(con_label);

end
