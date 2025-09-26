function [con_vec, con_label] = get_F_contrasts_hrf_left(conditions)

% define vectors for experimental factors
% remove the last 5 conditions as they are "nuisance" variables
n_nuis = 5;
test   = conditions(1:end-n_nuis);

% parse the conditions for exp factors (no RIGHT anymore)
stim = contains(test,'STM').*2-1;
move = (contains(test, '-Move') | contains(test, '_Move')).*2-1;
pain = contains(test, '_Pain').*2-1;

% main effects (only 3 factors now)
all_con = [pain; move; stim];
labels  = {{'Pain','Nopain'}, {'Move','Nomove'}, {'Stim','Cue'}}; % same order

% masks (drop RIGHT-based masks)
all_mask  = [(stim<0); (stim>0); (pain>0); (pain<0); (move>0); (move<0)];
m_labels  = {'CUE_', 'STM_', 'PAIN_', 'NOPAIN_', 'MOV_', 'NOMOV_'};

% what we want (example: stim-based and move-based main effects)
% NOTE: indices changed because no RIGHT
me      =  [3  3];          
me_mask = {[], []};       

% Example: 
% me      = [1 1 2 2 3 3]; 
% me_mask = {[1], [2], [1], [2], [1], [2]}; 
% would give Pain, Move, Stim contrasts across cue/stim masks

con_vec = {};
con_label = {};
i_fc = 1;

for i_me = 1:numel(me)
    % positive group
    s_f             = diag(prod([(all_con(me(i_me),:)>0); prod(all_mask(me_mask{i_me},:),1)]));
    con_vec{i_fc}   = [s_f(any(s_f),:) zeros(sum(any(s_f)),n_nuis)];
    %con_label{i_fc} = [strcat(m_labels{me_mask{i_me}}) labels{me(i_me)}{1}];
    i_fc = i_fc + 1;

    % negative group
    s_f             = diag(prod([(all_con(me(i_me),:)<0); prod(all_mask(me_mask{i_me},:),1)]));
    con_vec{i_fc}   = [s_f(any(s_f),:) zeros(sum(any(s_f)),n_nuis)];
    %con_label{i_fc} = [strcat(m_labels{me_mask{i_me}}) labels{me(i_me)}{2}];
    i_fc = i_fc + 1;

    % pooled F-test (all levels != 0)
    s_f             = diag(prod([(all_con(me(i_me),:)~=0); prod(all_mask(me_mask{i_me},:),1)]));
    con_vec{i_fc}   = [s_f(any(s_f),:) zeros(sum(any(s_f)),n_nuis)];
    %con_label{i_fc} = [strcat(m_labels{me_mask{i_me}}) labels{me(i_me)}{1} labels{me(i_me)}{2}];
    i_fc = i_fc + 1;
end

% make cellstr
%con_label = cellstr(con_label);

% quick hardcoded fix
con_vec = con_vec(1:2);
con_label = {'STM', 'CUE'};

end
