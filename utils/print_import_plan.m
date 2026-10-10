function print_import_plan(planned,skipped,default_task,acq_type)

% function print_import_plan(planned,skipped,default_task,acq_type)
% prints the table of what sbp_import_data imports (and what is not, and why)
%   planned, skipped  from plan_series (optionally with run, part, subtype from assign_runs)
%   default_task      vars.task
%   acq_type          vars.acq_type (or [])

fprintf('  ser     n  modality/type      task        run  part   img    -> file name\n');
for k = 1:numel(planned)
    p = planned(k);
    task = p.task;
    if isempty(task) && strcmp(p.mod,'func'), task = default_task; end
    fname = '';
    if isfield(p,'run') && ~isempty(p.run) && p.run > 0
        if strcmp(p.mod,'func') && ~isempty(acq_type) && numel(acq_type) == 1
            tag = acq_tag(acq_type{1},acq_type);
        elseif strcmp(p.mod,'func') && ~isempty(acq_type)
            tag = '<acq>';
        else
            tag = '';
        end
        fname = regexprep(bids_name(0,[],p.type,'task',task,'run',p.run,'acq',tag,'part',p.part),'^sub-00_',''); % sub/ses are the same for all
    end
    fprintf('  %3d  %4d  %-5s %-12s %-10s %3s  %-5s  %-5s  %s\n',p.ser,p.n,p.mod,p.type,task,num2str(field_or(p,'run','-')),field_or(p,'part',''),field_or(p,'subtype',''),fname);
end
for k = 1:numel(skipped)
    fprintf('  SKIPPED series %d (n = %d): %s\n',skipped(k).ser,skipped(k).n,skipped(k).reason);
end
end

function v = field_or(s,f,d)
if isfield(s,f) && ~isempty(s.(f)) && ~(isnumeric(s.(f)) && s.(f)==0), v = s.(f); else, v = d; end
end
