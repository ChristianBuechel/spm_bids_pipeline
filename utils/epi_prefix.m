function p = epi_prefix(vars,stage)

% function p = epi_prefix(vars,stage)
% prefix of the EPI images that a preprocessing STAGE reads, derived from the
% steps that were run before it (vars.preproc.steps, see epi_steps).
%
% Every step has a FIXED output prefix: nordic 'n', slicetiming 'a', realign (reslice) 'r',
% smooth 's'. SPM prepends the prefix, so the name of an image lists the steps in REVERSE
% order: NORDIC followed by slice timing gives 'a' + 'n' = 'an' (an<epi>.nii).
%
%   vars.preproc.steps = {'nordic','slicetiming','realign','smooth'};   % what was run on the EPIs,
%                                                    % in this order. 'realign' can be omitted (it is then
%                                                    % put after the other steps, before 'smooth')
% stage:
%   'nordic'       input prefix of sbp_nordic         (steps before nordic)
%   'slicetiming'  input prefix of sbp_slice_timing   (steps before slice timing)
%   'realign'      input prefix of realignment AND of everything that is built from the realignment.
%                  All these names are built from it (mean<p>, cmean<p>, vdm5_<p>, rp_<p>, ...,
%                  r<p> for the resliced images, see epi_file)
%   'smooth'       input prefix of sbp_smooth (r<p> if the images are smoothed after the realignment)
%   'analysis'     prefix of the images used for the first level analysis: after ALL steps
%                  (e.g. 'rn' or 'srn' if they are smoothed, see epi_file kind 'analysis')
%
% A single stage can be overridden explicitly: vars.input_prefix.<stage> = '...'
% Examples (steps -> realign input prefix / analysis prefix):
%   {'nordic','slicetiming'} -> 'an' / 'ran',  {'nordic','realign','smooth'} -> 'n' / 'srn',
%   {'slicetiming','realign'} -> 'a' / 'ra',  {} -> '' / 'r' (raw images)

out = struct('nordic','n','slicetiming','a','realign','r','smooth','s');

if isfield(vars,'input_prefix') && isstruct(vars.input_prefix) && isfield(vars.input_prefix,stage)
    p = vars.input_prefix.(stage);
    return
end

steps = epi_steps(vars);                         % checks the list

switch stage
    case {'nordic','slicetiming','realign','smooth'}
        k = find(strcmp(steps,stage));
        if isempty(k)
            error('The step "%s" is not in vars.preproc.steps, so it has no input prefix',stage);
        end
        p = prefix_after(steps(1:k-1),out);
    case 'analysis'
        p = prefix_after(steps,out);
    otherwise
        error('Unknown stage "%s" (nordic, slicetiming, realign, smooth, analysis)',stage);
end
end

function p = prefix_after(steps,out)
p = '';
for i = 1:numel(steps)
    p = [out.(steps{i}) p]; % SPM prepends
end
end
