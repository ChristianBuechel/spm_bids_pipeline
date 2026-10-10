function epi = get_ref_epi(BIDS,sub,varargin)

% function epi = get_ref_epi(BIDS,sub,varargin)
% returns the reference functional image (session 1, run 1) of a subject as 1x1 cell
% (brain only, see get_func_files for defaults and options)

epi = get_func_files(BIDS,sub,1,1,varargin{:});
end
