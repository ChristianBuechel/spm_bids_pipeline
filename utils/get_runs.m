function [ind,sess] = get_runs(scan)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ind = [0; find(diff(scan)>2*mean(diff(scan))); size(scan,1)];
sess = diff(ind);
end