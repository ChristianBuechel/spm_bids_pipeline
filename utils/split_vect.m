function chuckCell = split_vect(v,n)
% splits vector into number of n chunks of equal size
% based on lukas function, based on
% http://code.activestate.com/recipes/425044

chuckCell  = {};
%vectLength = numel(v);
vectLength = size(v,2);
splitsize  = 1/n*vectLength;

for i = 1:n
    idxs = [floor(round((i-1)*splitsize)):floor(round((i) * splitsize))-1]+1;
    chuckCell{end + 1} = v(:,idxs);
end

end