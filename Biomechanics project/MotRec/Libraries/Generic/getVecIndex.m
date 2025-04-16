function Index = getVecIndex(Name,Vector)
% GETVECINDEX find the position in the vector and returns the index,
% if it is not in vector returns zero 
%
%   Index = getVecIndex(Name,Vector)
%   Inputs:
%   Name = Is the name of Item.
%   Vector = Is the container of items.
%   Outputs:
%   Index = Index of the Item in container. Equal to 0 if it's not in container.

nItems = size(Vector,1);
for i= 1:nItems
    if iscell(Vector)
        if strcmpi(Name,Vector{i}.Name)
            Index = i;
            return
        end
    else
        if strcmpi(Name,Vector(i).Name)
            Index = i;
            return
        end
    end
end
Index = 0;