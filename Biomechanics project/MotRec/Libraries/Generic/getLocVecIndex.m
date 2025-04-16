function Index = getLocVecIndex(Name,Vector)
% GETLOCVECINDEX search the position in the local vector and returns the position,
% if it is not in vector returns zero 
%
%   Index = getLocVecIndex(Name,Vector)
%   Inputs:
%   Name = Is the name of Item.
%   Vector = Is the container of items.
%   Outputs:
%   Index = Pos the Item in container. Equal to 0 if it's not in container.

nItems = size(Vector,1);
for i= 1:nItems
    if strcmpi(Name,Vector(i).Point.Name)
        Index = i;
        return
    end
end
Index = 0;