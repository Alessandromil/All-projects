function ItemInVector = isInVector(Name,Vector)
% if Item is in vector retun 1 else 0
nItems = size(Vector,1);
ItemInVector = 0;
for i= 1:nItems
    if strcmpi(Name,Vector(i).Name)
        ItemInVector = 1;
        return
    end
end
