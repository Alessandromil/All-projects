function Coord = calcCoord(Vector)
% This function pass vector from XML form to cell array
% In
%   Vector = char with the form [ ; ; ]
% Out
%   Coord = char in the form cell{3,1} 

PosSeparator = findstr(';',Vector);
Coord{1,1} = Vector(2:(PosSeparator(1)-1));
Coord{2,1} = Vector((PosSeparator(1)+1):(PosSeparator(2)-1));
Coord{3,1} = Vector((PosSeparator(2)+1):(end-1));
end