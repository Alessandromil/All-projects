function NumCoord = calcNumCoord(CharCoord)
% This function pass vector from XML form to cell array
% In
%   CharCoord = char in the form cell{3,1}
% Out
%   NumCoord = double in form [3x1]

NumCoord(1) = str2double(CharCoord{1,1});
NumCoord(2) = str2double(CharCoord{2,1});
NumCoord(3) = str2double(CharCoord{3,1});
NumCoord = NumCoord';
