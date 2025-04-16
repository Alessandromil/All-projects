function BodyCS2_Pos_Points = changeCoordSys(BodyCS1_Pos_Points, d, R)

% CHANGECOORDSYS transforms the position vectors from 
% Coordinate System 1 (CS1) to Coordinate System 2 (CS2)
%
%   BodyCS2_Pos_Points = changeCoordSys(BodyCS1_Pos_Points, d, R)
%
%   Inputs:
%     + BodyCS1_Pos_Points are the coordinates of points referred to 
%       Coordinate system 1 (CS1) - double array(3 x nPoints)
%     + d is the position vector of the CS2 origin 
%       referred to the CS1 - double array(3x1)
%     + R is is the rotation matrix between CS1 and CS2 
%       - double array(3x3)
%   Outputs:
%     + BodyCS2_Pos_Points are the coordinates of points referred to 
%       Coordinate system 2 (CS2) - double array(3 x nPoints)

nPoints = size(BodyCS1_Pos_Points, 2);
BodyCS2_Pos_Points = R * (BodyCS1_Pos_Points - d * ones(1, nPoints));
