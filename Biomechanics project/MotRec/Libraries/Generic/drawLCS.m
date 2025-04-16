function drawLCS(BodyData, varargin)

% Plots local coordinates of points (markers & landmarks)
%
%    drawLCS(BodyData)
%    drawLCS(BodyData, CreateNewFigure)
%    drawLCS(BodyData, CreateNewFigure, FigureTitle)
%
%  + BodyData: struct with five fields:
%      BodyData.Name:         Body name
%      BodyData.MarkerNames:  cell with Marker names
%      BodyData.MarkerCoords: array (3 x nMarkers) with Markers local coordinates
%  + CreateNewFigure: 1 creates new figure, 0 no new figure
%  + FigureTitle: string with figure title
%
%  Example of body definition and usage of drawLCS:
%      AntebrazoData.Name         = 'Antebrazo';
%      AntebrazoData.MarkerNames  = {'M007', 'M008', 'M009', 'M010'};
%      AntebrazoData.MarkerCoords = [Antebrazo_Pos_M007, Antebrazo_Pos_M008, Antebrazo_Pos_M009, Antebrazo_Pos_M010];
%      drawLCS(AntebrazoData, 1);


% Body segment properties
BodyName     = BodyData.Name; 
MarkerNames  = BodyData.MarkerNames; 
MarkerCoords = BodyData.MarkerCoords;

if nargin == 1
    CreateNewFigure = 0;
    FigureTitle     = BodyName;    
elseif nargin == 2 
    CreateNewFigure = varargin{1};
    FigureTitle     = BodyName;
elseif nargin == 3
    CreateNewFigure = varargin{1};
    FigureTitle     = varargin{2};

end


FontSize = 12;
if CreateNewFigure
    figure('Name',FigureTitle,'NumberTitle','off','WindowStyle','docked'), axes('FontSize',FontSize)    
end
hold on, axis equal, grid on

% find scale
MaxCoord  = max(max(abs(MarkerCoords)));
LCSlength = 1.2 * MaxCoord;

% plot axes vectors
quiver3(0,0,0,LCSlength,0,0,0,'r','LineWidth',2);
quiver3(0,0,0,0,LCSlength,0,0,'g','LineWidth',2);
quiver3(0,0,0,0,0,LCSlength,0,'b','LineWidth',2);

% plot axes names
text(LCSlength,0,0,'X_L_C_S   ','FontSize',12,'FontWeight','Bold','Color','r','HorizontalAlignment','right'); % X-axis
text(0,LCSlength,0,'Y_L_C_S   ','FontSize',12,'FontWeight','Bold','Color','g','HorizontalAlignment','right'); % Y-axis
text(0,0,LCSlength,'Z_L_C_S   ','FontSize',12,'FontWeight','Bold','Color','b','HorizontalAlignment','right'); % Z-axis

% add title to figure
title(['Segment ',BodyName]);

% Plot calculated Markers referred to LCS
nMarkers = size(MarkerCoords,2);
for i = 1 : nMarkers
    plot3(MarkerCoords(1,i), MarkerCoords(2,i), MarkerCoords(3,i), 'bo');
    text(MarkerCoords(1,i), MarkerCoords(2,i), MarkerCoords(3,i), ['  ',MarkerNames{i}], 'FontSize', 10, 'Color', 'b');
    line([0, MarkerCoords(1,i)], [0,MarkerCoords(2,i)], [0,MarkerCoords(3,i)], 'Color', 'b');
end

view([155 10]);
