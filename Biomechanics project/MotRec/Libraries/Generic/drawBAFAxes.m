function drawBAFAxes(varargin)

% DRAWBAFAXES draws in a figure the skin-markers, landmarks, COM and joint 
% axes (if exist) of a body.
%
%    drawBAFAxes(BodyName, BodyBAF_Pos_Smark, BodyBAF_Pos_Lmark, BodyCM, ...
%                JointAxes, JointAxesName, PositionGUI)
%
%   Inputs:
%     + BodyName is the name of the body - string
%     + BodyBAF_Pos_Smark is a double array (3 x nSmark) with the 
%       coordinates of the skin-markers located on the body referred
%       to the Body BAF.
%     + BodyBAF_Pos_Lmark is a double array (3 x nSmark) with the 
%       coordinates of the landmarks located on the body referred
%       to the Body BAF.
%     + BodyCM is the double array(3x1) with the coordinates of the
%       centre of mass of the body.
%     + JointAxes is double array(3 x nJointAxes) with the coordinates
%       of the joint axes of the body referred to the body BAF
%     + JointAxesName is a cell (nJointAxes x 1) with the names of the
%       joint axes.
%     + OPTIONAL: PositionGUI is double array(4x1) with the position of the
%       GUI figure that has called this function. This vector contains:
%         PositionGUI(1): X coordinate of the left inferior corner
%         PositionGUI(2): Y coordinate of the left inferior corner
%         PositionGUI(3): Width of the figure
%         PositionGUI(4): Height of the figure
%   Outputs:
%     NONE


% extract input data
BodyName          = varargin{1}; 
BodyBAF_Pos_Smark = varargin{2}; 
BodyBAF_Pos_Lmark = varargin{3}; 
BodyCM            = varargin{4}; 
JointAxes         = varargin{5}; 
JointAxesName     = varargin{6}; 
if nargin == 7
    PositionGUI   = varargin{7};
end

if nargin == 7
    % Position of the GUI
    GUI_PosX   = PositionGUI(1);
    GUI_PosY   = PositionGUI(2);
    GUI_Width  = PositionGUI(3);
    GUI_Height = PositionGUI(4);

    % Dimension of the figure window where this function will draw.
    FigWidth  = 560;
    FigHeight = 420;

    % Create figure in the background
    h = figure('Visible','off','Position',[GUI_Width + 20, GUI_PosY-(FigHeight - GUI_Height), FigWidth, FigHeight]);

    % make figure visible
    figure(h);
else
    figure
end

% figure settings
hold on, axis equal

% Initialize variable
BodyProps = initBodyProps;

% sizes
nSmark = size(BodyBAF_Pos_Smark, 2);
nLmark = size(BodyBAF_Pos_Lmark, 2);

% get index body
IndexBody = getBodyIndex(BodyProps, BodyName);

% calculate figure scale
MaxSmark = max(max(abs(BodyBAF_Pos_Smark)));
MaxLmark = max(max(abs(BodyBAF_Pos_Lmark)));
MaxMarks = [MaxSmark MaxLmark];
ScaleBAF = 1.2 * max(MaxMarks);

% plot axes vectors
quiver3(0,0,0,ScaleBAF,0,0,0,'k'); 
quiver3(0,0,0,0,ScaleBAF,0,0,'k');
quiver3(0,0,0,0,0,ScaleBAF,0,'k');

% plot axes names
text(ScaleBAF,0,0,'X_B_A_F   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % X-axis
text(0,ScaleBAF,0,'Y_B_A_F   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Y-axis
text(0,0,ScaleBAF,'Z_B_A_F   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Z-axis

% Plot calculated skin-markers referred to BAF
for i = 1 : nSmark
    hold on;
    plot3(BodyBAF_Pos_Smark(1,i), BodyBAF_Pos_Smark(2,i), BodyBAF_Pos_Smark(3,i), 'bo');
    text(BodyBAF_Pos_Smark(1,i), BodyBAF_Pos_Smark(2,i), BodyBAF_Pos_Smark(3,i), ...
        ['  M', num2str(BodyProps{IndexBody, 2}(i))], 'FontSize', 10, 'Color', 'b');
    line([BodyBAF_Pos_Smark(1,i),0],[BodyBAF_Pos_Smark(2,i),0],[BodyBAF_Pos_Smark(3,i),0],'Color','b');
end

% Plot landmarks referred to BAF
for i = 1 : nLmark
    hold on;
    plot3(BodyBAF_Pos_Lmark(1,i), BodyBAF_Pos_Lmark(2,i), BodyBAF_Pos_Lmark(3,i), 'rx');
    text(BodyBAF_Pos_Lmark(1,i), BodyBAF_Pos_Lmark(2,i), BodyBAF_Pos_Lmark(3,i), ...
        ['L', num2str(BodyProps{IndexBody, 3}(i)), '  '], 'FontSize', 10, 'FontWeight', 'Bold', 'Color', 'r', ...
        'HorizontalAlignment', 'Right');
    line([BodyBAF_Pos_Lmark(1,i),0],[BodyBAF_Pos_Lmark(2,i),0],[BodyBAF_Pos_Lmark(3,i),0],'Color','r');
end

% Plot Center of Mass
if strcmp(BodyName,'LClav') | strcmp(BodyName,'RClav')
else
    plot3(BodyCM(1),BodyCM(2),BodyCM(3), 'g^');
    text(BodyCM(1),BodyCM(2),BodyCM(3), '  CM', 'FontSize', 10, 'FontWeight', 'Bold', 'Color', 'g');
    line([BodyCM(1),0],[BodyCM(2),0],[BodyCM(3),0],'Color','g');
end

% Plot Joint axes
if isempty(JointAxes)
        
elseif strcmp(BodyName,'RThigh')
    quiver3(BodyBAF_Pos_Lmark(1,4),BodyBAF_Pos_Lmark(2,4),BodyBAF_Pos_Lmark(3,4),JointAxes(1)*ScaleBAF,JointAxes(2)*ScaleBAF,JointAxes(3)*ScaleBAF);
    text((BodyBAF_Pos_Lmark(1,4) + JointAxes(1)*ScaleBAF),(BodyBAF_Pos_Lmark(2,4) + JointAxes(2)*ScaleBAF),(BodyBAF_Pos_Lmark(3,4) + JointAxes(3)*ScaleBAF),JointAxesName,...
        'FontSize',10,'FontWeight','Bold','Color','k','HorizontalAlignment','center');
elseif strcmp(BodyName,'RShank') | strcmp(BodyName,'LFarm') | strcmp(BodyName,'RFarm')
    quiver3(BodyBAF_Pos_Lmark(1,1),BodyBAF_Pos_Lmark(2,1),BodyBAF_Pos_Lmark(3,1),JointAxes(1)*ScaleBAF,JointAxes(2)*ScaleBAF,JointAxes(3)*ScaleBAF);
    text((BodyBAF_Pos_Lmark(1,1) + JointAxes(1)*ScaleBAF),(BodyBAF_Pos_Lmark(2,1) + JointAxes(2)*ScaleBAF),(BodyBAF_Pos_Lmark(3,1) + JointAxes(3)*ScaleBAF),JointAxesName,...
        'FontSize',10,'FontWeight','Bold','Color','k','HorizontalAlignment','center');
elseif strcmp(BodyName,'RUarm') | strcmp(BodyName,'LUarm')
    quiver3(BodyBAF_Pos_Lmark(1,4),BodyBAF_Pos_Lmark(2,4),BodyBAF_Pos_Lmark(3,4),JointAxes(1,1)*ScaleBAF,JointAxes(2,1)*ScaleBAF,JointAxes(3,1)*ScaleBAF);
    text((BodyBAF_Pos_Lmark(1,4) + JointAxes(1,1)*ScaleBAF),(BodyBAF_Pos_Lmark(2,4) + JointAxes(2,1)*ScaleBAF),(BodyBAF_Pos_Lmark(3,4) + JointAxes(3,1)*ScaleBAF),JointAxesName{1},...
        'FontSize',10,'FontWeight','Bold','Color','k','HorizontalAlignment','center');
    quiver3(BodyBAF_Pos_Lmark(1,4),BodyBAF_Pos_Lmark(2,4),BodyBAF_Pos_Lmark(3,4),JointAxes(1,2)*ScaleBAF,JointAxes(2,2)*ScaleBAF,JointAxes(3,2)*ScaleBAF);
    text((BodyBAF_Pos_Lmark(1,4) + JointAxes(1,2)*ScaleBAF),(BodyBAF_Pos_Lmark(2,4) + JointAxes(2,2)*ScaleBAF),(BodyBAF_Pos_Lmark(1,4) + JointAxes(3,2)*ScaleBAF),JointAxesName{2},...
        'FontSize',10,'FontWeight','Bold','Color','k','HorizontalAlignment','center');
end

title(['Local axes of the ', BodyName, '. Skin-markers in blue; Landmarks in red'], 'FontWeight', 'Bold', 'FontSize', 12);

% define points of view depending on the body
if strcmp(BodyName,'Thor')
    view([140.45 -40]);
elseif strcmp(BodyName,'Pel')
    view([119.27 35.17]);
elseif strcmp(BodyName,'LUarm')
    view([111.73 -41.89]);
elseif strcmp(BodyName,'RUarm')
    view([111.73 -41.89]);
elseif strcmp(BodyName,'LFarm')
    view([117.31 -47.66]);
elseif strcmp(BodyName,'RFarm')
    view([117.31 -47.66]);
elseif strcmp(BodyName,'LHand')
    view([126.9 -55.27]);
elseif strcmp(BodyName,'RHand')
    view([126.9 -55.27]);
elseif strcmp(BodyName,'RThigh')
    view([105.51 38.96]);
elseif strcmp(BodyName,'RShank')
    view([99.51 18.76]);
elseif strcmp(BodyName,'RFoot')
    view([111.24 58.22]);
end
