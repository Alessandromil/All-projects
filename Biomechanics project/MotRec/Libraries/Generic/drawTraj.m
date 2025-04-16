function drawTraj(Markers, varargin)

MarkerNames = fieldnames(Markers);

if nargin == 2
    Markers_Filt = varargin{1};
end

nMarkers = length(MarkerNames);

FontSize  = 14;
LineWidth =  2;
for i=1:nMarkers
    if nargin == 1
        figure('Name',['Marker ',MarkerNames{i}], 'NumberTitle','off','WindowStyle','docked');
    elseif nargin == 2
        figure('Name',['Filt Marker ',MarkerNames{i}], 'NumberTitle','off','WindowStyle','docked');
    end
    axes('FontSize',FontSize), hold on;
    xlabel('Numero de frame','FontSize',FontSize);
    ylabel('Coordenada [m]', 'FontSize',FontSize);
    if nargin == 1
        plot(Markers.(MarkerNames{i}), 'LineWidth',LineWidth);  
        legend('X coord','Y coord','Z coord');        
    elseif nargin == 2
        plot(Markers.(MarkerNames{i}), 'LineWidth',LineWidth);  
        plot(Markers_Filt.(MarkerNames{i}), ':', 'LineWidth',LineWidth);  
        legend('orig X coord','orig Y coord','orig Z coord',...
               'filt X coord','filt Y coord','filt Z coord');        
    end
end    

