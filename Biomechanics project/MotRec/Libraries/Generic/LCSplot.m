function h=LCSplot(ax,X,Y,Z,SegmentName)

h=gobjects(3,2);

% find scale
MaximumCoords=zeros(1,3);
MaximumCoords(1,1)  = max(max(abs(X(:,:))));
MaximumCoords(1,2)  = max(max(abs(Y(:,:))));
MaximumCoords(1,3)  = max(max(abs(Z(:,:))));
MaxCoord = max(MaximumCoords);
LCSlength = 1.3 * MaxCoord;

% plot axes vectors
h(1,1)=quiver3(ax,0,0,0,LCSlength,0,0,0,'r','LineWidth',2);
h(2,1)=quiver3(ax,0,0,0,0,LCSlength,0,0,'g','LineWidth',2);
h(3,1)=quiver3(ax,0,0,0,0,0,LCSlength,0,'b','LineWidth',2);

% plot axes names
axisname=strings(1,3);
if strcmp(SegmentName,'Ground')==1 %Enters here if it is drawing the GCS
    axisname(1,1)=strcat('X_{',SegmentName,'}   ');
    axisname(1,2)=strcat('Y_{',SegmentName,'}   ');
    axisname(1,3)=strcat('Z_{',SegmentName,'}   ');
    h(1,2)=text(ax,LCSlength,0,0,axisname(1,1),'FontSize',12,'FontWeight','Bold','Color','r','HorizontalAlignment','right','BackgroundColor','none'); % X-axis
    h(2,2)=text(ax,0,LCSlength,0,axisname(1,2),'FontSize',12,'FontWeight','Bold','Color','g','HorizontalAlignment','right','BackgroundColor','none'); % Y-axis
    h(3,2)=text(ax,0,0,LCSlength,axisname(1,3),'FontSize',12,'FontWeight','Bold','Color','b','HorizontalAlignment','right','BackgroundColor','none'); % Z-axis
else %Enters here if it is drawing a LCS
    axisname(1,1)=strcat('X_{',SegmentName,'}   ');
    axisname(1,2)=strcat('Y_{',SegmentName,'}   ');
    axisname(1,3)=strcat('Z_{',SegmentName,'}   ');
    h(1,2)=text(ax,LCSlength,0,0,axisname(1,1),'FontSize',10,'FontWeight','normal','Color','r','HorizontalAlignment','right','BackgroundColor','none'); % X-axis
    h(2,2)=text(ax,0,LCSlength,0,axisname(1,2),'FontSize',10,'FontWeight','normal','Color','g','HorizontalAlignment','right','BackgroundColor','none'); % Y-axis
    h(3,2)=text(ax,0,0,LCSlength,axisname(1,3),'FontSize',10,'FontWeight','normal','Color','b','HorizontalAlignment','right','BackgroundColor','none'); % Z-axis
end