function Point = getPointInSeg(Segments,SegName,PointName)
%GETPOINTINSEG seek point in the local point list in the Segment list.
% and return a handle to this point.

Point = [];
NSegments = size(Segments,1);

for i= 1:NSegments
    if strcmpi(SegName,Segments(i).Name)
        NPoints = size(Segments(i).LocalPoints,1);
        for j = 1:NPoints
            if strcmpi(PointName,Segments(i).LocalPoints(j).Point.Name)
                Point = Segments(i).LocalPoints(j);
                return;
            end
        end
    end
end

%error('Segment "%s" is not defined for this model or Point "%s" does not belong to the segment',SegName,PointName)


