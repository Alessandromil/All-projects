function Marker = getMarkerInSeg(Segments,SegName,MarkerName)

%GETMARKERINSEG seek marker in the local markert list in the Segment list.
% and return a handle to this marker.

Marker = [];
NSegments = size(Segments,1);

for i= 1:NSegments
    if strcmpi(SegName,Segments(i).Name)
        NPoints = size(Segments(i).LocalMarkers,1);
        for j = 1:NPoints
            if strcmpi(MarkerName,Segments(i).LocalMarkers(j).Name)
                Marker = Segments(i).LocalMarkers(j);
                return;
            end
        end
    end
end



%error('Segment "%s" is not defined for this model or Point "%s" does not belong to the segment',SegName,PointName)


