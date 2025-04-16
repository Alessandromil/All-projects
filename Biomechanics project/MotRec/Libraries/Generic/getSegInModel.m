function SegIndex = getSegInModel(Segments,SegName)

%GETSEGINMODEL seek Segment in the Segment list & return the segment index

SegIndex = [];
NSegments = size(Segments,1);

for i= 1:NSegments
    if strcmpi(SegName,Segments(i).Name)
        SegIndex = i;
    end
end
