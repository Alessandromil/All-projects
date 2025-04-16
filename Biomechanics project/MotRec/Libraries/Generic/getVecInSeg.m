function Vector = getVecInSeg(Segments,SegName,VectorName)

% GETVECINSEG seek vector in the Local point list in SegmentList of Human.
% and return a handle to this vector.

Vector = [];
NSegments = size(Segments,1);

for i= 1:NSegments
    if strcmpi(SegName,Segments(i).Name)
        NVectors = size(Segments(i).LocalVectors,1);
        for j = 1:NVectors
            if strcmpi(VectorName,Segments(i).LocalVectors(j).Vector.Name)
                Vector = Segments(i).LocalVectors(j);
                return;
            end
        end
    end
end

%error('Vector %s does not belong to Segment %s or Segment %s does not exist',VectorName,SegName,SegName)
