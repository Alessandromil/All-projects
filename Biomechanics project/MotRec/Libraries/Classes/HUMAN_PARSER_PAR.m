classdef HUMAN_PARSER_PAR < handle
    %HUMAN_PARSER_PAR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        HumanDOM        % Document Object Model node representing the parsed document
    end
    
    methods
        function HPPar = HUMAN_PARSER_PAR(CalibPath,CalibFile)
            % Constructor of class
            
            HPPar.HumanDOM = xmlread([CalibPath,CalibFile]);
        end
        function parseMarkers(HPPar,Human)
            HumanSegments = HPPar.HumanDOM.getElementsByTagName('Segment');
            for i = 1:HumanSegments.getLength
                % Add Segmets
                SegmentListItem = HumanSegments.item(i-1);
                SegmentName     = char(SegmentListItem.getAttribute('Name'));
                SegmentPoints   = SegmentListItem.getElementsByTagName('Point');
                SegmentMarkers  = SegmentListItem.getElementsByTagName('Marker');
                SegmentIndex = getVecIndex(SegmentName,Human.Segments);
                if SegmentIndex==0
                    error('The segment "%s" does not belong to the model',SegmentName)
                end
                for j = 1:SegmentPoints.getLength
                    PointsListItem = SegmentPoints.item(j-1);
                    PointName = char(PointsListItem.getAttribute('Name'));
                    PointLCoord = char(PointsListItem.getAttribute('LocCoord'));
                    PointIndex = getLocVecIndex(PointName,Human.Segments(SegmentIndex).LocalPoints);
                    if isempty(Human.Segments(SegmentIndex).LocalPoints(PointIndex).LocCoord)
                        Human.Segments(SegmentIndex).LocalPoints(PointIndex).LocCoord = str2num(PointLCoord)*0.001;
                    end
                end
                for j = 1:SegmentMarkers.getLength
                    MarkerListItem = SegmentMarkers.item(j-1);
                    MarkerName = char(MarkerListItem.getAttribute('Name'));
                    MarkerLCoord = char(MarkerListItem.getAttribute('LocCoord'));
                    if ~isempty(MarkerLCoord)
                        NumMarkerCoord = str2num(MarkerLCoord)*0.001;
                    end
                    % add marker in list of markers in human and asign its local coordinate value
                    AddedMarker = Human.addMarker(MarkerName);
                    
                    %  Add Local Point to Segment                                                  
                    Human.Segments(SegmentIndex).addMarker(AddedMarker);
                    NSegmentMarkers = size(Human.Segments(SegmentIndex).LocalMarkers,1);
                    Human.Segments(SegmentIndex).LocalMarkers(NSegmentMarkers).LocCoord = NumMarkerCoord;
                    
                end
            end
            % add markers in q
            Human.addMarkerInq();
        end

    end
    
end

