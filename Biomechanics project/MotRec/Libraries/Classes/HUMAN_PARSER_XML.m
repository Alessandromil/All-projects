classdef HUMAN_PARSER_XML < handle
    %HumanParseXML read XML human definition and fills HumanSym class.
    %   
    
    
    properties
        HumanDOM        % Document Object Model node representing the parsed document
    end
    
    methods
        function HPxml = HUMAN_PARSER_XML()
            % Constructor of class
        end
        function [OK, Message] = checkInvalidChars(HPxml, String, InvalidChars)
            Message = [];
            OK = 1;
            nChars = length(InvalidChars);
            for i=1:nChars
                InvalidChar = InvalidChars(i);
                if ~isempty(strfind(String, InvalidChar))
                    if strcmp(' ',InvalidChar)
                        InvalidCharText = '"blank space"';                        
                    else
                        InvalidCharText = InvalidChar;
                    end
                    Message = [Message, sprintf(['     String contains invalid character ',InvalidCharText,'\n'])];
                    OK = 0;                    
                end                
            end
        end
        function isValidDrawSeq(HPxml, SegmentName, DrawSeq)
            [OK, Message] = HPxml.checkInvalidChars(DrawSeq,';[] ');
            if ~OK
                Header = sprintf(['Segment "',SegmentName,'" -> "DrawSeq" contains invalid character(s):\n']);
                error([Header,Message]);
            end
        end
        function [OK, Message] = isStrArray3x1(HPxml, CharArray)
            % checks if string CharArray is a valid 3x1 vector:
            % - if Yes OK=1, if No OK=0
            % - spaces at the beginning or end are not a problem for function
            %   str2num, which will transform the string into numeric array
            % -
            CharArray = strtrim(CharArray); % Remove leading and trailing white space from string
            % check start and end brakets
            FirstChar = CharArray(1); % must be [
            EndChar   = CharArray(end); % must be ]
            Message = [];
            OK = 1;
            if ~strcmp(FirstChar,'[')
                Message = [Message, sprintf('  Vector must start with symbol [\n')];
                OK = 0;
            end
            if ~strcmp(EndChar,']')
                Message = [Message, sprintf('  Vector must end with symbol ]\n')];
                OK = 0;
            end
            % check number of characters ; (there must be 2)
            if length(strfind(CharArray(2:end-1),';')) ~= 2
                Message = [Message, sprintf('  Each coordinate must be separated by ;\n')];
                OK = 0;
            end
            % check decimal character and separation character between coords
            if ~isempty(strfind(CharArray(2:end-1),','))
                Message = [Message, sprintf('  Character , is not allowed. Decimal symbol is . and coordinates must be separated by ;\n')];
                OK = 0;
            end
            % check number of coordinates
            NumArray = str2num(CharArray);
            if length(NumArray) ~= 3
                Message = [Message, sprintf('  Three coordinates must be defined\n')];
                OK = 0;
            end
        end
        function isValidCoMCoord(HPxml, SegmentName, CoMLCoord)
            % checks if string CoMLCoord are valid CoM coordinates
            [OK, Message] = HPxml.isStrArray3x1(CoMLCoord);
            if ~OK
                Header = sprintf(['Segment ',SegmentName,' -> CoM coordinates have not a valid format:\n']);
                error([Header,Message]);
            end            
        end
        function isValidGraphicTras(SegmentName,GraphicTras)
            % checks if string GraphicTras is a valid translation
            [OK, Message] = HPxml.isStrArray3x1(GraphicTras);
            if ~OK
                Header = sprintf(['Segment ',SegmentName,' -> Graphic definition -> translation has not a valid format:\n']);
                error([Header,Message]);
            end            
        end        
        function isValidMarkerCoord(HPxml, SegmentName, MarkerName, MarkerLCoord)
            % checks if string MarkerLCoord are valid marker coordinates
            [OK, Message] = HPxml.isStrArray3x1(MarkerLCoord);
            if ~OK
                Header = sprintf(['Segment ',SegmentName,' -> Marker ',MarkerName,' -> coordinates have not a valid format:\n']);
                error([Header,Message]);
            end
        end
        function isValidVectorCoord(HPxml, SegmentName, VectorName, VectorLCoord)
            % checks if string VectorLCoord are valid vector coordinates
            [OK, Message] = HPxml.isStrArray3x1(VectorLCoord);
            if ~OK
                Header = sprintf(['Segment ',SegmentName,' -> Vector ',VectorName,' -> coordinates have not a valid format:\n']);
                error([Header,Message]);
            end
        end
        
        function isValidPointCoord(HPxml, SegmentName, PointName, PointCoord)
            % checks if string PointCoord are valid point coordinates
            [OK, Message] = HPxml.isStrArray3x1(PointCoord);
            if ~OK
                Header = sprintf(['Segment ',SegmentName,' -> Point ',PointName,' -> coordinates have not a valid format:\n']);
                error([Header,Message]);
            end
        end
        
        function isValidPointName(HPxml,PointName)
            if strcmpi(PointName,'A')
                error(['You have choose "A" as point name in segment ',SegmentName,'\n.       You must choose other name because the A name is incompatible with the symbolic toolbox'] )
            end            
        end
        
        function parseAngles(HPxml,HumanModel)
            HumanAngles = HPxml.HumanDOM.getElementsByTagName('Angle');
            for i = 1:HumanAngles.getLength
                AngleListItem = HumanAngles.item(i-1);
                AngleName1 = char(AngleListItem.getAttribute('Name1'));
                AngleName2 = char(AngleListItem.getAttribute('Name2'));
                AngleName3 = char(AngleListItem.getAttribute('Name3'));
                AngleJoint = char(AngleListItem.getAttribute('Joint'));
                AngleVec1Seg1 = char(AngleListItem.getAttribute('Seg1VecY'));
                AngleVec2Seg1 = char(AngleListItem.getAttribute('Seg1VecZ'));
                AngleVec1Seg2 = char(AngleListItem.getAttribute('Seg2VecY'));
                AngleVec2Seg2 = char(AngleListItem.getAttribute('Seg2VecZ'));
                AngleSeg1VecRef = char(AngleListItem.getAttribute('Seg1VecRef'));
                AngleSeg2VecRef = char(AngleListItem.getAttribute('Seg2VecRef'));
                AngleRotSeq = char(AngleListItem.getAttribute('RotSeq'));
                HumanModel.addAngle(AngleName1,AngleName2,AngleName3,AngleJoint,AngleVec1Seg1,AngleVec1Seg2,...
                                    AngleVec2Seg1,AngleVec2Seg2,AngleSeg1VecRef,AngleSeg2VecRef,AngleRotSeq);
            end
        end
        
        function parseJoints(HPxml,HumanModel)
            % Get joint of XML and set to Model
            HumanJoints = HPxml.HumanDOM.getElementsByTagName('Joint');
            for i = 1:HumanJoints.getLength
                JointListItem = HumanJoints.item(i-1);
                JointName = char(JointListItem.getAttribute('Name'));
                JointType = char(JointListItem.getAttribute('Type'));
                JointSeg1 = char(JointListItem.getAttribute('Seg1'));
                JointSeg2 = char(JointListItem.getAttribute('Seg2'));
                JointPoint1 = char(JointListItem.getAttribute('Point1'));
                JointPoint2 = char(JointListItem.getAttribute('Point2'));
                JointVecSeg1 = char(JointListItem.getAttribute('Seg1Axis'));
                JointVecSeg2 = char(JointListItem.getAttribute('Seg2Axis'));
                V1_Ang_V2    = char(JointListItem.getAttribute('AxesAng'));
                
                % Check if common Attributes for all joints are defined
                if isempty(JointName), error('There is a "Joint" element without attribute "Name"'); end
                if isempty(JointType),
                    error(['Joint "',JointName,'" -> Attribute "Type" not defined. Options: FLOAT, REV, UNI, SPH'])
                end

                % Attributes specific of each joint are checked inside function addJoint
                HumanModel.addJoint(JointName,JointType,JointSeg1,JointSeg2,JointPoint1,JointPoint2,JointVecSeg1,JointVecSeg2,V1_Ang_V2);
            end
        end
        
        function parseName(HPxml,HumanModel)%,FileXML
           % Get name of XMl,Coord type and Model Type and set to Model
           NHumans = HPxml.HumanDOM.getElementsByTagName('HumanModel');
           HumanModel.ModelName = char(NHumans.item(0).getAttribute('Name'));
           HumanModel.CoordType = char(NHumans.item(0).getAttribute('CoordType'));
           HumanModel.ModelType = char(NHumans.item(0).getAttribute('ModelType'));
           if ~strcmp(HumanModel.CoordType,'Natural')
               error('Attribute "CoordType" must be "Natural"')
           end
           if isempty(HumanModel.ModelName)
               error('Attribute "Name" must be provided for element "HumanModel"')
           end
%            if ~strcmp(HumanModel.ModelName,FileXML(1:(end-4)))
%                error('The name of the xml file and the model name defined inside the file must be the same')
%            end
        end
        
        function parseSegments(HPxml,HumanModel,ModelPath)
            global PathBar PathBar2
            % Get segment and his point and vectors of XML and set to Model
            HumanSegments = HPxml.HumanDOM.getElementsByTagName('Segment');
            FixedBodies = 0;
            
            for i = 1:HumanSegments.getLength
                
                % Initialize vars for each segment
                NPoints_withLCoord = 0; % For current segment total number of points
                                        % for which the local coordinates are provided.
                NVectors_withLCoord = 0; % For current segment total number of vectors
                                         % for which the local coordinates are provided.
                
                % Add Segmets
                SegmentListItem = HumanSegments.item(i-1);
                SegmentName = char(SegmentListItem.getAttribute('Name'));
                HumanModel.addSegment(SegmentName);
                SegmentMass = char(SegmentListItem.getAttribute('Mass'));
                HumanModel.Segments(i).Mass = str2num(SegmentMass);
                
                % Check attribute Type="Fixed" if provided
                SegmentType = char(SegmentListItem.getAttribute('Type'));
                if strcmp(SegmentType,'Fixed')
                    FixedBodies = FixedBodies + 1;
                    if FixedBodies > 1
                        error(['Segment "',SegmentName,'" -> Type="Fixed". There is already one "Fixed" segment']);
                    end
                    % Set segment as Fixed
                    HumanModel.Segments(i).Fixed = 1;
                elseif ~isempty(SegmentType)
                    error(['Segment "',SegmentName,'" -> attribute "Type" can take only one value: "Fixed"']);
                end
                
                % Calc number of points and vector in segments
                SegmentPoints    = SegmentListItem.getElementsByTagName('Point');
                SegmentMarkers   = SegmentListItem.getElementsByTagName('Marker');
                SegmentPointsRel = SegmentListItem.getElementsByTagName('PointRel');
                SegmentCoM       = SegmentListItem.getElementsByTagName('CoM');
                SegmentI         = SegmentListItem.getElementsByTagName('MoI');
                SegmentVectors   = SegmentListItem.getElementsByTagName('Vector');
                SegmentGraphics  = SegmentListItem.getElementsByTagName('Graphic');
                
                for j = 1:SegmentPoints.getLength
                    % Add points and his propertis in Human and in Segment
                    PointListItem = SegmentPoints.item(j-1);
                    PointName = char(PointListItem.getAttribute('Name'));
                    PointLCoord = char(PointListItem.getAttribute('LocCoord'));
                    
                    % Check if all Attributes are provided
                    if isempty(PointName)
                        error(['Segment "',SegmentName,'" -> one "Point" element has not attribute "Name"']);
                    end
                    if isempty(PointLCoord)
                        error(['Segment "',SegmentName,'" -> Point "',PointName,'" -> attribute "LocCoord" is not defined']);
                    end
                    % Attribute LocCoord is optional in DHErgo because local coordinates 
                    % are provided later with the subject paremeter file. 
                    % However, this version is for Biomechanics and includes this check !!!!!
                    
                    % check if point belongs already to the segment
                    HPxml.isValidPointName(PointName);
                    PointIndex = getPointInSeg(HumanModel.Segments,SegmentName,PointName); 
                    if ~isempty(PointIndex)
                        error(['Segment "',SegmentName,'" -> Point "',PointName,'" -> this point is defined twice']);
                    end
                    
                    % add point to Point list and Segment                    
                    AddedPoint = HumanModel.addPoint(PointName);
                    HumanModel.Segments(i).addPoint(AddedPoint);
                    
                    % If point is in a fixed segment it's fixed
                    if HumanModel.Segments(i).Fixed == 1
                        AddedPoint.Fixed = 1;
                    end
                    
                    % Add local coord to point
                    if ~isempty(PointLCoord)
                        HPxml.isValidPointCoord(SegmentName,PointName,PointLCoord)
                        NumPointCoord = str2num(PointLCoord);
                        HumanModel.Segments(i).LocalPoints(j).LocCoord = NumPointCoord;
                        NPoints_withLCoord = NPoints_withLCoord + 1;
                    end
                    
                    % Add global coord if the segmetn is of Type=Fixed
                    if HumanModel.Segments(i).Fixed == 1
                    %if strcmpi(HumanModel.Segments(i).Name,'Ground')
                        HumanModel.Segments(i).LocalPoints(j).Point.GlobalCoord = ...
                            HumanModel.Segments(i).LocalPoints(j).LocCoord;
                    end
                end
                
                for j = 1:SegmentCoM.getLength
                    % Add the CoM in Human and in Segment
                    if j>1
                        error(['Segment ',SegmentName,' -> There is more than 1 CoM in the segment']);
                    end
                    CoMListItem = SegmentCoM.item(j-1);
                    AddedCoM = HumanModel.addCoM(SegmentName);
                    HumanModel.Segments(i).CoM = LOCAL_POINT(AddedCoM,SegmentName);
                    % Add local coord to CoM
                    CoMLCoord = char(CoMListItem.getAttribute('LocCoord'));
                     if ~isempty(CoMLCoord)
                        HPxml.isValidCoMCoord(SegmentName,CoMLCoord)
                        NumCoMCoord = str2num(CoMLCoord);
                        HumanModel.Segments(i).CoM(j).LocCoord = NumCoMCoord;
                    end
                end
                
                for j = 1:SegmentI.getLength
                    % Add the MoI in Human and in Segment
                    if j>1
                        error(['Segment ',SegmentName,' -> There is more than 1 I in the segment']);
                    end
                    MoIListItem = SegmentI.item(j-1);
                    % Add local coord to MoI
                    Ixx = str2num(char(MoIListItem.getAttribute('Ixx')));
                    Ixy = str2num(char(MoIListItem.getAttribute('Ixy')));
                    Ixz = str2num(char(MoIListItem.getAttribute('Ixz')));
                    Iyx = str2num(char(MoIListItem.getAttribute('Iyx')));
                    Iyy = str2num(char(MoIListItem.getAttribute('Iyy')));
                    Iyz = str2num(char(MoIListItem.getAttribute('Iyz')));
                    Izx = str2num(char(MoIListItem.getAttribute('Izx')));
                    Izy = str2num(char(MoIListItem.getAttribute('Izy')));
                    Izz = str2num(char(MoIListItem.getAttribute('Izz')));
                    if isempty(Ixx), error(['Segment ',SegmentName,' -> Ixx has a not valid value']); end
                    if isempty(Ixy), error(['Segment ',SegmentName,' -> Ixy has a not valid value']); end
                    if isempty(Ixz), error(['Segment ',SegmentName,' -> Ixz has a not valid value']); end
                    if isempty(Iyx), error(['Segment ',SegmentName,' -> Iyx has a not valid value']); end
                    if isempty(Iyy), error(['Segment ',SegmentName,' -> Iyy has a not valid value']); end
                    if isempty(Iyz), error(['Segment ',SegmentName,' -> Iyz has a not valid value']); end
                    if isempty(Izx), error(['Segment ',SegmentName,' -> Izx has a not valid value']); end
                    if isempty(Izy), error(['Segment ',SegmentName,' -> Izy has a not valid value']); end
                    if isempty(Izz), error(['Segment ',SegmentName,' -> Izz has a not valid value']); end
                    HumanModel.Segments(i).I = [Ixx,Ixy,Ixz;Iyx,Iyy,Iyz;Izx,Izy,Izz];
                end
                
                for j = 1:SegmentMarkers.getLength
                    % Add markers and his propertis in Human and in Segment
                    MarkerListItem = SegmentMarkers.item(j-1);
                    MarkerName     = char(MarkerListItem.getAttribute('Name'));
                    MarkerLCoord   = char(MarkerListItem.getAttribute('LocCoord'));
                    
                    % Check if all Attributes are provided
                    if isempty(MarkerName)
                        error(['Segment "',SegmentName,'" -> one "Marker" element has not attribute "Name"']);
                    end
                    if isempty(MarkerLCoord)
                        error(['Segment "',SegmentName,'" -> Marker "',MarkerName,'" -> attribute "LocCoord" is not defined']);
                    end
                    % Attribute LocCoord is optional in DHErgo because local coordinates 
                    % are provided later with the subject paremeter file. 
                    % However, this version is for Biomechanics and includes this check !!!!!
                    
                    % Attribute LocCoord is optional. In DHErgo coordinates are provided later.
                    
                    % check if Marker belongs already to the segment
                    MarkerIndex = getMarkerInSeg(HumanModel.Segments,SegmentName,MarkerName);
                    if ~isempty(MarkerIndex)
                        error(['Segment "',SegmentName,'" -> Marker "',MarkerName,'" -> this marker is defined twice']);
                    end

                    % add marker to Marker list and Segment                    
                    AddedMarker = HumanModel.addMarker(MarkerName);
                    HumanModel.Segments(i).addMarker(AddedMarker);
                    
                    % Add local coord to marker
                    if ~isempty(MarkerLCoord)
                        HPxml.isValidMarkerCoord(SegmentName,MarkerName,MarkerLCoord)
                        NumMarkerCoord = str2num(MarkerLCoord);
                        HumanModel.Segments(i).LocalMarkers(j).LocCoord = NumMarkerCoord;
                    end
                end
                
                for j = 1:SegmentPointsRel.getLength
                    PointRelItem = SegmentPointsRel.item(j-1);
                    PointRelName = char(PointRelItem.getAttribute('Name'));
                    AddedPointRel = HumanModel.addRelPoint(PointRelName);
                    HumanModel.Segments(i).addRelPoint(AddedPointRel);                    
                end
                
                for j = 1:SegmentVectors.getLength

                    % Add vectors and his properties in Human and in Segment
                    VectorListItem = SegmentVectors.item(j-1);
                    VectorName = char(VectorListItem.getAttribute('Name'));
                    VectorLCoord = char(VectorListItem.getAttribute('LocCoord'));

                    % Check if all Attributes are provided
                    if isempty(VectorName)
                        error(['Segment "',SegmentName,'" -> one "Vector" element has not attribute "Name"']);
                    end
                    if isempty(VectorLCoord)
                        error(['Segment "',SegmentName,'" -> Vector "',VectorName,'" -> attribute "LocCoord" is not defined']);
                    end
                    % Attribute LocCoord is optional in DHErgo because local coordinates 
                    % are provided later with the subject paremeter file. 
                    % However, this version is for Biomechanics and includes this check !!!!!
                   
                    % check if vector belongs already to the segment
                    VectorIndex = getVecInSeg(HumanModel.Segments,SegmentName,VectorName);
                    if ~isempty(VectorIndex)
                        error(['Segment "',SegmentName,'" -> Vector "',VectorName,'" -> this vector is defined twice']);
                    end

                    % add vector to Vector list and Segment                    
                    AddedVector = HumanModel.addVector(VectorName);
                    HumanModel.Segments(i).addVector(AddedVector);
                    
                    % If vecttor is in fixed segment is fixed
                    if HumanModel.Segments(i).Fixed == 1
                        AddedVector.Fixed = 1;
                    end
                    
                    % Add local coord to vector
                    if ~isempty(VectorLCoord)
                        HPxml.isValidVectorCoord(SegmentName,VectorName,VectorLCoord)
                        NumVectorCoord = str2num(VectorLCoord);
                        HumanModel.Segments(i).LocalVectors(j).LocCoord = NumVectorCoord;
                        NVectors_withLCoord = NVectors_withLCoord + 1;
                    end
                    
                    % Add global coord if is Ground
                    if HumanModel.Segments(i).Fixed == 1 % If is Ground
                    %if strcmpi(HumanModel.Segments(i).Name,'Ground')
                        HumanModel.Segments(i).LocalVectors(j).Vector.GlobalCoord = ...
                            HumanModel.Segments(i).LocalVectors(j).LocCoord;
                    end
                end
                
                for j = 1:SegmentGraphics.getLength
                    
                    GraphicListItem = SegmentGraphics.item(j-1);
                    
                    % DrawSeq is a list of points belonging to the Segment. Format: DrawSeq="P1,P2,P3,P1"
                    DrawSeq = char(GraphicListItem.getAttribute('DrawSeq'));
                    % GraphicFile is a .x file with graphic objects
                    GraphicFile = char(GraphicListItem.getAttribute('File'));                    
                    % Check if at least one of them (DrawSeq or File) has been defined. If not give error
                    if isempty(GraphicFile) && isempty(DrawSeq)
                        error(['Segment "',SegmentName,'" -> "Graphic" element defined but not "DrawSeq" or "File" attribute found'])
                    end
                    
                    % ----------------------------------------------------
                    % If a wireframe is defined by DrawSeq
                    % ----------------------------------------------------
                    if ~isempty(DrawSeq) % add wireframe to Segment 
                        
                        Radius = char(GraphicListItem.getAttribute('Radius'));
                        Color = char(GraphicListItem.getAttribute('Colour'));
                        
                        % Check if all Attributes are provided
                        if isempty(Radius)
                            error(['Segment "',SegmentName,'" -> "Graphic DrawSeq" -> Attribute "Radius" not defined.'])
                        end
                        if isempty(Color)
                            error(['Segment "',SegmentName,'" -> "Graphic DrawSeq" -> Attribute "Colour" not defined.'])
                        end                        
                        
                        % Check if DrawSeg has the right format e.g. DrawSeq="P1,P2,P3,P1"
                        DrawSeq = strtrim(DrawSeq); % Remove leading and trailing white space from string
                        HPxml.isValidDrawSeq(SegmentName, DrawSeq);
                        % Extract points in DrawSeq
                        PointsInDrawSeq = strread(DrawSeq,'%s','delimiter',',');
                        nPointsInDrawSeq = length(PointsInDrawSeq);
                        
                        % Check if Points or Markers in DrawSeq belong to Segment
                        PointsInSeg  = HumanModel.Segments(i).LocalPoints;
                        MarkersInSeg = HumanModel.Segments(i).LocalMarkers;
                        for k=1:nPointsInDrawSeq
                            PointIndex = getLocVecIndex(PointsInDrawSeq{k}, PointsInSeg);
                            MarkerIndex = getLocVecIndex(PointsInDrawSeq{k}, MarkersInSeg);
                            if PointIndex == 0 && MarkerIndex == 0
                                error(['Segment "',SegmentName,'" -> "DrawSeq" contains Point or Marker name "',PointsInDrawSeq{k},'" not defined in segment'])
                            end
                        end
                        
                        % add graphic 2 model
                        HumanModel.Segments(i).setGraphicWireFrame(DrawSeq,Radius,Color);
                    end
                    
                    % ----------------------------------------------------
                    % If a graphic x-file is associated to the segment
                    % ----------------------------------------------------
                    if ~isempty(GraphicFile) % add graphic file to Segment
                        % check if graphic file "GraphiFile" exists and has the right extension
                        checkFileAndPath([ModelPath,'graphics',PathBar], GraphicFile, {'x';'stl'});    %<===== Ignacio
                        % get Tras, Rot, Scal
                        GraphicTras = char(GraphicListItem.getAttribute('Tras')); % symbolic coords are substituted by their values in getSubject(EXPERIMENT)
                        GraphicRot  = char(GraphicListItem.getAttribute('Rot'));
                        GraphicScal = char(GraphicListItem.getAttribute('Scal'));
                        
                        % Check if all Attributes are provided
                        if isempty(GraphicTras)
                            error(['Segment "',SegmentName,'" -> "Graphic File" -> Attribute "Tras" not defined.'])
                        end
                        if isempty(GraphicRot)
                            error(['Segment "',SegmentName,'" -> "Graphic File" -> Attribute "Rot" not defined.'])
                        end
                        if isempty(GraphicScal)
                            error(['Segment "',SegmentName,'" -> "Graphic File" -> Attribute "Scal" not defined.'])
                        end
                        
                        % check if                        
                        NumGraphic{1} = HumanModel.Segments(i).setGraphicCell(GraphicTras);
                        NumGraphic{2} = HumanModel.Segments(i).setGraphicCell(GraphicRot); 
                        NumGraphic{3} = HumanModel.Segments(i).setGraphicCell(GraphicScal); 
                        
                        % check if all point names of vectors included in NumGraphic{?} have been
                        % replaced. If there is some error in the point or vector name definition
                        % the whole name or part of it is left in NumGraphic{?}
                        Attrib = {'Tras';'Rot';'Scal'};
                        for k=1:3
                            if sum(isletter(NumGraphic{k})) > 0 % there are letters left in "NumGraphic"
                                error(['Segment "',SegmentName,'" -> Graphic File="',GraphicFile,...
                                    '" -> Attribute "',Attrib{k},'" ->\n',...
                                    '    It contains a symbolic point or vector incorrectly defined\n',...
                                    '    Check spelling (Upper case, lower case, etc.)']);
                            end
                        end
                        
                        
                        % add graphic file to Segment
                        HumanModel.Segments(i).setGraphicFiles(GraphicFile,GraphicTras,GraphicRot,GraphicScal);
                    end
                    
                end
                
                % ----------------------------------------------------------------
                % Check if the segment base is properly defined.
                %  THIS CAN BE DONE ONLY WHEN LocCoord ARE DEFINED for all points
                %  and vectors of each segment. (Not possible in DHErgo)
                % ----------------------------------------------------------------
                % Check if the number of points & vectors is correct to define a base.
                % Two options are possible:
                %     1) Three vectors(not aligned) and one point
                %     2) Two points and two vectors(defining three directions not aligned)
                % Then error should be given when:
                %     1) The segment has ONLY: 2 points & 1 vector
                %                              OR 1 point & 2 vectors
                %     2) The directions define by the points and vectors are very close
                %        to each other.

                % Only if the local coordinates of all points and vectors of 
                % the current segment are defined, the base is checked
                if (NPoints_withLCoord == SegmentPoints.getLength) && ...
                        (NVectors_withLCoord == SegmentVectors.getLength)
                    % 1) Check if the number of points and vectors is correct
                    if (SegmentPoints.getLength <= 2 && SegmentVectors.getLength <= 1) || ...
                       (SegmentPoints.getLength <= 1 && SegmentVectors.getLength <= 2)
                        error(['Segment "',SegmentName,'" -> The ',num2str(NPoints_withLCoord),...
                            ' points & ',num2str(NVectors_withLCoord),' vectors provided are not enough to define the segment.\n' ...
                            '        At least THREE vectors and ONE point OR TWO vectors and TWO points are need.']);
                    end
                    
                    % 2) Check if directions defined by points an vectors are not aligned
                    if (SegmentPoints.getLength >= 2 && SegmentVectors.getLength >= 2)
                        P1   = HumanModel.Segments(i).LocalPoints(1).LocCoord;
                        P2   = HumanModel.Segments(i).LocalPoints(2).LocCoord;
                        Dir1 = P2-P1;
                        Dir2 = HumanModel.Segments(i).LocalVectors(1).LocCoord;
                        Dir3 = HumanModel.Segments(i).LocalVectors(2).LocCoord;
                    elseif (SegmentPoints.getLength == 1 && SegmentVectors.getLength >= 3)
                        Dir1 = HumanModel.Segments(i).LocalVectors(1).LocCoord;
                        Dir2 = HumanModel.Segments(i).LocalVectors(2).LocCoord;
                        Dir3 = HumanModel.Segments(i).LocalVectors(3).LocCoord;
                    end
                    
                    % Volumen defined by the three directions defining the base
                    %  vol = abs( dot( cross(Dir1,Dir2),Dir3 ) );
                    vol = abs( (Dir1(2)*Dir2(3)-Dir1(3)*Dir2(2)) * Dir3(1) +  ...
                        (Dir1(3)*Dir2(1)-Dir1(1)*Dir2(3)) * Dir3(2) + ...
                        (Dir1(1)*Dir2(2)-Dir1(2)*Dir2(1)) * Dir3(3)  );
                    par = vol/(norm(Dir1)*norm(Dir2)*norm(Dir3));
                    if par < 2.5882e-1   % sin(15) = 2.5882e-1;
                        error(['Segment "',SegmentName,'" -> Directions defining the segment are not adequate.\n', ...
                            '        A segment can be defined by:\n', ...
                            '           1) THREE vectors (not aligned) and ONE point\n', ...
                            '           2) TWO vectors and TWO points (defining three directions not aligned)']);
                    end
                end
            end
        end
        
        % -------------
        function parseSensors(HPxml,HumanModel)
            HumanSensors = HPxml.HumanDOM.getElementsByTagName('Sensor');
            for i = 1:HumanSensors.getLength
                SensorList = HumanSensors.item(i-1);
                SensorName = char(SensorList.getAttribute('Name'));
                SensorType = char(SensorList.getAttribute('Type'));
                InitPoint =  char(SensorList.getAttribute('InitPoint'));
                EndPoint =  char(SensorList.getAttribute('EndPoint'));
                Seg1 = char(SensorList.getAttribute('Seg1'));
                Seg2 = char(SensorList.getAttribute('Seg2'));
                RotSeq = char(SensorList.getAttribute('RotSeq'));
                
                % Check if Attributes are correct
                if isempty(SensorName), error('There is one "Sensor" element without attribute "Name"'); end
                if isempty(SensorType),
                    error(['Sensor "',SensorName,'" -> Attribute "Type" not defined. Options: TRS, SPH'])
                end
                
                if (strcmpi(SensorType,'TRS'))
                    if isempty(InitPoint)
                        error(['Sensor "',SensorName,' -> Attribute "InitPoint" not defined.'])                        
                    end
                    if isempty(EndPoint)
                        error(['Sensor "',SensorName,' -> Attribute "EndPoint" not defined.'])                        
                    end
                    
                elseif (strcmpi(SensorType,'SPH'))
                    if isempty(Seg1)
                        error(['Sensor "',SensorName,' -> Attribute "Seg1" not defined.'])                        
                    end
                    if isempty(Seg2)
                        error(['Sensor "',SensorName,' -> Attribute "Seg2" not defined.'])                        
                    end
                    if isempty(RotSeq)
                        error(['Sensor "',SensorName,' -> Attribute "RotSeq" not defined.'])
                    else
                        if strcmpi(RotSeq,'XYZ')
                            RotSeq = '123';
                        elseif strcmpi(RotSeq,'ZYX')
                            RotSeq = '321';
                        else
                            error(['Sensor "',SensorName,'" -> RotSeq="',RotSeq,'" is not valid. Options: XYZ or ZYX'])
                        end
                    end
                else
                    error(['Sensor "',SensorName,'" -> Type ="',SensorType,'" is not valid. Options: TRS, SPH'])
                end
                
                % If provided process Permutations
                Perm1 = char(SensorList.getAttribute('Perm1'));
                Perm2 = char(SensorList.getAttribute('Perm2'));
                Joint = char(SensorList.getAttribute('Joint'));
                if ~isempty(Perm1)
                    PosSep = strfind(Perm1,',');
                    Perm1x = str2num(Perm1(1:(PosSep-1)));
                    if isempty(Perm1x) % is char
                        Perm1x = Perm1(1:(PosSep-1));
                    end
                    Perm1y = str2num(Perm1((PosSep+1):end));
                    if isempty(Perm1y) % is char
                        Perm1y = Perm1((PosSep+1):end);
                    end                    
                else
                    Perm1x = [];
                    Perm1y = [];
                end
                if ~isempty(Perm2)
                    PosSep = strfind(Perm2,',');
                    Perm2x = str2num(Perm2(1:(PosSep-1)));
                    if isempty(Perm2x) % is char
                        Perm2x = Perm2(1:(PosSep-1));
                    end
                    Perm2y = str2num(Perm2((PosSep+1):end));
                    if isempty(Perm2y) % is char
                        Perm2y = Perm2((PosSep+1):end);
                    end
                else
                    Perm2x = [];
                    Perm2y = [];
                end
                
                % Add Sensor to the model
                HumanModel.addSensor(SensorName,SensorType,InitPoint,EndPoint,Seg1,Seg2,RotSeq,Perm1x,Perm1y,Perm2x,Perm2y);
                if ~isempty(Joint)
                    HumanModel.addSensorToJoint(Joint,SensorName);
                end
                
            end
        end
        % -------------
        function readxml(HPxml,PathXML,FileXML)
            HPxml.HumanDOM = xmlread([PathXML,FileXML]);
        end
        
        
    end
    
end

