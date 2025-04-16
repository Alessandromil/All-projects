classdef HUMAN_MODEL < handle 
    %HumanModel kinematic model description
    
    properties
        ModelName           % Name of the model                         char[1x1]
        CoordType           % Type of coord to define the model         char[1x1] Natural
        ModelType           % Type of model to write the results        char[1x1] Ramsis/PAM
        SubjectName         % Name of the Subject                                               char[1x1]
        Points              % Vector with all points which are not markers, AL nor Com.          POINT[]
        PointsRel           % Vector with all points rel of model       POINT[]
        Markers             % Vector with all markers of model.         POINT[]
        CoMs                % Vector with all center of mass in model.  POINT[]
        ALs                 % Vector with all anatomical landmarks.     POINT[]
        Vectors             % Vector with all vectors of model.         VECTOR[]
        Angles              % Vector with all angles of model.          ANGLE[]
        Segments            % Vector of Segments.                       SEGMENT[]
        Sensors             % Cell of Sensors.                          SENSOR{}
        SubAdditPar         % Subject additional Parameters not palpated  struct
                            % SubAdditPar.(ParameterName)
        q                   % Cell of all model variables.              char{}
        Joints = {};        % Cell of Joints.                           JOINT{}
        Phi                 % Vector of Phi sym/char ecuations          sym[]/char{NPhix1}
        Phiq                % structure with information for sparce matrix:
                            % Phiq.s: value of Phiq                     char{}
                            % Phiq.row:  row index vector               double[]
                            % Phiq.cols: col index vector               double[]
        Mass                % The mass of the subject                   double[]
        Gender              % The gender of the subject                 char{}
        SubjectParameters   % stucture vector with subject par.information( concretar con lo que es):
                            % SubjectParameters(i).Keyword              char
                            % SubjectParameters(i).Value                char
                            % SubjectParameters(i).Acronym              char
        MeasVAxis           % Is the Vertical Axis during the measurement double(3x1)
       
    end
    
    methods

        function HM = HUMAN_MODEL()
        end
        function AL = addAL(HM,ALName)
            % get position of the AL in vector, if AL isn't in vector,
            % add AL to vector. 
            ALIndex = getVecIndex(ALName,HM.ALs);
            if ALIndex == 0
                AL = POINT(ALName);
                HM.ALs = [HM.ALs;AL];
            else
                AL = HM.ALs(ALIndex);
            end
        end
        function addAngle(HM, Name1, Name2, Name3, Joint, Vec1Seg1, Vec1Seg2, Vec2Seg1, Vec2Seg2,Seg1VecRef,Seg2VecRef,RotSeq)
            NJoints = size(HM.Joints,1);
            for i= 1:NJoints
                if strcmpi(Joint,HM.Joints{i,1}.Name)
                    Seg1 = HM.Joints{i,1}.Seg1;
                    Seg2 = HM.Joints{i,1}.Seg2;
                    if(strcmpi(HM.Joints{i,1}.Type,'SPH'))
                        % Return the LOCAL_VECTORs of Segments
                        LV1Seg1 = getVecInSeg(HM.Segments,Seg1,Vec1Seg1);
                        LV1Seg2 = getVecInSeg(HM.Segments,Seg2,Vec1Seg2);
                        LV2Seg1 = getVecInSeg(HM.Segments,Seg1,Vec2Seg1);
                        LV2Seg2 = getVecInSeg(HM.Segments,Seg2,Vec2Seg2);
                        
                        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        % FUNCTION getVecInSeg has sido modificada. Ahora devuelve la variables vacias
                        % si no la encuentra en lugar de dar error. Es lo mismo que se ha hecho para getPointInSeg
                        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        
                        % Add 
                        HM.Angles = [HM.Angles; ANGLE(Name1,Name2,Name3,HM.Joints{i,1},LV1Seg1,LV1Seg2,LV2Seg1,LV2Seg2,RotSeq)];
                    elseif (strcmpi(HM.Joints{i,1}.Type,'UNI'))
                        % Return the LOCAL_VECTORs of Segments
                        LV1Seg1 = getVecInSeg(HM.Segments,Seg1,Seg1VecRef);
                        LV1Seg2 = getVecInSeg(HM.Segments,Seg2,Seg2VecRef);
                        % Add 
                        HM.Angles = [HM.Angles; ANGLE(Name1,Name2,Name3,HM.Joints{i,1},LV1Seg1,LV1Seg2,0,0,0)];
                    elseif(strcmpi(HM.Joints{i,1}.Type,'REV'))
                         % Return the LOCAL_VECTORs of Segments
                        LV1Seg1 = getVecInSeg(HM.Segments,Seg1,Seg1VecRef);
                        LV1Seg2 = getVecInSeg(HM.Segments,Seg2,Seg2VecRef);
                        % Add 
                        HM.Angles = [HM.Angles; ANGLE(Name1,Name2,Name3,HM.Joints{i,1},LV1Seg1,LV1Seg2,0,0,0)];
                    end
                    return;
                end
            end
        end
        function CoM = addCoM(HM,CoMName)
            % Calc position the point in vector, if point isn't in vector,
            % add point to vector. 
            CoMIndex = getVecIndex(CoMName,HM.CoMs);
            if CoMIndex == 0
                CoM = POINT(CoMName);
                HM.CoMs = [HM.CoMs;CoM];
            else
                CoM = HM.CoMs(CoMIndex);
            end
        end
        %---------
        function addJoint(HM,Name,Type,Seg1,Seg2,Point1,Point2,Vector1,Vector2,V1_Ang_V2)
            
            % According to the type of joint, create one and add to the list of joints            
            if strcmpi(Type,'UNI')
                
                % Check if all Attributes for a UNI joint are defined
                if isempty(Seg1), error(['Joint "',Name,'" -> Attribute "Seg1" not defined.']); end
                if isempty(Seg2), error(['Joint "',Name,'" -> Attribute "Seg2" not defined.']); end
                if isempty(Point1), error(['Joint "',Name,'" -> Attribute "Point1" not defined.']); end
                if isempty(Point2), error(['Joint "',Name,'" -> Attribute "Point2" not defined.']); end
                if isempty(Vector1), error(['Joint "',Name,'" -> Attribute "Seg1Axis" not defined.']); end
                if isempty(Vector2), error(['Joint "',Name,'" -> Attribute "Seg2Axis" not defined.']); end
                if isempty(V1_Ang_V2), error(['Joint "',Name,'" -> Attribute "AxesAng" not defined.']); end

                
                % Return the LOCAL_POINTs of segments
                LPoint1 = getPointInSeg(HM.Segments,Seg1,Point1); 
                LPoint2 = getPointInSeg(HM.Segments,Seg2,Point2); 
                % Return the LOCAL_VECTORs of Segments
                LVector1 = getVecInSeg(HM.Segments,Seg1,Vector1);
                LVector2 = getVecInSeg(HM.Segments,Seg2,Vector2);
                
                % Check if segment and point exists in the model
                if isempty(LPoint1)
                    error(['Joint "',Name,'" -> Point1="',Point1,'" does not belong to segment "',Seg1,...
                           '"\n                    or segment "',Seg1,'" is not defined.']);
                end
                if isempty(LPoint2)
                    error(['Joint "',Name,'" -> Point2="',Point2,'" does not belong to segment "',Seg2,...
                           '"\n                    or segment "',Seg2,'" is not defined.']);
                end

                % Check if segment and vector exists in the model
                if isempty(LVector1)
                    error(['Joint "',Name,'" -> Vector1="',Vector1,'" does not belong to segment "',Seg1,...
                           '"\n                    or segment "',Seg1,'" is not defined.']);
                end
                if isempty(LVector2)
                    error(['Joint "',Name,'" -> Vector2="',Vector2,'" does not belong to segment "',Seg2,...
                           '"\n                    or segment "',Seg2,'" is not defined.']);
                end                
                
                % Add universal joint to the model
                NJoints = size(HM.Joints,1);
                HM.Joints{NJoints+1,1} = UNIJOINT(Name,Type,Seg1,Seg2,LPoint1,LPoint2,LVector1,LVector2,V1_Ang_V2);
                
            elseif strcmpi(Type,'SPH')
                
                % Check if all Attributes for a SPH joint are defined
                if isempty(Seg1), error(['Joint "',Name,'" -> Attribute "Seg1" not defined.']); end
                if isempty(Seg2), error(['Joint "',Name,'" -> Attribute "Seg2" not defined.']); end
                if isempty(Point1), error(['Joint "',Name,'" -> Attribute "Point1" not defined.']); end
                if isempty(Point2), error(['Joint "',Name,'" -> Attribute "Point2" not defined.']); end
                
                % Return the LOCAL_POINTs of segments
                LPoint1 = getPointInSeg(HM.Segments,Seg1,Point1);
                LPoint2 = getPointInSeg(HM.Segments,Seg2,Point2);
                
                % Check if segment and point exists in the model
                if isempty(LPoint1)
                    error(['Joint "',Name,'" -> Point1="',Point1,'" does not belong to segment "',Seg1,...
                        '"\n                    or segment "',Seg1,'" is not defined.']);
                end
                if isempty(LPoint2)
                    error(['Joint "',Name,'" -> Point2="',Point2,'" does not belong to segment "',Seg2,...
                        '"\n                    or segment "',Seg2,'" is not defined.']);
                end
                
                % Add spherical joint to the model
                NJoints = size(HM.Joints,1);
                HM.Joints{NJoints+1,1} = SPHJOINT(Name,Type,Seg1,Seg2,LPoint1,LPoint2);
                
            elseif strcmpi(Type,'REV')
                
                % Check if all Attributes for a REV joint are defined
                if isempty(Seg1), error(['Joint "',Name,'" -> Attribute "Seg1" not defined.']); end
                if isempty(Seg2), error(['Joint "',Name,'" -> Attribute "Seg2" not defined.']); end
                if isempty(Point1), error(['Joint "',Name,'" -> Attribute "Point1" not defined.']); end
                if isempty(Point2), error(['Joint "',Name,'" -> Attribute "Point2" not defined.']); end
                if isempty(Vector1), error(['Joint "',Name,'" -> Attribute "Seg1Axis" not defined.']); end
                if isempty(Vector2), error(['Joint "',Name,'" -> Attribute "Seg2Axis" not defined.']); end
                
                % Return the LOCAL_POINTs of segments
                LPoint1 = getPointInSeg(HM.Segments,Seg1,Point1);
                LPoint2 = getPointInSeg(HM.Segments,Seg2,Point2);
                % Return the LOCAL_VECTORs of Segments
                LVector1 = getVecInSeg(HM.Segments,Seg1,Vector1);
                LVector2 = getVecInSeg(HM.Segments,Seg2,Vector2);
                
                % Check if segment and point exists in the model
                if isempty(LPoint1)
                    error(['Joint "',Name,'" -> Point1="',Point1,'" does not belong to segment "',Seg1,...
                           '"\n                    or segment "',Seg1,'" is not defined.']);
                end
                if isempty(LPoint2)
                    error(['Joint "',Name,'" -> Point2="',Point2,'" does not belong to segment "',Seg2,...
                           '"\n                    or segment "',Seg2,'" is not defined.']);
                end

                % Check if segment and vector exists in the model
                if isempty(LVector1)
                    error(['Joint "',Name,'" -> Vector1="',Vector1,'" does not belong to segment "',Seg1,...
                           '"\n                    or segment "',Seg1,'" is not defined.']);
                end
                if isempty(LVector2)
                    error(['Joint "',Name,'" -> Vector2="',Vector2,'" does not belong to segment "',Seg2,...
                           '"\n                    or segment "',Seg2,'" is not defined.']);
                end
                
                % Add revolution joint to the model
                NJoints = size(HM.Joints,1);
                HM.Joints{NJoints+1,1} = REVJOINT(Name,Type,Seg1,Seg2,LPoint1,LPoint2,LVector1,LVector2);
                
            elseif strcmpi(Type,'FLOAT')
                
                % Check if all Attributes for a FLOAT joint are defined
                if isempty(Seg1), error(['Joint "',Name,'" -> Attribute "Seg1" not defined.']); end
                if isempty(Seg2), error(['Joint "',Name,'" -> Attribute "Seg2" not defined.']); end
                
                % Check is segments belong to the model
                Seg1Index = getSegInModel(HM.Segments,Seg1);
                Seg2Index = getSegInModel(HM.Segments,Seg2);

                % Check if segment and point exists in the model
                if isempty(Seg1Index)
                    error(['Joint "',Name,'" -> Seg1="',Seg1,'" not defined for this model']);
                end
                if isempty(Seg2Index)
                    error(['Joint "',Name,'" -> Seg2="',Seg2,'" not defined for this model']);
                end                
                
                % Add float joint to the model
                NJoints = size(HM.Joints,1);
                HM.Joints{NJoints+1,1} = FLOATJOINT(Name,Type,Seg1,Seg2);
                
            else           
                error(['The Joint "',Name,'" is not a correct type of Joint'])
            end
        end
        
        %---------
        function addModelJoints(HM)
            % get segments index
            PelvisIndex = getVecIndex('Pelvis',HM.Segments);
            RThighIndex = getVecIndex('RightThigh',HM.Segments);
            LThighIndex = getVecIndex('LeftThigh',HM.Segments);
            RShankIndex = getVecIndex('RightShank',HM.Segments);
            LShankIndex = getVecIndex('LeftShank',HM.Segments);
            RFootIndex = getVecIndex('RightFoot',HM.Segments);            
            LFootIndex = getVecIndex('LeftFoot',HM.Segments);
            ThoraxIndex = getVecIndex('Thorax',HM.Segments);
            RClavicleIndex = getVecIndex('RightClavicle',HM.Segments);
            RScapulaIndex = getVecIndex('RightScapula',HM.Segments);
            RHumerusIndex = getVecIndex('RightHumerus',HM.Segments);
            RForearmIndex = getVecIndex('RightForearm',HM.Segments);
            RUlnaIndex = getVecIndex('RightUlna',HM.Segments);
            RRadiusIndex = getVecIndex('RightRadius',HM.Segments);
            RHandIndex = getVecIndex('RightHand',HM.Segments);
            LClavicleIndex = getVecIndex('LeftClavicle',HM.Segments);
            LScapulaIndex = getVecIndex('LeftScapula',HM.Segments);
            LHumerusIndex = getVecIndex('LeftHumerus',HM.Segments);
            LForearmIndex = getVecIndex('LeftForearm',HM.Segments);
            LUlnaIndex = getVecIndex('LeftUlna',HM.Segments);
            LRadiusIndex = getVecIndex('LeftRadius',HM.Segments);
            LHandIndex = getVecIndex('LeftHand',HM.Segments);
            % add points to the model
            if PelvisIndex ~= 0
                % Righ Hip Joint
                AddedPoint = HM.addPoint('RHJC');
                HM.Segments(PelvisIndex).addPoint(AddedPoint);
                % Left Hip Joint
                AddedPoint = HM.addPoint('LHJC');
                HM.Segments(PelvisIndex).addPoint(AddedPoint);
                % Medial Lumbar Joint (for the moment only point is added)
                AddedPoint = HM.addPoint('MLJC');
                HM.Segments(PelvisIndex).addPoint(AddedPoint);
            end
            if RThighIndex ~= 0 
                % Righ Hip Joint
                AddedPoint = HM.addPoint('RHJC');
                HM.Segments(RThighIndex).addPoint(AddedPoint);
                % Righ Knee Joint 
                AddedPoint = HM.addPoint('RKJC');
                HM.Segments(RThighIndex).addPoint(AddedPoint);
            end
            if RShankIndex ~= 0
                % Righ Knee Joint
                AddedPoint = HM.addPoint('RKJC');
                HM.Segments(RShankIndex).addPoint(AddedPoint);
                % Righ Ankle Joint
                AddedPoint = HM.addPoint('RAJC');
                HM.Segments(RShankIndex).addPoint(AddedPoint);
            end
            if RFootIndex ~= 0
                % Righ Ankle Joint
                AddedPoint = HM.addPoint('RAJC');
                HM.Segments(RFootIndex).addPoint(AddedPoint);
            end
            if LThighIndex ~= 0 
                % Left Hip Joint
                AddedPoint = HM.addPoint('LHJC');
                HM.Segments(LThighIndex).addPoint(AddedPoint);
                % Left Knee Joint 
                AddedPoint = HM.addPoint('LKJC');
                HM.Segments(LThighIndex).addPoint(AddedPoint);
            end
            if LShankIndex ~= 0
                % Left Knee Joint
                AddedPoint = HM.addPoint('LKJC');
                HM.Segments(LShankIndex).addPoint(AddedPoint);
                % Left Ankle Joint
                AddedPoint = HM.addPoint('LAJC');
                HM.Segments(LShankIndex).addPoint(AddedPoint);
            end
            if LFootIndex ~= 0
                % Left Ankle Joint
                AddedPoint = HM.addPoint('LAJC');
                HM.Segments(LFootIndex).addPoint(AddedPoint);
            end
            % Add joints to the model
            if PelvisIndex ~= 0 && LThighIndex ~= 0
                HM.addJoint('LHJC','SPH','Pelvis','LeftThigh','LHJC','LHJC',[],[],[]);
            end
            if LThighIndex ~= 0 && LShankIndex ~= 0
                HM.addJoint('LKJC','SPH','LeftThigh','LeftShank','LKJC','LKJC',[],[],[]);
            end
            if LShankIndex ~= 0 && LFootIndex ~= 0
                HM.addJoint('LAJC','SPH','LeftShank','LeftFoot','LAJC','LAJC',[],[],[]);
            end
            if PelvisIndex ~= 0 && RThighIndex ~= 0
                HM.addJoint('RHJC','SPH','Pelvis','RightThigh','RHJC','RHJC',[],[],[]);
            end
            if RThighIndex ~= 0 && RShankIndex ~= 0
                HM.addJoint('RKJC','SPH','RightThigh','RightShank','RKJC','RKJC',[],[],[]);
            end
            if RShankIndex ~= 0 && RFootIndex ~= 0
                HM.addJoint('RAJC','SPH','RightShank','RightFoot','RAJC','RAJC',[],[],[]);
            end
                
%                 if LThighIndex ~= 0
%                      HM.Segments(LThighIndex).addPoint(AddedPoint);
%                      HM.addJoint('LHJC','SPH','Pelvis','LeftThigh','LHJC','LHJC',[],[],[],[]);
%                      % Left Knee Joint 
%                     AddedPoint = HM.addPoint('LKJC');
%                     HM.Segments(LThighIndex).addPoint(AddedPoint);
%                     if LShankIndex ~= 0
%                         HM.Segments(LShankIndex).addPoint(AddedPoint);
%                         HM.addJoint('LKJC','SPH','LeftThigh','LeftShank','LKJC','LKJC',[],[],[],[]);
%                         % Left Ankle Joint
%                         AddedPoint = HM.addPoint('LAJC');
%                         HM.Segments(LShankIndex).addPoint(AddedPoint);
%                         if LFootIndex ~= 0
%                             HM.Segments(LFootIndex).addPoint(AddedPoint);
%                             HM.addJoint('LAJC','SPH','LeftShank','LeftFoot','LAJC','LAJC',[],[],[],[]);
%                         end
%                     end
%                 end
                
            if ThoraxIndex ~= 0
                % Right Sternoclavicular joint
                AddedPoint = HM.addPoint('RSCJC');
                HM.Segments(ThoraxIndex).addPoint(AddedPoint);
                if RClavicleIndex ~= 0
                    HM.Segments(RClavicleIndex).addPoint(AddedPoint);
                    HM.addJoint('RSCJC','SPH','Thorax','RightClavicle','RSCJC','RSCJC',[],[],[]);
                    % Right Acromioclavicular joint
                    AddedPoint = HM.addPoint('RACJC');
                    HM.Segments(RClavicleIndex).addPoint(AddedPoint);
                    if RScapulaIndex ~= 0
                        HM.Segments(RScapulaIndex).addPoint(AddedPoint);
                        HM.addJoint('RACJC','SPH','RightClavicle','RightScapula','RACJC','RACJC',[],[],[]);
                        % Right glenohumeral joint
                        AddedPoint = HM.addPoint('RGHJC');
                        HM.Segments(RScapulaIndex).addPoint(AddedPoint);
                        if RHumerusIndex ~= 0
                            HM.Segments(RHumerusIndex).addPoint(AddedPoint);
                            HM.addJoint('RGHJC','SPH','RightScapula','RightHumerus','RGHJC','RGHJC',[],[],[]);
                            % Right HumeroUlnar joint centre
                            AddedPoint = HM.addPoint('RHUJC');
                            HM.Segments(RHumerusIndex).addPoint(AddedPoint);
                            if RUlnaIndex ~= 0
                                HM.Segments(RUlnaIndex).addPoint(AddedPoint);
                                HM.addJoint('RHUJC','SPH','RightHumerus','RightUlna','RHUJC','RHUJC',[],[],[]);
                                % Right Radio-ulnar joint centre
                                AddedPoint = HM.addPoint('RRUJC');
                                HM.Segments(RUlnaIndex).addPoint(AddedPoint);
                                if RRadiusIndex ~= 0
                                    HM.Segments(RRadiusIndex).addPoint(AddedPoint);
                                    HM.addJoint('RRUJC','SPH','RightUlna','RightRadius','RRUJC','RRUJC',[],[],[]);
                                    AddedPoint = HM.addPoint('RWJC');
                                    HM.Segments(RRadiusIndex).addPoint(AddedPoint);
                                end
                            end
                            % Right Elbow joint centre
                            AddedPoint = HM.addPoint('REJC');
                            HM.Segments(RHumerusIndex).addPoint(AddedPoint);
                            if RForearmIndex ~= 0
                                HM.Segments(RForearmIndex).addPoint(AddedPoint);
                                HM.addJoint('REJC','SPH','RightHumerus','RightForearm','REJC','REJC',[],[],[]);
                                AddedPoint = HM.addPoint('RWJC');
                                HM.Segments(RForearmIndex).addPoint(AddedPoint);
                            end
                            if RHandIndex ~= 0
                                % Right Wrist joint centre
                                AddedPoint = HM.addPoint('RWJC');
                                HM.Segments(RHandIndex).addPoint(AddedPoint);
                                if RRadiusIndex ~= 0
                                    HM.addJoint('RWJC','SPH','RightRadius','RightHand','RWJC','RWJC',[],[],[]);
                                end
                                if RForearmIndex ~= 0
                                    HM.addJoint('RWJC','SPH','RightForearm','RightHand','RWJC','RWJC',[],[],[]);
                                end
                            end
                                
                        end
                    end
                end
                % Left Sternoclavicular joint
                AddedPoint = HM.addPoint('LSCJC');
                HM.Segments(ThoraxIndex).addPoint(AddedPoint);
                if LClavicleIndex ~= 0
                    HM.Segments(LClavicleIndex).addPoint(AddedPoint);
                    HM.addJoint('LSCJC','SPH','Thorax','LeftClavicle','LSCJC','LSCJC',[],[],[]);
                    % Left Acromioclavicular joint
                    AddedPoint = HM.addPoint('LACJC');
                    HM.Segments(LClavicleIndex).addPoint(AddedPoint);
                    if LScapulaIndex ~= 0
                        HM.Segments(LScapulaIndex).addPoint(AddedPoint);
                        HM.addJoint('LACJC','SPH','LeftClavicle','LeftScapula','LACJC','LACJC',[],[],[]);
                        % Left glenohumeral joint
                        AddedPoint = HM.addPoint('LGHJC');
                        HM.Segments(LScapulaIndex).addPoint(AddedPoint);
                        if LHumerusIndex ~= 0
                            HM.Segments(LHumerusIndex).addPoint(AddedPoint);
                            HM.addJoint('LGHJC','SPH','LeftScapula','LeftHumerus','LGHJC','LGHJC',[],[],[]);
                            % Left HumeroUlnar joint centre
                            AddedPoint = HM.addPoint('LHUJC');
                            HM.Segments(LHumerusIndex).addPoint(AddedPoint);
                            if LUlnaIndex ~= 0
                                HM.Segments(LUlnaIndex).addPoint(AddedPoint);
                                HM.addJoint('LHUJC','SPH','LeftHumerus','LeftUlna','LHUJC','LHUJC',[],[],[]);
                                % Left Radio-ulnar joint centre
                                AddedPoint = HM.addPoint('LRUJC');
                                HM.Segments(LUlnaIndex).addPoint(AddedPoint);
                                if LRadiusIndex ~= 0
                                    HM.Segments(LRadiusIndex).addPoint(AddedPoint);
                                    HM.addJoint('LRUJC','SPH','LeftUlna','LeftRadius','LRUJC','LRUJC',[],[],[]);
                                    AddedPoint = HM.addPoint('LWJC');
                                    HM.Segments(LRadiusIndex).addPoint(AddedPoint);
                                end
                            end
                            % Left Elbow joint centre
                            AddedPoint = HM.addPoint('LEJC');
                            HM.Segments(LHumerusIndex).addPoint(AddedPoint);
                            if LForearmIndex ~= 0
                                HM.Segments(LForearmIndex).addPoint(AddedPoint);
                                HM.addJoint('LEJC','SPH','LeftHumerus','LeftForearm','LEJC','LEJC',[],[],[]);
                                AddedPoint = HM.addPoint('LWJC');
                                HM.Segments(LForearmIndex).addPoint(AddedPoint);
                            end
                            if LHandIndex ~= 0
                                % Left Wrist joint centre
                                AddedPoint = HM.addPoint('LWJC');
                                HM.Segments(LHandIndex).addPoint(AddedPoint);
                                if LRadiusIndex ~= 0
                                    HM.addJoint('LWJC','SPH','LeftRadius','LeftHand','LWJC','LWJC',[],[],[]);
                                end
                                if LForearmIndex ~= 0
                                    HM.addJoint('LWJC','SPH','LeftForearm','LeftHand','LWJC','LWJC',[],[],[]);
                                end
                            end
                                
                        end
                    end
                end
                
            end
        end
        function Marker = addMarker(HM,PointName)
            % Calc position the point in vector, if point isn't in vector,
            % add point to vector. 
            MarkerIndex = getVecIndex(PointName,HM.Markers);
            if MarkerIndex == 0
                Marker = POINT(PointName);
                HM.Markers = [HM.Markers;Marker];
            else
                Marker = HM.Markers(MarkerIndex);
            end
        end
        function addMarkerInq(HM)
            NMarkers = size(HM.Markers,1);
            for i=1:NMarkers
                Nq = size(HM.q,1);
                MarkerName = HM.Markers(i).CoordName;
                HM.Markers(i).PosInq = Nq+1;
                HM.q{Nq+1,1} = MarkerName{1};
                HM.q{Nq+2,1} = MarkerName{2};
                HM.q{Nq+3,1} = MarkerName{3};
            end
        end
        function Point = addPoint(HM,PointName)
            % Calc position the point in vector, if point isn't in vector,
            % add point to vector. 
            PointIndex = getVecIndex(PointName,HM.Points);
            if PointIndex == 0
                Point = POINT(PointName);
                HM.Points = [HM.Points;Point];
            else
                Point = HM.Points(PointIndex);
            end
        end
        function PointRel = addRelPoint(HM,PointRelName)
            % Calc position the point in vector, if point isn't in vector,
            % add point to vector. 
            PointIndex = getVecIndex(PointRelName,HM.PointsRel);
            if PointIndex == 0
                PointRel = POINT(PointRelName);
                HM.PointsRel = [HM.PointsRel;PointRel];
            else
                PointRel = HM.PointsRel(PointIndex);
            end 
        end
        function addSegment(HM,SegmentName)
            HM.Segments =[HM.Segments; SEGMENT(SegmentName)];
        end
        function addSensor(HM,Name,Type,InitPoint,EndPoint,Seg1,Seg2,RotSeq,Perm1x,Perm1y,Perm2x,Perm2y)
            NSensors = size(HM.Sensors,1);
            if (strcmpi(Type,'TRS'))
                InitPIndex = getVecIndex(InitPoint,HM.Points);
                EndPIndex  = getVecIndex(EndPoint,HM.Points);
                if InitPIndex==0, error(['Sensor "',Name,'" -> Point "',InitPoint,'" does not belong to any segment.']); end
                if EndPIndex==0, error(['Sensor "',Name,'" -> Point "',EndPoint,'" does not belong to any segment.']); end                
                HM.Sensors{NSensors+1,1} = TRS_SENSOR(Name,Type,HM.Points(InitPIndex),HM.Points(EndPIndex));
            elseif (strcmpi(Type,'SPH'))
                Seg1Index = getVecIndex(Seg1,HM.Segments);
                Seg2Index = getVecIndex(Seg2,HM.Segments);
                if Seg1Index==0, error (['Sensor "',Name,'" -> Segment "',Seg1,'" is not defined for this model.']); end
                if Seg2Index==0, error (['Sensor "',Name,'" -> Segment "',Seg2,'" is not defined for this model.']); end
                HM.Sensors{NSensors+1,1} = SPH_SENSOR(Name,Type,HM.Segments(Seg1Index),HM.Segments(Seg2Index),RotSeq,...
                                                        Perm1x,Perm1y,Perm2x,Perm2y);
            else
                error([Type,' is an incorrec type of sensor']);
            end
            
        end
        function addSensorToJoint(HM,JointName,SensorName)
            SensorIndex  = getVecIndex(SensorName,HM.Sensors);
            JointIndex   = getVecIndex(JointName,HM.Joints);
            if strcmpi(HM.Joints{JointIndex,1}.Type,'Float') && strcmpi(HM.Sensors{SensorIndex}.Type,'TRS')
                HM.Joints{JointIndex,1}.TrsSensor = HM.Sensors{SensorIndex};
            else
                HM.Joints{JointIndex,1}.Sensor = HM.Sensors{SensorIndex};
            end
            
        end
        function Vector  = addVector(HM,VectorName)
            % Calc position the vector in vector, if vector isn't in vector,
            % add vector to vector.            
            VectorIndex = getVecIndex(VectorName,HM.Vectors);
            if VectorIndex == 0
                Vector = VECTOR(VectorName);
                HM.Vectors = [HM.Vectors;Vector];
            else
                Vector = HM.Vectors(VectorIndex);
            end
        end 
        function addVector2Segment (HM,SegName,VecName,LocCoord)
            VectorIndex  = getVecIndex(VecName,HM.Vectors);
            SegmentIndex = getVecIndex(SegName,HM.Segments);
            if VectorIndex == 0
                error ('%s does not belong to model.',VecName);
            elseif SegmentIndex == 0
                error ('%s does not belong to model.',SegName);
            else
               Vector = HM.Vectors(VectorIndex);
               HM.Segments(SegmentIndex).addVector2Segment(Vector,LocCoord);
               
            end
                
        end
        function calcALLCS_T_RamsisLCS(HM)
            NSegments = size(HM.Segments,1);
            for i=1:NSegments
                HM.Segments(i).calcALLCS_T_RamsisLCS(HM.MeasVAxis);
            end
        end
        function calcInertiaPars(HM) 
           NSegments = size(HM.Segments,1); 
           for i=1:NSegments
               CoMName = HM.Segments(i).Name;
               CoM = addCoM(HM,CoMName);
               HM.Segments(i).calcInertiaPars(HM.Mass,CoM);
            end
        end
        function calcLCS(HM)
            NSegments = size(HM.Segments,1);
            for i=1:NSegments
                HM.Segments(i).calcLCS();
            end
        end
        function calcTCS(HM)
            NSegments = size(HM.Segments,1);
            % To calculate TCS
            for i=1:NSegments
                HM.Segments(i).calcTCS();
            end
        end
        function calcTransferedPoint(HM)
            % To calculate transferece point for LCS
            PelvisIndex = getVecIndex('Pelvis',HM.Segments);
            LShankIndex = getVecIndex('LeftShank',HM.Segments);
            LFootIndex = getVecIndex('LeftFoot',HM.Segments);
            LThighIndex = getVecIndex('LeftThigh',HM.Segments);
            RShankIndex = getVecIndex('RightShank',HM.Segments);
            RFootIndex = getVecIndex('RightFoot',HM.Segments);
            RThighIndex = getVecIndex('RightThigh',HM.Segments);
            ThoraxIndex = getVecIndex('Thorax',HM.Segments);
            RClavicleIndex = getVecIndex('RightClavicle',HM.Segments);
            LClavicleIndex = getVecIndex('LeftClavicle',HM.Segments);
            RHumerusIndex = getVecIndex('RightHumerus',HM.Segments);
            LHumerusIndex = getVecIndex('LeftHumerus',HM.Segments);
            RForearmIndex = getVecIndex('RightForearm',HM.Segments);
            LForearmIndex = getVecIndex('LeftForearm',HM.Segments);
            RUlnaIndex = getVecIndex('RightUlna',HM.Segments);
            LUlnaIndex = getVecIndex('LeftUlna',HM.Segments);
            RRadiusIndex = getVecIndex('RightRadius',HM.Segments);
            LRadiusIndex = getVecIndex('LeftRadius',HM.Segments);
            RHandIndex = getVecIndex('RightHand',HM.Segments);
            LHandIndex = getVecIndex('LeftHand',HM.Segments);
            TSegsIndices = [];
            Post_Pos_TOSegMarkers = [];
            TCSTO_Pos_TOSegMarkers= [];
            TCSTO_Pos_TOSegALs    = [];
            TOLocalALs = [];
            PosturNames = [];
            
            % First posture in ALM file is used for transfering data to next segment
            PostureName = HM.Segments(1).LocalMarkers(1).Point.Postures(1).Name;
            
            if RThighIndex ~=0 && PelvisIndex ~=0
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(RThighIndex,PelvisIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if LThighIndex ~=0 && PelvisIndex ~=0
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(LThighIndex,PelvisIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if RThighIndex ~=0 && RShankIndex ~=0
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(RShankIndex,RThighIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if RShankIndex ~=0 && RFootIndex ~=0
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(RFootIndex,RShankIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if LThighIndex ~=0 && LShankIndex ~=0
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(LShankIndex,LThighIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if LShankIndex ~=0 && LFootIndex ~=0
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(LFootIndex,LShankIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if ThoraxIndex ~= 0 && RClavicleIndex ~= 0
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(RClavicleIndex,ThoraxIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if ThoraxIndex ~= 0 && LClavicleIndex ~= 0
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(LClavicleIndex,ThoraxIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if RHumerusIndex ~= 0 && RForearmIndex ~= 0
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(RHumerusIndex,RForearmIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(RForearmIndex,RHumerusIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if LHumerusIndex ~= 0 && LForearmIndex ~= 0
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(LHumerusIndex,LForearmIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(LForearmIndex,LHumerusIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if RUlnaIndex ~=0 && RHumerusIndex ~=0
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(RUlnaIndex,RHumerusIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if LUlnaIndex ~=0 && LHumerusIndex ~=0
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(LUlnaIndex,LHumerusIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if RRadiusIndex ~=0 && RHumerusIndex ~=0 && RUlnaIndex ~=0
                TransfOrSegIndex = [RHumerusIndex;RUlnaIndex];
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(RRadiusIndex,TransfOrSegIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if LRadiusIndex ~=0 && LHumerusIndex ~=0 && LUlnaIndex ~=0
                TransfOrSegIndex = [LHumerusIndex;LUlnaIndex];
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(LRadiusIndex,TransfOrSegIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if RHandIndex ~=0 && RRadiusIndex ~=0
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(RHandIndex,RRadiusIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            if LHandIndex ~=0 && LRadiusIndex ~=0
                [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = HM.getDataToTransfPoints(LHandIndex,LRadiusIndex,PostureName,...
                                                                                        Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames);
            end
            for i=TSegsIndices
                HM.Segments(i).calcTransferedALs(Post_Pos_TOSegMarkers(i,:),TCSTO_Pos_TOSegMarkers(i,:),TCSTO_Pos_TOSegALs(i,:),TOLocalALs(i,:),PosturNames(i));
            end
        end
        function [Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames] = getDataToTransfPoints(HM,TransfSegIndex,TransfOrSegIndex,PostureName,...
                                                                                                            Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,TSegsIndices,PosturNames)
        % TSegsIndices     = All the segments to transfer points. Vector [1x NTSegsIndices]
        % TransfSegIndex   = The index of Segment to transfer the points. double
        % TransfOrSegIndex = The indices of Segment from transfered points. Vector [NTransfOrSeg x 1]
        % Post_Pos_TOSegMarkers  = Position of the markers of transfer origin segment in the Posture (Global) coordinate system. Vector[3xNMarkersTOSeg] 
        % TCSTO_Pos_TOSegMarkers = Position of the markers of transfer origin segment in the transfer origin segment TCS. Vector[3x,NMarkersTOSeg]
        % TCSTO_Pos_TOSegALs     = Position of the AL of transfer origin segment in the transfer origin segment TCS. Vector[3xNALsTOSeg]
        % TOLocalALs             = Vector with all the ALs in the transfer origin segment. Vector[NALsTOSegx1]
            NTransfOrSeg = size(TransfOrSegIndex,1);
            PostureIndex = getVecIndex(PostureName,HM.Segments(TransfOrSegIndex(1)).LocalMarkers(1).Point.Postures);
            PosturNames{TransfSegIndex} = PostureName;
            for i=1:NTransfOrSeg
                NMarkersTOSeg = size(HM.Segments(TransfOrSegIndex(i)).LocalMarkers,1);
                for j=1:NMarkersTOSeg
                    Post_Pos_TOSegMarkers_j(:,j) = HM.Segments(TransfOrSegIndex(i)).LocalMarkers(j).Point.Postures(PostureIndex).Glob;
                end
                Post_Pos_TOSegMarkers{TransfSegIndex,i}  = Post_Pos_TOSegMarkers_j;
                TCSTO_Pos_TOSegMarkers{TransfSegIndex,i} = HM.Segments(TransfOrSegIndex(i)).TCS_Pos_Markers;
                TCSTO_Pos_TOSegALs{TransfSegIndex,i}     = HM.Segments(TransfOrSegIndex(i)).TCS_Pos_Landmark;
                TOLocalALs{TransfSegIndex,i}             = HM.Segments(TransfOrSegIndex(i)).LocalALs;
            end
            TSegsIndices = [TSegsIndices,TransfSegIndex];
        end
        function checkSubjectParsData(HM,File)
            % CHECKSUBJECTPARSDATA cheks that all parameters of the model are in the Subject Parameter file
            NSegments = size(HM.Segments,1);
            NSensors = size(HM.Sensors,1);
            for i=1:NSegments
                HM.Segments(i).checkSubjectParsData(File)
            end
            for i=1:NSensors
                if strcmp(HM.Sensors{i}.Type,'SPH')
                    HM.Sensors{i}.checkSubjectParsData(File)
                end
                
            end
        end
        function drawLCS(HM)
            NSegments = size(HM.Segments,1);
            for i=1:NSegments
                %if ~strcmpi(HM.Segments(i).Name,'Ground')
                if HM.Segments(i).Fixed ~= 1 % If is not Ground                    
                    HM.Segments(i).drawLCS;
                end
            end
        end
        function fillq(HM)
            % fill q with problem variables
            % initialize
            NCoords = 1;
            NPoints  = size(HM.Points,1);
            NPointsRel = size(HM.PointsRel,1);
            NVectors = size(HM.Vectors,1);
            NAngles  = size(HM.Angles,1);
            % fill q with points
            for i = 1:NPoints
                % Test if point is fixed 
                if HM.Points(i).Fixed == 0
                    HM.Points(i).PosInq = NCoords;
                    HM.q{NCoords,1}   = HM.Points(i).CoordName{1};
                    HM.q{NCoords+1,1} = HM.Points(i).CoordName{2};
                    HM.q{NCoords+2,1} = HM.Points(i).CoordName{3};
                    NCoords = NCoords + 3;
                end
            end
            for i = 1:NPointsRel
                HM.PointsRel(i).PosInq = NCoords;
                HM.q{NCoords,1}   = HM.PointsRel(i).CoordName{1};
                HM.q{NCoords+1,1} = HM.PointsRel(i).CoordName{2};
                HM.q{NCoords+2,1} = HM.PointsRel(i).CoordName{3};
                NCoords = NCoords + 3;
            end
            % fill q with vectors
            for i = 1:NVectors
                % Test if vector is fixed 
                if HM.Vectors(i).Fixed == 0
                    HM.Vectors(i).PosInq = NCoords;
                    HM.q{NCoords,1}   = HM.Vectors(i).CoordName{1};
                    HM.q{NCoords+1,1} = HM.Vectors(i).CoordName{2};
                    HM.q{NCoords+2,1} = HM.Vectors(i).CoordName{3};
                    NCoords = NCoords + 3;
                end
            end
            % fill q with guided angles
            for i = 1:NAngles
                if(strcmpi(HM.Angles(i).Joint.Type,'SPH'))
                    HM.Angles(i).PosInq = [NCoords,NCoords+1,NCoords+2];
                    HM.q{NCoords,1}   = [HM.Angles(i).Name1];
                    HM.q{NCoords+1,1} = [HM.Angles(i).Name2];
                    HM.q{NCoords+2,1} = [HM.Angles(i).Name3];
                    NCoords = NCoords + 3;
                elseif(strcmpi(HM.Angles(i).Joint.Type,'UNI'))
                    HM.Angles(i).PosInq = [NCoords,NCoords+1];
                    HM.q{NCoords,1}   = [HM.Angles(i).Name1];
                    HM.q{NCoords+1,1} = [HM.Angles(i).Name2];
                    NCoords = NCoords + 2;
                elseif(strcmpi(HM.Angles(i).Joint.Type,'REV'))
                    HM.Angles(i).PosInq = NCoords;
                    HM.q{NCoords,1}   = [HM.Angles(i).Name1];
                    NCoords = NCoords + 1;
                else
                    error('Incorrect type of Joint')
                end
            end
        end
        function calcAnglesData(HM,q_t,Deltat,Settings)
            NJoints = size(HM.Joints,1);
            for i=1:NJoints
                if isempty(HM.Joints{i}.Sensor)
                    error(['The model ',HM.ModelName,' does not have associated sensor for the join ',HM.Joints{i}.Name])
                end
                Sensor = HM.Joints{i}.Sensor.getSensorData(q_t);
%                 for j=1:size(Sensor.Val,1)
%                     if abs(Sensor.Val(j,1))< 0.001
%                         Sensor.Val(j,1) = 0;
%                     end
%                     if abs(Sensor.Val(j,2))< 0.001
%                         Sensor.Val(j,2) = 0;
%                     end
%                     if abs(Sensor.Val(j,3))< 0.001
%                         Sensor.Val(j,3) = 0;
%                     end
%                 end
%                 Angles = Sensor.Val*pi/180;
%                 HM.Joints{i}.Angles = Sensor.Val*pi/180;
                if ~isempty(Settings.Smoothing.Method)
                    HM.Joints{i}.Angles = filtTraj(Sensor.Val*pi/180, {'butter';7; 1/Deltat});
                else
                    HM.Joints{i}.Angles = Sensor.Val*pi/180;
                end
                if (max(abs(HM.Joints{i}.Angles(:,1)))) < 0.0001
                    HM.Joints{i}.Angles(:,1) = 0;
                end
                if (max(abs(HM.Joints{i}.Angles(:,2)))) < 0.0001
                    HM.Joints{i}.Angles(:,2) = 0;
                end
                if (max(abs(HM.Joints{i}.Angles(:,3)))) < 0.0001
                    HM.Joints{i}.Angles(:,3) = 0;
                end
%                 [Anglesdot,Angles2dot] = calcVelAcc(Angles,Deltat);
                [HM.Joints{i}.Anglesdot,HM.Joints{i}.Angles2dot] = calcVelAcc(HM.Joints{i}.Angles,Deltat);
                
                
%                 if (max(abs(HM.Joints{i}.Anglesdot(:,1)))) < 0.0001
%                     HM.Joints{i}.Anglesdot(:,1) = 0;
%                 end
%                 if (max(abs(HM.Joints{i}.Anglesdot(:,2)))) < 0.0001
%                     HM.Joints{i}.Anglesdot(:,2) = 0;
%                 end
%                 if (max(abs(HM.Joints{i}.Anglesdot(:,3)))) < 0.0001
%                     HM.Joints{i}.Anglesdot(:,3) = 0;
%                 end
%                 
%                 if (max(abs(HM.Joints{i}.Angles2dot(:,1)))) < 0.0001
%                     HM.Joints{i}.Angles2dot(:,1) = 0;
%                 end
%                 if (max(abs(HM.Joints{i}.Angles2dot(:,2)))) < 0.0001
%                     HM.Joints{i}.Angles2dot(:,2) = 0;
%                 end
%                 if (max(abs(HM.Joints{i}.Angles2dot(:,3)))) < 0.0001
%                     HM.Joints{i}.Angles2dot(:,3) = 0;
%                 end
%                 NFrames = size(HM.Joints{i}.Angles,1);
%                 tend = Deltat * (NFrames-1);
%                 time = 0.0:Deltat:tend;
%                 figure('Name',[HM.Joints{i}.Name,'_Pos1'],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%                 grid on 
%                 plot(time,HM.Joints{i}.Angles(:,1),'r')
%                 plot(time,Angles(:,1),'b')
%                 figure('Name',[HM.Joints{i}.Name,'_Pos2'],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%                 grid on 
%                 plot(time,HM.Joints{i}.Angles(:,2),'r')
%                 plot(time,Angles(:,2),'b')
%                 figure('Name',[HM.Joints{i}.Name,'_Pos3'],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%                 grid on 
%                 plot(time,HM.Joints{i}.Angles(:,3),'r')
%                 plot(time,Angles(:,3),'b')
%                 figure('Name',[HM.Joints{i}.Name,'_Vel1'],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%                 grid on 
%                 plot(time,HM.Joints{i}.Anglesdot(:,1),'r')
% %                 plot(time,Anglesdot(:,1),'b')
%                 figure('Name',[HM.Joints{i}.Name,'_Vel2'],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%                 grid on 
%                 plot(time,HM.Joints{i}.Anglesdot(:,2),'r')
% %                 plot(time,Anglesdot(:,2),'b')
%                 figure('Name',[HM.Joints{i}.Name,'_Vel3'],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%                 grid on 
%                 plot(time,HM.Joints{i}.Anglesdot(:,3),'r')
% %                 plot(time,Anglesdot(:,3),'b')
%                 figure('Name',[HM.Joints{i}.Name,'_Acc1'],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%                 grid on 
%                 plot(time,HM.Joints{i}.Angles2dot(:,1),'r')
% %                 plot(time,Angles2dot(:,1),'b')
%                 figure('Name',[HM.Joints{i}.Name,'_Acc2'],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%                 grid on 
%                 plot(time,HM.Joints{i}.Angles2dot(:,2),'r')
% %                 plot(time,Angles2dot(:,2),'b')
%                 figure('Name',[HM.Joints{i}.Name,'_Acc3'],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%                 grid on 
%                 plot(time,HM.Joints{i}.Angles2dot(:,3),'r')
% %                 plot(time,Angles2dot(:,3),'b')
%                 HM.Joints{i}.Angles2dot = zeros(size(Sensor.Val,1),3);
            end
        end
        function EndSegments = getEndSegments(HM)
            NSeg = size(HM.Segments,1);
            EndSegments = [];
            for i=1:NSeg
                if isempty(HM.Segments(i).Distals)
                    EndSegments = [EndSegments;HM.Segments(i)];
                end
            end
        end
        function InitSegment = getInitSegment(HM)
            NSeg = size(HM.Segments,1);
            InitSegment = [];
            for i=1:NSeg
                if isempty(HM.Segments(i).Proximal)
                    InitSegment = [InitSegment;HM.Segments(i)];
                end
            end
        end
        function getLocalPointsCoords(HM)
            NSegments = size(HM.Segments,1);
            for i=1:NSegments
                % Not include Ground, there are not markers
                if HM.Segments(i).Fixed ~= 1 % If is not Ground                    
                %if ~strcmpi(HM.Segments(i).Name,'Ground')
                    HM.Segments(i).getLocalPointsCoords();
%                     HM.Segments(i).getGraphicData();
                end
            end
            
        end
        function setModelHierarchy(HM)
            NJoints = size(HM.Joints,1);
            NSegments = size(HM.Segments,1);
            for i=1:NJoints
                ParentName = HM.Joints{i}.Seg1;
                ChildName  = HM.Joints{i}.Seg2;
                ParentIndex = getVecIndex(ParentName,HM.Segments);
                ChildIndex  = getVecIndex(ChildName,HM.Segments);
                HM.Segments(ParentIndex).addDistal(HM.Segments(ChildIndex),HM.Joints{i});
                HM.Segments(ChildIndex).addProximal(HM.Segments(ParentIndex),HM.Joints{i});
            end
            for i=1:NSegments
                if size(HM.Segments(i).Distals)>1
                    HM.Segments(i).Root = 1;
                end
            end
        end
        function Phi = mkAdditionalCtrs(HM)
            NSegments = size(HM.Segments,1);
            Phi = {};
            for i=1:NSegments
                if(HM.Segments(i).Fixed == 0)
                    Phi = [Phi;HM.Segments(i).mkAdditionalCtrs()];
                end
            end
            HM.Phi = [HM.Phi;Phi];
        end
        function mkModelCtrs(HM)
            HM.mkRigidSegmentCtrs();
%             HM.mkRelPointCtrs();
            HM.mkJointCtrs();
            HM.mkAngleCtrs();
            % Make char Phi and Phiq
            % initialize time counter
            %c=tic;
            
            verChar = version('-release');
            verNum  = str2num(verChar(1:4));
            
            if verNum >= 2020
                % release 2020a produce a CLEAN output for function jacobian()
                HM.mkPhiq_2020a();
            else            
                % release 2019a and older produce a DIRTY output for function jacobian()
                HM.mkPhiq();
            end
            
            NPhi = size(HM.Phi,1);
            for i=1:NPhi
                Phi_tmp{i,1} = char(HM.Phi(i));
            end
            HM.Phi = Phi_tmp;
           %[H, MI, S] = second2HMS(toc(c));
           % printElapsedTime(H, MI, S, 3, 'char time: ');
        end
        function mkAngleCtrs(HM)
            NAngles  = size(HM.Angles,1);
            for i = 1:NAngles
                if(strcmpi(HM.Angles(i).Joint.Type,'SPH'))
                    HM.Phi = [HM.Phi; mkCtrSPHAngle(HM.Angles(i))];
                elseif(strcmpi(HM.Angles(i).Joint.Type,'UNI'))
                    HM.Phi = [HM.Phi; mkCtrUNIAngle(HM.Angles(i))];
                elseif(strcmpi(HM.Angles(i).Joint.Type,'REV'))
                    HM.Phi = [HM.Phi; mkCtrREVAngle(HM.Angles(i))];
                end
            end
        end
        function mkJointCtrs(HM)
            NJoints = size(HM.Joints,1);
            for i= 1:NJoints
                
                if strcmpi(HM.Joints{i,1}.Type,'UNI')
                    if (HM.Joints{i,1}.CAngle == pi/2)
                        HM.Phi = [HM.Phi; mkCtrUjoint(HM.Joints{i,1}.Vector1,HM.Joints{i,1}.Vector2)];
                    else
                        HM.Phi = [HM.Phi; mkCtrUjoint(HM.Joints{i,1}.Vector1,HM.Joints{i,1}.Vector2,HM.Joints{i,1}.CAngle)];
                    end
%                 elseif strcmpi(HM.Joints{i,1}.Type,'SPH')
%                     if ~strcmpi(HM.Joints{i,1}.Point1.Point.Name,HM.Joints{i,1}.Point2.Point.Name)
%                         
%                     end
                elseif strcmpi(HM.Joints{i,1}.Type,'REV')
                    if ~strcmpi(HM.Joints{i,1}.Vector1.Vector.Name,HM.Joints{i,1}.Vector2.Vector.Name)
                        HM.Phi = [HM.Phi; mkCtrRevJoint(HM.Joints{i,1}.Vector1,HM.Joints{i,1}.Vector2)];
                    end
                end
            end
        end
        function mkRelPointCtrs(HM)
%             NSegments = size(HM.Segments,1);
%             for i = 1:NSegments
%                 HM.Phi = [HM.Phi; HM.Segments(i).mkRelPointCtr()];
%             end
        end
        function mkRigidSegmentCtrs(HM)
            % initialize time counter
%             b=tic;
            HM.setTypeOfBases();
            NSegments = size(HM.Segments,1);
            for i = 1:NSegments
                % Test if is fixed segment
                if (HM.Segments(i).Fixed)== 0 % It is not Fixed (Ground)
                  HM.Phi = [HM.Phi; HM.Segments(i).mkCtrRigidSegment()];
                end
            end
            % display simulation time
%             [H, MI, S] = second2HMS(toc(b));
%             printElapsedTime(H, MI, S, 3, 'constraint time: ');

        end
        function Phi = mkMarkerCtrs(HM)
            NSegments = size(HM.Segments,1);
            Phi = {};
            for i=1:NSegments
                if(HM.Segments(i).Fixed == 0)
                    Phi = [Phi;HM.Segments(i).mkMarkerCtrs()];
                end
            end
            HM.Phi = [HM.Phi;Phi];
        end
        function calcMassMatrices(HM)
           NSegments = size(HM.Segments,1); 
           for i=1:NSegments
                if(HM.Segments(i).Fixed == 0)
                    J=[5 7 8; 12 4 5; 4 9 12];
                    HM.Segments(i).calcMassMatrix(J);
                end
            end
        end
        function mkPhiq(HM)
            % Convert symbolic jacobian to text
            PhiqChar_tmp  = char(jacobian(HM.Phi,sym(HM.q,'real')));
            PhiqChar_tmp2 = PhiqChar_tmp(10:end-1);  % Remove 'matrix([[' in first and ')' at the end.
            PhiqChar      = strread(PhiqChar_tmp2, '%s', 'delimiter', ',[]'); % Return vector with Phiq values and two empty values inter rows

            % Convert text to sparse format
            i= 1; % row
            j= 1; % col
            
            % Variables to define Sparse matrix
            HM.Phiq.rows = []; % row index vector of non-zero elements
            HM.Phiq.cols = []; % col index vector of non-zero elements
            HM.Phiq.s = {};
            
            k = 1;
            nElements = length(PhiqChar);
            while k <= nElements
                
                %If element is not zero store position and expression
                if ~strcmp(PhiqChar{k},'0') && ~isempty(PhiqChar{k})
                    HM.Phiq.rows = [HM.Phiq.rows; i];
                    HM.Phiq.cols = [HM.Phiq.cols; j];
                    HM.Phiq.s =[HM.Phiq.s; PhiqChar{k}];
                end
                
                j=j+1; % Increment column counter
                
                if isempty(PhiqChar{k}) % If true start of new row, reset counters
                    k = k+1; % needed because two consecutive white space are generated by strread
                    i = i+1;
                    j = 1;
                end
                
                k = k+1; % Increment while counter
            end
        end
        function mkPhiq_2020a(HM)
            % Convert symbolic jacobian to text
            PhiqChar_tmp  = char(jacobian(HM.Phi,sym(HM.q,'real')));            
            % Matlab 2020a (and later) produce a clean output for jacobian
            
            % Separo la jacobiana en filas
            PhiqChar_tmp2 = strread(PhiqChar_tmp(2:end-1), '%s', 'delimiter', ';');
            % PhiqChar_tmp2 es un cell columna y cada elemento contiene toda una fila de la jacobiana de Phi
            
            % Vamos a transforma en un cell nRows x nCols donde cada elemento es un elemento de la jacobiana
            nRows = length(HM.Phi);
            nCols = length(HM.q);
            PhiqChar = cell(nRows, nCols);
            for i = 1:nRows
                PhiqChar(i,:) = strread(PhiqChar_tmp2{i}, '%s', 'delimiter', ',')';
            end
            
            % Variables to define Sparse matrix
            HM.Phiq.rows = []; % row index vector of non-zero elements
            HM.Phiq.cols = []; % col index vector of non-zero elements
            HM.Phiq.s = {};            
            
            for i= 1:nRows
                for j= 1:nCols                    
                    if ~strcmp(PhiqChar{i,j},'0')                     
                        HM.Phiq.rows = [HM.Phiq.rows; i];
                        HM.Phiq.cols = [HM.Phiq.cols; j];
                        HM.Phiq.s =[HM.Phiq.s; PhiqChar{i,j}];
                    end                
                end
            end
        end
        function mkPointPhiq(HM,PhiPoint)
            NPhi = size(HM.Phi,1);
            Nq   = size(HM.q,1);
            NPhiPoint = size(PhiPoint,1);
            NIniPhi = NPhi-NPhiPoint;
            % Type of Phi
            % A-B-(N1)*C-(N2)*D-(N3)*(E-B)
            %   where:
            %      A is a 3x1 vector
            for i=1:NPhiPoint
                for j=1:Nq
                    Posq = findstr(HM.q{j},PhiPoint{i});
                    PosLoc = findstr(['_',HM.q{j}],PhiPoint{i});
                    Numq_Loc = size(Posq,2);
                    Num_Loc  = size(PosLoc,2);
                    Numq = Numq_Loc-Num_Loc;
                    PosMulti = findstr(')* (',PhiPoint{i});
                    if Numq>0
                        % q is in the Phi{i} constraint
                        if Numq == 1
                            % Choose Pos in q and not Local Point
                            for k=1:Numq_Loc
                                if k==Numq_Loc
                                    Posq = Posq(k);
                                elseif Posq(k)-PosLoc(k) ~= 1
                                    Posq = Posq(k);
                                    break;
                                end
                            end
                            % Calculate the derivative of A = 1
                            if(Posq == 1)
                                HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                HM.Phiq.cols = [HM.Phiq.cols; j];
                                HM.Phiq.s =[HM.Phiq.s; '1'];
                                % Calculate the derivative of C = -N1
                            elseif (Posq-PosMulti(1))==4
                                MulSigMul = PhiPoint{i}(1:PosMulti(1));
                                PosIntMinus  = findstr('- (',MulSigMul);
                                IniN = PosIntMinus;
                                HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                HM.Phiq.cols = [HM.Phiq.cols; j];
                                HM.Phiq.s =[HM.Phiq.s; PhiPoint{i}(IniN:(PosMulti(1)))];
                                % Calculate the derivative of D = -N2
                            elseif (Posq-PosMulti(2))==4
                                MulSigMul = PhiPoint{i}(PosMulti(1)+1:PosMulti(2));
                                PosIntMinus  = findstr('- (',MulSigMul);
                                IniN = PosMulti(1)+PosIntMinus;
                                HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                HM.Phiq.cols = [HM.Phiq.cols; j];
                                HM.Phiq.s =[HM.Phiq.s;PhiPoint{i}(IniN:(PosMulti(2)))]; 
                                % Calculate the derivative of E = -N3
                            elseif (Posq-PosMulti(3))==4
                                MulSigMul = PhiPoint{i}(PosMulti(2)+1:PosMulti(3));
                                PosIntMinus  = findstr('- (',MulSigMul);
                                IniN = PosMulti(2)+PosIntMinus;
                                HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                HM.Phiq.cols = [HM.Phiq.cols; j];
                                HM.Phiq.s =[HM.Phiq.s;PhiPoint{i}(IniN:(PosMulti(3)))];
                            else
%                                 error('Diferent type of equation...')
                                  warning(['Phi is: ',PhiPoint{i},' and q is: ',HM.q{j}])
                            end
                            % Calculate the derivative of B = -1+N3
                        elseif Numq == 2
                            MulSigMul = PhiPoint{i}(PosMulti(2)+1:PosMulti(3));
                            PosIntMinus  = findstr('- (',MulSigMul);
                            IniN = PosMulti(2)+PosIntMinus+2;
                            HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                            HM.Phiq.cols = [HM.Phiq.cols; j];
                            HM.Phiq.s =[HM.Phiq.s;['-1+',PhiPoint{i}(IniN:(PosMulti(3)))]];
                        else
%                             error('Diferent type of equation...')
                              warning(['Phi is: ',PhiPoint{i},' and q is: ',HM.q{j}])
                        end
                    end
                end
            end
        end
        function mkVectorPhiq(HM,PhiVector)
            NPhi = size(HM.Phi,1);
            Nq   = size(HM.q,1);
            NPhiVector = size(PhiVector,1);
            NIniPhi = NPhi-NPhiVector;
            for i=1:NPhiVector
                for j=1:Nq
                    Posq = findstr(HM.q{j},PhiVector{i});
                    PosLoc = findstr(['_',HM.q{j}],PhiVector{i});
                    Numq_Loc = size(Posq,2);
                    Num_Loc  = size(PosLoc,2);
                    Numq = Numq_Loc-Num_Loc;
                    PosMulti = findstr(')* (',PhiVector{i});
                    PosMinus = findstr(' - ',PhiVector{i});
                    % q is in the Phi{i} constraint
                    if Numq>0
                         % Choose Pos in q and not Local Point
                         for k=1:Numq_Loc
                             if k==Numq_Loc
                                 Posq = Posq(k);
                             elseif Posq(k)-PosLoc(k) ~= 1
                                 Posq = Posq(k);
                                 break;
                             end
                         end
                        % Type of Phi
                        % A(-)N1*B(-)N2*C(-)N3*(D-E)
                        if (size(PosMulti,2) == 3)
                            % q is in the Phi{i} constraint
                            % Calculate the derivative of A = 1
                            if(Posq == 1)
                                HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                HM.Phiq.cols = [HM.Phiq.cols; j];
                                HM.Phiq.s =[HM.Phiq.s; '1'];
                            % Calculate the derivative of B = -N1
                            elseif (Posq-PosMulti(1))==4
                                MulSigMul = PhiVector{i}(1:PosMulti(1));
                                PosIntMinus  = findstr('- (',MulSigMul);
                                IniN = PosIntMinus;
                                HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                HM.Phiq.cols = [HM.Phiq.cols; j];
                                HM.Phiq.s =[HM.Phiq.s; PhiVector{i}(IniN:(PosMulti(1)))];
                            % Calculate the derivative of C = -N2
                            elseif (Posq-PosMulti(2))==4
                                MulSigMul = PhiVector{i}(PosMulti(1)+1:PosMulti(2));
                                PosIntMinus  = findstr('- (',MulSigMul);
                                IniN = PosMulti(1)+PosIntMinus;
                                HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                HM.Phiq.cols = [HM.Phiq.cols; j];
                                HM.Phiq.s =[HM.Phiq.s;PhiVector{i}(IniN:(PosMulti(2)))];
                            % Calculate the derivative of D = -N3
                            elseif (Posq-PosMulti(3))==5
                                MulSigMul = PhiVector{i}(PosMulti(2)+1:PosMulti(3));
                                PosIntMinus  = findstr('- (',MulSigMul);
                                IniN = PosMulti(2)+PosIntMinus;
                                HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                HM.Phiq.cols = [HM.Phiq.cols; j];
                                HM.Phiq.s =[HM.Phiq.s;PhiVector{i}(IniN:(PosMulti(3)))];
                            % Calculate the derivative of E = +N3
                            elseif (Posq-PosMulti(3))>5
                                MulSigMul = PhiVector{i}(PosMulti(2)+1:PosMulti(3));
                                PosIntMinus  = findstr('- (',MulSigMul);
                                IniN = PosMulti(2)+PosIntMinus+2;
                                HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                HM.Phiq.cols = [HM.Phiq.cols; j];
                                HM.Phiq.s =[HM.Phiq.s;PhiVector{i}(IniN:(PosMulti(3)))];
                            end
                        elseif (size(PosMulti,2) == 2)
                            % % Type of Phi
                            % A(-)N1*B(-)N2*C
                            if isempty(PosMinus)
                                % Calculate the derivative of A = 1
                                if(Posq == 1)
                                    HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                    HM.Phiq.cols = [HM.Phiq.cols; j];
                                    HM.Phiq.s =[HM.Phiq.s; '1'];
                                % Calculate the derivative of B = (-)N1
                                elseif (Posq-PosMulti(1))==4
                                    MulSigMul = PhiVector{i}(1:PosMulti(1));
                                    PosIntMinus  = findstr('- (',MulSigMul);
                                    IniN = PosIntMinus;
                                    HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                    HM.Phiq.cols = [HM.Phiq.cols; j];
                                    HM.Phiq.s =[HM.Phiq.s; PhiVector{i}(IniN:(PosMulti(1)))];
                                % Calculate the derivative of C = (-)N2
                                elseif (Posq-PosMulti(2))==4
                                    MulSigMul = PhiVector{i}(PosMulti(1)+1:PosMulti(2));
                                    PosIntMinus  = findstr('- (',MulSigMul);
                                    IniN = PosMulti(1)+PosIntMinus;
                                    HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                    HM.Phiq.cols = [HM.Phiq.cols; j];
                                    HM.Phiq.s =[HM.Phiq.s;PhiVector{i}(IniN:(PosMulti(2)))];
                                end
                            else
                            % Type of Phi
                            % A(-)N1*B(-)N2*(C-D)
                                % Calculate the derivative of A = 1
                                if(Posq == 1)
                                    HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                    HM.Phiq.cols = [HM.Phiq.cols; j];
                                    HM.Phiq.s =[HM.Phiq.s; '1'];
                                % Calculate the derivative of B = (-)N1
                                elseif (Posq-PosMulti(1))==4
                                    MulSigMul = PhiVector{i}(1:PosMulti(1));
                                    PosIntMinus  = findstr('- (',MulSigMul);
                                    IniN = PosIntMinus;
                                    HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                    HM.Phiq.cols = [HM.Phiq.cols; j];
                                    HM.Phiq.s =[HM.Phiq.s; PhiVector{i}(IniN:(PosMulti(1)))];
                                % Calculate the derivative of C = (-)N2
                                elseif (Posq-PosMulti(2))==5
                                    MulSigMul = PhiVector{i}(PosMulti(1)+1:PosMulti(2));
                                    PosIntMinus  = findstr('- (',MulSigMul);
                                    IniN = PosMulti(1)+PosIntMinus;
                                    HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                    HM.Phiq.cols = [HM.Phiq.cols; j];
                                    HM.Phiq.s =[HM.Phiq.s;PhiVector{i}(IniN:(PosMulti(2)))];
                                % Calculate the derivative of E = (+)N2
                                elseif (Posq-PosMulti(2))>5
                                    MulSigMul = PhiVector{i}(PosMulti(1)+1:PosMulti(2));
                                    PosIntMinus  = findstr('- (',MulSigMul);
                                        IniN = PosMulti(1)+PosIntMinus+2;
                                        HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                        HM.Phiq.cols = [HM.Phiq.cols; j];
                                        HM.Phiq.s =[HM.Phiq.s;PhiVector{i}(IniN:(PosMulti(2)))];
                                end
                            end
                        elseif (size(PosMulti,2) == 1)
                           % Type of Phi
                           % A(-)N1*B
                           if isempty(PosMinus)
                               % Calculate the derivative of A = 1
                               if(Posq == 1)
                                   HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                   HM.Phiq.cols = [HM.Phiq.cols; j];
                                   HM.Phiq.s =[HM.Phiq.s; '1'];
                                   % Calculate the derivative of B = (-)N1
                               elseif (Posq-PosMulti(1))==4
                                   MulSigMul = PhiVector{i}(1:PosMulti(1));
                                   PosIntMinus  = findstr('- (',MulSigMul);
                                   IniN = PosIntMinus;
                                   HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                   HM.Phiq.cols = [HM.Phiq.cols; j];
                                   HM.Phiq.s =[HM.Phiq.s; PhiVector{i}(IniN:(PosMulti(1)))];
                               end
                           % Type of Phi
                           % A(-)N1*(B-C)
                           else
                               % Calculate the derivative of A = 1
                               if(Posq == 1)
                                   HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                   HM.Phiq.cols = [HM.Phiq.cols; j];
                                   HM.Phiq.s =[HM.Phiq.s; '1'];
                                % Calculate the derivative of B = (-)N1
                                elseif (Posq-PosMulti(1))==5
                                    MulSigMul = PhiVector{i}(1:PosMulti(1));
                                    PosIntMinus  = findstr('- (',MulSigMul);
                                    IniN = PosMulti(1)+PosIntMinus;
                                    HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                    HM.Phiq.cols = [HM.Phiq.cols; j];
                                    HM.Phiq.s =[HM.Phiq.s;PhiVector{i}(IniN:(PosMulti(1)))];
                                % Calculate the derivative of E = (+)N3
                                elseif (Posq-PosMulti(1))>5
                                    MulSigMul = PhiVector{i}(1:PosMulti(1));
                                    PosIntMinus  = findstr('- (',MulSigMul);
                                    IniN = PosMulti(1)+PosIntMinus+2;
                                    HM.Phiq.rows = [HM.Phiq.rows; i+NIniPhi];
                                    HM.Phiq.cols = [HM.Phiq.cols; j];
                                    HM.Phiq.s =[HM.Phiq.s;PhiVector{i}(IniN:(PosMulti(1)))];
                                end
                           end
                        else
                            error('Incorrect Type of equation');
                        end
                    end
                end
            end
                        
        
            
        end
        function parseAlmFile(HM,Path,File)
            % -------------------------------------------------------------------------------------------------------------------
            % read the file
            % -------------------------------------------------------------------------------------------------------------------
            % Initialize counters
            FileRow = 1;  % row number of the file where the function is reading the file
            % Open the calibration file for reading
            FidCalib = fopen([Path,File]);
            if FidCalib == -1
                error(['The file ',File,' can not be opened in path (',Path,'). ']);
            end
            Cont = 0;
            NumPar = 2000;
            VAxis = 0;
            while feof(FidCalib) == 0
                % Get the whole line
                Line = fgets(FidCalib);
                if (FileRow >1 && FileRow<NumPar)
                    % All the subject parameters
                    if strcmp(deblank(Line),'ENVIRONMENTAL_PARAMETER')
                        NumPar = FileRow;
                        VAxis = 1;
                    else
                        [Acronym, Remainder] = strtok(Line,',');
                        [Value,Remainder]    = strtok(Remainder,',');
                        [Keyword,Remainder]  = strtok(Remainder,',');
                        Cont = Cont +1;
                        HM.SubjectParameters(Cont).Keyword = deblank(Keyword);
                        HM.SubjectParameters(Cont).Value   = Value;
                        HM.SubjectParameters(Cont).Acronym = Acronym;
                        if strcmp(deblank(Keyword),'Weight') %|| strcmp(Acronym,'AM25')
                            HM.Mass = str2num(Value);
                        elseif strcmp(deblank(Keyword),'Gender') %|| strcmp(Acronym,'AM33')
                            HM.Gender = deblank(Value);
                        end
                    end
                elseif FileRow >NumPar && VAxis == 1
                    [Keyword, Remainder] = strtok(Line,',');
                    if strcmp(Keyword,'Measurement Vertical Axis')
                        HM.MeasVAxis = str2num(Remainder)';
                    end
                    if strcmp(deblank(Line),'MARKERS')
                        MarkerTable = 1;
                        ALTable = 0;
                        KeyPosture = 0;
                        VAxis = 0;
                    end
                elseif FileRow >NumPar && MarkerTable == 1
                    % Complete the marker table. For each segment put all the Markers belongig to it
                    [Keyword, Remainder] = strtok(Line,',');
                    if strcmpi(Keyword,'SEGMENT')
                        % When the table is completed pass to the AL Table and store the marker position in M.
                        MarkerTable = 0;
                        ALTable =1;
                        Line = Line(29:end);
                        Pos = strfind(Line,',,,');
                        NMarkers = size(Pos,2);
                        M.Name = Line(2:Pos(1)-1);
                        for i=1:NMarkers-1
                            Mtmp.Name = Line(Pos(i)+3:Pos(i+1)-1);
                            M =[M;Mtmp];
                        end
                        SegNameTmp='A';
                        j=0;
                    elseif strcmp(deblank(Line),'ANATOMICAL_LANDMARKS')
                        % We only want to pass the line
                    else
                        % Create segments and markers belongin to them.
                        k = strfind(Line,',');
                        NMarkers = size(k,2);
                        SegmentName = Line(1:k(1)-1);
                        HM.addSegment(SegmentName);
                        SegmentPos = getVecIndex(SegmentName,HM.Segments);
                        for i=1:NMarkers-1
                            MarkerName = Line(k(i)+1:k(i+1)-1);
                            Marker = HM.addMarker(MarkerName);
                            HM.Segments(SegmentPos).addMarker(Marker);
                        end
%                         MarkerName = Line(k(NMarkers)+1:end-2);
%                         MarkerName = Line(k(NMarkers)+1:end-1);
                        MarkerName = deblank(Line(k(NMarkers)+1:end));
                        Marker = HM.addMarker(MarkerName);
                        HM.Segments(SegmentPos).addMarker(Marker);
                    end
                elseif FileRow > NumPar && ALTable ==1
                    % Put the global position of the AL. And global positions of marker in the segment
                    if strcmpi(Line(1:7),'POSTURE')
                        % When AL Table is completed pass to posture line
                        ALTable =0;
                        KeyPosture = 1;
                        Head =1;
                    else
                        % Obtain the AL globlal position and add AL to the segment
                        k = strfind(Line,',');
                        SegmentName = Line(1:k(1)-1);
                        ALName = Line(k(1)+1:k(2)-1);
                        ALx = str2double(Line(k(2)+1:k(3)-1));
                        ALy = str2double(Line(k(3)+1:k(4)-1));
                        ALz = str2double(Line(k(4)+1:k(5)-1));
                        AL = HM.addAL(ALName);
                        SegmentPos = getVecIndex(SegmentName,HM.Segments);
                        HM.Segments(SegmentPos).addAL(AL);
                        ALPos = getVecIndex(ALName,HM.ALs);
                        HM.ALs(ALPos).GlobalCoord = [ALx;ALy;ALz];
                        NLocalMarkers = size(HM.Segments(SegmentPos).LocalMarkers,1);
                        if strcmpi(SegmentName,SegNameTmp)
                            j = j+1;
                        else
                            j=1;
                        end
                        SegNameTmp = SegmentName;
                        % For each AL the markers have diferent measured positions
                        for i=1:NLocalMarkers
                            MarkerName = HM.Segments(SegmentPos).LocalMarkers(i).Point.Name;
                            MarkerPos = getVecIndex(MarkerName,M);
                            Mx = str2double(Line(k((3*MarkerPos-2)+4)+1:k((3*MarkerPos-2)+5)-1));
                            My = str2double(Line(k((3*MarkerPos-1)+4)+1:k((3*MarkerPos-1)+5)-1));
                            % To add the last marker position
                            if ((3*MarkerPos)+5) > size(k,2)
                                Mz = str2double(Line(k((3*MarkerPos)+4)+1:end));
                            else
                                Mz = str2double(Line(k((3*MarkerPos)+4)+1:k((3*MarkerPos)+5)-1));
                            end
                            HM.Segments(SegmentPos).LocalMarkers(i).Point.MeasuredCoord(:,j) = [Mx;My;Mz];
                        end
                        % To add the global position of the parent segment
                        if strcmpi (SegmentName,'LeftShank')|| strcmpi (SegmentName,'LeftFoot')
                            if strcmpi (SegmentName,'LeftShank')
                                PreSegPos = getVecIndex('LeftThigh',HM.Segments);
                            elseif strcmpi (SegmentName,'LeftFoot')
                                PreSegPos = getVecIndex('LeftShank',HM.Segments);
                            else
                                error ('Incorrect segment Name')
                            end
                            NPreSegMarkers = size(HM.Segments(PreSegPos).LocalMarkers,1);
                            for i=1:NPreSegMarkers
                                MarkerName = HM.Segments(PreSegPos).LocalMarkers(i).Point.Name;
                                MarkerPos = getVecIndex(MarkerName,M);
                                Mx = str2double(Line(k((3*MarkerPos-2)+4)+1:k((3*MarkerPos-2)+5)-1));
                                My = str2double(Line(k((3*MarkerPos-1)+4)+1:k((3*MarkerPos-1)+5)-1));
                                % To add the last marker position
                                if ((3*MarkerPos)+5) > size(k,2)
                                    Mz = str2double(Line(k((3*MarkerPos)+4)+1:end));
                                else
                                    Mz = str2double(Line(k((3*MarkerPos)+4)+1:k((3*MarkerPos)+5)-1));
                                end
                                HM.Segments(SegmentPos).PreSegMarkPos(:,i) = [Mx;My;Mz];
                            end
                        end
                    end
                elseif FileRow > NumPar && KeyPosture ==1
                    % Store Posture position in markers global positions.
                    if Head ==1
                        Head = 0;
                    else
                        [PosturName, Remainder] = strtok(Line,',');
                        Posture.Name = PosturName;
                        k = strfind(Remainder,',');
                        NMarkers = size(HM.Markers,1);
                        for i=1:NMarkers
                            MarkerName = HM.Markers(i).Name;
                            MarkerPos = getVecIndex(MarkerName,M);
                            Mx = str2double(Remainder(k(3*MarkerPos-2)+1:k((3*MarkerPos-2)+1)-1));
                            My = str2double(Remainder(k(3*MarkerPos-1)+1:k((3*MarkerPos-1)+1)-1));
%                             if MarkerPos == NMarkers                            
                            if 3*MarkerPos == size(k,2)
                                Mz = str2double(Remainder(k(3*MarkerPos)+1:end));
                            else
                                Mz = str2double(Remainder(k(3*MarkerPos)+1:k((3*MarkerPos)+1)-1));
                            end
                            Posture.Glob = [Mx;My;Mz];
                            HM.Markers(i).Postures = [HM.Markers(i).Postures;Posture]; 
                        end
                    end
                end
                FileRow = FileRow + 1;
            end
            % Check the Mass and gender
            if isempty(HM.Mass)
                error('Incorrect keyword for the weight');
            end
            if isempty(HM.Gender)
                error('Incorrect keyword for the gender');
            end
            % Pass the gender to the segments
            NSegments = size(HM.Segments,1);
            for i=1:NSegments
                HM.Segments(i).Gender = HM.Gender;
            end
            fclose(FidCalib);
        end
        function parseExpFile(HM,Path,File)
            if strcmpi(File(end-2:end),'exp')
                FidExp = fopen([Path,File]);
                SecMI = 0;
                while feof(FidExp) == 0
                    % Get the whole line
                    Line = fgets(FidExp);
                    if strcmp(deblank(Line),'SECTION MOMENT_OF_INERTIA')
                        SecMI = 1;
                    end 
                end
                fclose(FidExp);
                if SecMI == 0
                    FidPSP = fopen([Path,File(1:end-3),'rsp']);
                    while feof(FidPSP) == 0
                       % Get the whole line
                       Line = fgets(FidPSP);
                       if strcmp(deblank(Line),'SECTION MOMENT_OF_INERTIA')
                           LineMI{1} = Line;
                           RowMI = 1;
                           SecMI = 1;
                       elseif SecMI == 1;
                           RowMI = RowMI +1;
                           LineMI{RowMI} = Line;
                       end 
                    end
                    fclose(FidPSP);
                    FidExp = fopen([Path,File],'a');
                    NRowMI = size(LineMI,2);
                    fprintf(FidExp,'\n');
                    for i=1:NRowMI
                        fprintf(FidExp,'%s', LineMI{1,i});
                    end
                    fclose(FidExp);
                end
            end
            ParserExp = HUMAN_PARSER_EXP();
            ParserExp.readExp(Path,File);
            ParserExp.checkPointsData(HM.Points,File);
            ParserExp.parsePoints(HM.Points);
            ParserExp.parseVectors(HM.Vectors);
            ParserExp.parseAngles(HM.Angles);
            ParserExp.parseXMat(HM);
            ParserExp.parseMarkers(HM);
            ParserExp.parseSegmentsParameters(HM);
            HM.getLocalPointsCoords();
%             Escribir_puntos_locales;
        end
        function parseMatFile(HM,PathMat,File)
            ParserMat = HUMAN_PARSER_MAT();
            ParserMat.readMat(PathMat,File);
            ParserMat.parsePoints(HM.Segments);
            ParserMat.parseMarkers(HM);
            ParserMat.parseXMat(HM);
        end
        function parseXmlModelFile(HM,Path,File) 
            ParserXml = HUMAN_PARSER_XML();
            ParserXml.readxml(Path,File)
            ParserXml.parseName(HM);
            ParserXml.parseSegments(HM,Path); 
            ParserXml.parseJoints(HM);
            ParserXml.parseSensors(HM);
            ParserXml.parseAngles(HM);
            
        end
        function parseSubjectPars(HM,Path,File)
            if strcmpi(File((end-2):end),'EXP')||strcmpi(File((end-2):end),'rsp')
                HM.parseExpFile(Path,File);
            elseif strcmpi(File((end-2):end),'mat')
                HM.parseMatFile(Path,File);
            elseif strcmpi(File((end-2):end),'psp')
                HM.parsePspFile(Path,File); 
            elseif strcmpi(File((end-2):end),'par')
                HM.parseParFile(Path,File);
            end
            
        end
        function parseParFile(HM,PathPar,FilePar)
            ParsePar = HUMAN_PARSER_PAR(PathPar,FilePar);
            ParsePar.parseMarkers(HM);
            NSegments = size(HM.Segments,1);
            NAngles  = size(HM.Angles,1);
            for i=1:NSegments
                NVectors = size(HM.Segments(i).LocalVectors,1);
                NPoints  = size(HM.Segments(i).LocalPoints,1);
                NMarkers = size(HM.Segments(i).LocalMarkers,1);
                for j=1:NVectors
                    HM.Segments(i).LocalVectors(j).Vector.GlobalCoord = HM.Segments(i).LocalVectors(j).LocCoord;
                end
                for j=1:NPoints
                    if HM.Segments(i).Fixed == 1 % If is Ground
                    %if strcmpi(HM.Segments(i).Name,'Ground')
                        HM.Segments(i).LocalPoints(j).Point.GlobalCoord = HM.Segments(i).LocalPoints(j).LocCoord;
                    elseif strcmpi(HM.Segments(i).Name,'Pelvis')
                        HM.Segments(i).LocalPoints(1).Point.GlobalCoord = HM.Segments(i).LocalPoints(1).LocCoord;
                        HM.Segments(i).LocalPoints(j).Point.GlobalCoord = HM.Segments(i).LocalPoints(j).LocCoord + HM.Segments(i).LocalPoints(1).Point.GlobalCoord;
                    else
                        HM.Segments(i).LocalPoints(j).Point.GlobalCoord = HM.Segments(i).LocalPoints(j).LocCoord + HM.Segments(i).LocalPoints(1).Point.GlobalCoord;
                    end
                end 
                for j=1:NMarkers
                    if HM.Segments(i).Fixed == 1 % If is Ground
                    %if strcmpi(HM.Segments(i).Name,'Ground')
                        HM.Segments(i).LocalMarkers(j).Point.GlobalCoord = HM.Segments(i).LocalPoints(j).LocCoord;
                    else
                        HM.Segments(i).LocalMarkers(j).Point.GlobalCoord = HM.Segments(i).LocalMarkers(j).LocCoord + HM.Segments(i).LocalPoints(1).Point.GlobalCoord;
                    end
                end 
            end
            for i=1:NAngles
                HM.Angles(i).a1 = 0;
                HM.Angles(i).a2 = 0;
                HM.Angles(i).a3 = 0;
            end
        end
        function parsePspFile(HM,Path,File)
            ParsePsp = HUMAN_PARSER_PSP();
            ParsePsp.readxml(Path,File);
            ParsePsp.parseSegmentsInfo(HM);
            ParsePsp.parseSensor(HM);
            HM.checkSubjectParsData(File);
            ParsePsp.parseXMat(HM);
        end
        function setSubAdditPar(HM,SubAdditPar)
            HM.SubAdditPar = SubAdditPar;
            NSegments = size(HM.Segments,1);
            for i=1:NSegments
                HM.Segments(i).SubAdditPar = SubAdditPar;
            end
        end
        function setTypeOfBases(HM)
            NSegments = size(HM.Segments,1);
            for i = 1:NSegments
                if HM.Segments(i).Fixed ~= 1 % If is not Fixed(Ground)
                    HM.Segments(i).setTypeOfBasis;
                end
            end
        end
        function setGraphicPars(HM)
            NSegments = size(HM.Segments,1);
            for i=1:NSegments
                HM.Segments(i).setGraphicPars();
            end
        end
        function writeRSP(HM,FileName,Path)
            NSegments = size(HM.Segments,1);
            FileNameExt = [FileName,'.rsp'];
            
            % get segments index
            PelvisIndex = getVecIndex('Pelvis',HM.Segments);
            RThighIndex = getVecIndex('RightThigh',HM.Segments);
            LThighIndex = getVecIndex('LeftThigh',HM.Segments);
            RShankIndex = getVecIndex('RightShank',HM.Segments);
            LShankIndex = getVecIndex('LeftShank',HM.Segments);
            RFootIndex = getVecIndex('RightFoot',HM.Segments);            
            LFootIndex = getVecIndex('LeftFoot',HM.Segments);
            
            if PelvisIndex ~= 0
                RHJCIndex = getLocVecIndex('RHJC',HM.Segments(PelvisIndex).LocalPoints);
                LHJCIndex = getLocVecIndex('LHJC',HM.Segments(PelvisIndex).LocalPoints);
                PointTranslator(1,:) = {'MLJC','lumbar-sacrum-joint','GLK'};
                PointTranslator(2,:) = {'LHJC','hip-joint-l','GHUL'};
                PointTranslator(3,:) = {'RHJC','hip-joint-r','GHUR'};
                TranslatorIndex = 3+1;
                if RThighIndex ~= 0
                    RKJCIndex = getLocVecIndex('RKJC',HM.Segments(RThighIndex).LocalPoints);
                    PointTranslator(TranslatorIndex,:) = {'RKJC','knee-joint-r','GKNR'};
                    TranslatorIndex = TranslatorIndex+1;
                    if RShankIndex ~= 0
                        RAJCIndex = getLocVecIndex('RAJC',HM.Segments(RShankIndex).LocalPoints);
                        PointTranslator(TranslatorIndex,:) = {'RAJC','ankle-joint-r','GSPR'};
                        if RFootIndex ~= 0
                            RFM2Index = getLocVecIndex('RFM2',HM.Segments(RFootIndex).LocalALs);
                        end
                    end
                end
                if LThighIndex ~= 0
                    LKJCIndex = getLocVecIndex('LKJC',HM.Segments(LThighIndex).LocalPoints);
                    PointTranslator(TranslatorIndex,:) = {'LKJC','knee-joint-l','GKNL'};
                    TranslatorIndex = TranslatorIndex+1;
                    if LShankIndex ~= 0
                        LAJCIndex = getLocVecIndex('LAJC',HM.Segments(LShankIndex).LocalPoints);
                        PointTranslator(TranslatorIndex,:) = {'LAJC','ankle-joint-l','GSPL'};
                        if LFootIndex
                            LFM2Index = getLocVecIndex('LFM2',HM.Segments(LFootIndex).LocalALs);
                        end
                    end
                end    
            end
%             TranslatorIndex = 1;
%             if PelvisIndex ~= 0
%                 RHJCIndex = getLocVecIndex('RHJC',HM.Segments(PelvisIndex).LocalPoints);
%                 LHJCIndex = getLocVecIndex('LHJC',HM.Segments(PelvisIndex).LocalPoints);
%                 PointTranslator(1,:) = {'MLJC','lumbar-sacrum-joint','GLK'};
%                 PointTranslator(2,:) = {'LHJC','hip-joint-l','GHUL'};
%                 PointTranslator(3,:) = {'RHJC','hip-joint-r','GHUR'};
%                 TranslatorIndex = 3+1;
%             end
%             if LThighIndex ~= 0
%                 LKJCIndex = getLocVecIndex('LKJC',HM.Segments(LThighIndex).LocalPoints);
%                 PointTranslator(TranslatorIndex,:) = {'LKJC','knee-joint-l','GKNL'};
%                 TranslatorIndex = TranslatorIndex+1;
%             end
%             if LShankIndex ~= 0
%                 LAJCIndex = getLocVecIndex('LAJC',HM.Segments(LShankIndex).LocalPoints);
%                 PointTranslator(TranslatorIndex,:) = {'LAJC','ankle-joint-l','GSPL'};
%             end
%             if LFootIndex
%                 LFM2Index = getLocVecIndex('LFM2',HM.Segments(LFootIndex).LocalALs);
%             end
            Angle.GHZ_a1  = 0*pi/180; Angle.GHZ_a2  = 0*pi/180; Angle.GHZ_a3  = 0*pi/180;
            Angle.GHUL_a3 = 0*pi/180; Angle.GHUL_a2 = 0*pi/180; Angle.GHUL_a1 = 0*pi/180;
            Angle.GKNL_a2 = 0*pi/180; Angle.GKNL_a1 = 0*pi/180; Angle.GKNL_a3 = 0*pi/180;
            Angle.GSPL_a3 = 0*pi/180; Angle.GSPL_a1 = 0*pi/180; Angle.GSPL_a2 = 0*pi/180;
            Angle.GHUR_a3 = 0*pi/180; Angle.GHUR_a2 = 0*pi/180; Angle.GHUR_a1 = 0*pi/180;
            Angle.GKNR_a2 = 0*pi/180; Angle.GKNR_a1 = 0*pi/180; Angle.GKNR_a3 = 0*pi/180;
            Angle.GSPR_a3 = 0*pi/180; Angle.GSPR_a1 = 0*pi/180; Angle.GSPR_a2 = 0*pi/180;
            d0 = zeros(3,1);
            
%             Angle.GHZ_a1  = -92.7*pi/180; Angle.GHZ_a2  = -2.138*pi/180; Angle.GHZ_a3  = -1.554*pi/180;
%             Angle.GHUL_a3 = -2.61*pi/180; Angle.GHUL_a2 = 0*pi/180; Angle.GHUL_a1 = 0*pi/180;
%             Angle.GKNL_a2 = -12*pi/180; Angle.GKNL_a1 = 13*pi/180; Angle.GKNL_a3 = 0*pi/180;
%             Angle.GSPL_a3 = 0*pi/180; Angle.GSPL_a1 = 0*pi/180; Angle.GSPL_a2 = 0*pi/180;
            
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
            % ILARIA FOR 10_VH_BMW2
%             Angle.GHZ_a1  = 88.1*pi/180; Angle.GHZ_a2  = -0.2*pi/180; Angle.GHZ_a3  = -38.0*pi/180;
%             Angle.GHUL_a3 = 75.3*pi/180; Angle.GHUL_a2 = -8.1*pi/180; Angle.GHUL_a1 = 9.7*pi/180;
%             Angle.GKNL_a2 = 78.4*pi/180; Angle.GKNL_a1 = 6.0*pi/180; Angle.GKNL_a3 = 0*pi/180;
%             Angle.GSPL_a3 = 112.7*pi/180; Angle.GSPL_a1 = -10.0*pi/180; Angle.GSPL_a2 = 0*pi/180;
%             d0 = zeros(3,1);
% -------------------------------------------------------------------------
%             % ILARIA FOR 05_SL_BMW1
%             Angle.GHZ_a1  = 87.3*pi/180; Angle.GHZ_a2  = 4.9*pi/180; Angle.GHZ_a3  = -19.0*pi/180;
%             Angle.GHUL_a3 = 87.7*pi/180; Angle.GHUL_a2 = -5.2*pi/180; Angle.GHUL_a1 = 11.4*pi/180;
%             Angle.GKNL_a2 = 71.5*pi/180; Angle.GKNL_a1 = 3.5*pi/180; Angle.GKNL_a3 = 0*pi/180;
%             Angle.GSPL_a3 = 109.2*pi/180; Angle.GSPL_a1 = -9.7*pi/180; Angle.GSPL_a2 = 0*pi/180;
%             d0 = zeros(3,1);
% -------------------------------------------------------------------------
%             % ILARIA FOR 05_SL_BMW1
%             Angle.GHZ_a1  = 90.2*pi/180; Angle.GHZ_a2  = -2.6*pi/180; Angle.GHZ_a3  = -44.7*pi/180;
%             Angle.GHUL_a3 = 85.6*pi/180; Angle.GHUL_a2 = 8.5*pi/180; Angle.GHUL_a1 = -6.2*pi/180;
%             Angle.GKNL_a2 = 94.8*pi/180; Angle.GKNL_a1 = 15*pi/180; Angle.GKNL_a3 = 0*pi/180;
%             Angle.GSPL_a3 = 86.7*pi/180; Angle.GSPL_a1 = -21.9*pi/180; Angle.GSPL_a2 = 0*pi/180;
%             d0 = zeros(3,1);
% -------------------------------------------------------------------------
%             % ILARIA FOR 05_SL_BMW1
%             Angle.GHZ_a1  = 89.6*pi/180; Angle.GHZ_a2  = 1.2*pi/180; Angle.GHZ_a3  = -46.1*pi/180;
%             Angle.GHUL_a3 = 69*pi/180; Angle.GHUL_a2 = -2.8*pi/180; Angle.GHUL_a1 = -2.2*pi/180;
%             Angle.GKNL_a2 = 79*pi/180; Angle.GKNL_a1 = -1.4*pi/180; Angle.GKNL_a3 = 0*pi/180;
%             Angle.GSPL_a3 = 115.4*pi/180; Angle.GSPL_a1 = -13.7*pi/180; Angle.GSPL_a2 = 0*pi/180;
%             d0 = zeros(3,1);
% -------------------------------------------------------------------------
%             % ILARIA FOR 02_LZ_PCA2
%             Angle.GHZ_a1  = 91.2*pi/180; Angle.GHZ_a2  = 1.5*pi/180; Angle.GHZ_a3  = -39.6*pi/180;
%             Angle.GHUL_a3 = 74.82*pi/180; Angle.GHUL_a2 = 6.1*pi/180; Angle.GHUL_a1 = 1.14*pi/180;
%             Angle.GKNL_a2 = 82.4*pi/180; Angle.GKNL_a1 = 3.4*pi/180; Angle.GKNL_a3 = 0*pi/180;
%             Angle.GSPL_a3 = 96.1*pi/180; Angle.GSPL_a1 = -12.5*pi/180; Angle.GSPL_a2 = 0*pi/180;
%             d0 = zeros(3,1);
% -------------------------------------------------------------------------
%             % ILARIA FOR 08_ZS_PCA1
%             Angle.GHZ_a1  = 91.05*pi/180; Angle.GHZ_a2  = -2.9*pi/180; Angle.GHZ_a3  = -48.8*pi/180;
%             Angle.GHUL_a3 = 54.5*pi/180; Angle.GHUL_a2 = 1.89*pi/180; Angle.GHUL_a1 = 12.14*pi/180;
%             Angle.GKNL_a2 = 75.65*pi/180; Angle.GKNL_a1 = -2.25*pi/180; Angle.GKNL_a3 = 0*pi/180;
%             Angle.GSPL_a3 = 85.62*pi/180; Angle.GSPL_a1 = -22.15*pi/180; Angle.GSPL_a2 = 0*pi/180;
%             d0 = zeros(3,1);
% -------------------------------------------------------------------------
%             % ILARIA FOR 03_FD_BMW1
%             Angle.GHZ_a1  = 88.65*pi/180; Angle.GHZ_a2  = 1.4*pi/180; Angle.GHZ_a3  = -46.15*pi/180;
%             Angle.GHUL_a3 = 72.2*pi/180; Angle.GHUL_a2 = 4.9*pi/180; Angle.GHUL_a1 = -9.25*pi/180;
%             Angle.GKNL_a2 = 77.8*pi/180; Angle.GKNL_a1 = 5.85*pi/180; Angle.GKNL_a3 = 0*pi/180;
%             Angle.GSPL_a3 = 95.5*pi/180; Angle.GSPL_a1 = 3.7*pi/180; Angle.GSPL_a2 = 0*pi/180;
%             d0 = zeros(3,1);
% -------------------------------------------------------------------------
%           ILARIA FOR 09_HB_BMW2
%             Angle.GHZ_a1  = 90.47*pi/180; Angle.GHZ_a2  = -7.7*pi/180; Angle.GHZ_a3  = -42.65*pi/180;
%             Angle.GHUL_a3 = 81.9*pi/180; Angle.GHUL_a2 = 6.8*pi/180; Angle.GHUL_a1 = -0.3*pi/180;
%             Angle.GKNL_a2 = 93.5*pi/180; Angle.GKNL_a1 = 26.3*pi/180; Angle.GKNL_a3 = 0*pi/180;
%             Angle.GSPL_a3 = 104.5*pi/180; Angle.GSPL_a1 = -27.8*pi/180; Angle.GSPL_a2 = 0*pi/180;
%             d0 = zeros(3,1);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
            
            GHZ_R_GHZ = rot123s(Angle.GHZ_a1,Angle.GHZ_a2,Angle.GHZ_a3);
            GHUL_R_GHUL = rot321s(Angle.GHUL_a3,Angle.GHUL_a2,Angle.GHUL_a1);
            GKNL_R_GKNL = rot21(Angle.GKNL_a2,Angle.GKNL_a1);
            GSPL_R_GSPL = rot31(Angle.GSPL_a3,Angle.GSPL_a1);
            GHUR_R_GHUR = rot321s(Angle.GHUL_a3,Angle.GHUL_a2,Angle.GHUL_a1);
            GKNR_R_GKNR = rot21(Angle.GKNL_a2,Angle.GKNL_a1);
            GSPR_R_GSPR = rot31(Angle.GSPL_a3,Angle.GSPL_a1);
            
            Glob_R_GHZ = rot_y(-pi/2)*rot_x(pi/2);
            GHZ_R_BEC = eye(3);
            BEC_R_GHUL = rot_y(pi);
            GHUL_R_OSL = eye(3);
            OSL_R_GKNL = rot_x(-pi/2);
            GKNL_R_USL = rot_x(pi/2);
            USL_R_GSPL = eye(3);
            GSPL_R_FUL = eye;
            FUL_R_GFBL= eye(3);
            BEC_R_GHUR = rot_y(pi);
            GHUR_R_OSR = eye(3);
            OSR_R_GKNR = rot_x(-pi/2);
            GKNR_R_USR = rot_x(pi/2);
            USR_R_GSPR = eye(3);
            GSPR_R_FUR = eye(3);
            FUR_R_GFBR= eye(3);
            BEC_R_GLK = eye(3);
            
            
            Glob_R_BEC = Glob_R_GHZ*GHZ_R_GHZ*GHZ_R_BEC;
            Glob_T_BEC = [Glob_R_BEC,d0;0,0,0,1];
            Glob_R_OSL = Glob_R_BEC*BEC_R_GHUL*GHUL_R_GHUL*GHUL_R_OSL;
            Glob_R_USL = Glob_R_OSL*OSL_R_GKNL*GKNL_R_GKNL*GKNL_R_USL;
            Glob_R_FUL = Glob_R_USL*USL_R_GSPL*GSPL_R_GSPL*GSPL_R_FUL;
            Glob_R_OSR = Glob_R_BEC*BEC_R_GHUR*GHUR_R_GHUR*GHUR_R_OSR;
            Glob_R_USR = Glob_R_OSR*OSR_R_GKNR*GKNR_R_GKNR*GKNR_R_USR;
            Glob_R_FUR = Glob_R_USR*USR_R_GSPR*GSPR_R_GSPR*GSPR_R_FUR;
            
            if PelvisIndex ~= 0
                Pelvis_Pos_GHZ = 0.5*(HM.Segments(PelvisIndex).LocalPoints(RHJCIndex).LocCoord + ...
                                    HM.Segments(PelvisIndex).LocalPoints(LHJCIndex).LocCoord);
                Glob_Pos_GHZ = Glob_T_BEC*HM.Segments(PelvisIndex).RamsisLCS_T_ALLCS*[Pelvis_Pos_GHZ;1];
            end
            % open file
            fid = fopen([Path,FileNameExt],'w');
            % Header
            fprintf(fid,'RAMSIS_MODEL_DATA_EXPORT_FILE V2.00\n');
            % Points
            fprintf(fid,'SECTION SKELETON_POINTS\n');
            fprintf(fid,['     GHZ	          hip-center\t',num2str(Glob_Pos_GHZ(1:3)','%-11.3f'),'\n']);
            NLPoints = size(HM.Segments(PelvisIndex).LocalPoints);
            for i=1:NLPoints
                PointName = HM.Segments(PelvisIndex).LocalPoints(i).Point.Name;
                PLocCoord = HM.Segments(PelvisIndex).LocalPoints(i).LocCoord;
                Glob_Pos_Point = Glob_T_BEC*HM.Segments(PelvisIndex).RamsisLCS_T_ALLCS*[PLocCoord;1];
                if strcmpi(PointName,'LHJC')
                    Glob_Pos_GHUL = Glob_Pos_Point;
                elseif strcmpi(PointName,'RHJC')
                    Glob_Pos_GHUR = Glob_Pos_Point;
                elseif strcmpi(PointName,'MLJC')
                    Glob_Pos_GLK = Glob_Pos_Point;                    
                end
                PointIndex = getCellIndex(PointTranslator,PointName);
                fprintf(fid,['\t',PointTranslator{PointIndex,3},'\t\t\t',PointTranslator{PointIndex,2},'\t\t',num2str(Glob_Pos_Point(1:3)','%-11.3f'),'\n']);
            end
            if LThighIndex ~= 0
                Glob_T_OSL = [Glob_R_OSL,Glob_Pos_GHUL(1:3);0,0,0,1];
                LCS_Pos_GKNL = HM.Segments(LThighIndex).LocalPoints(LKJCIndex).LocCoord;
                Glob_Pos_GKNL = Glob_T_OSL*HM.Segments(LThighIndex).RamsisLCS_T_ALLCS*[LCS_Pos_GKNL;1];
                PointIndex = getCellIndex(PointTranslator,'LKJC');
                fprintf(fid,['\t',PointTranslator{PointIndex,3},'\t\t\t',PointTranslator{PointIndex,2},'\t\t',num2str(Glob_Pos_GKNL(1:3)','%-11.3f'),'\n']);
                if LShankIndex ~= 0
                    Glob_T_USL = [Glob_R_USL,Glob_Pos_GKNL(1:3);0,0,0,1];
                    LCS_Pos_GSPL = HM.Segments(LShankIndex).LocalPoints(LAJCIndex).LocCoord;
                    Glob_Pos_GSPL = Glob_T_USL*HM.Segments(LShankIndex).RamsisLCS_T_ALLCS*[LCS_Pos_GSPL;1];
                    PointIndex = getCellIndex(PointTranslator,'LAJC');
                    fprintf(fid,['\t',PointTranslator{PointIndex,3},'\t\t\t',PointTranslator{PointIndex,2},'\t\t',num2str(Glob_Pos_GSPL(1:3)','%-11.3f'),'\n']);
                    if LFootIndex ~= 0
                        Glob_T_FUL = [Glob_R_FUL,Glob_Pos_GSPL(1:3);0,0,0,1];
                        %             LCS_Pos_GFBL = HM.Segments(FootIndex).LocalALs(LFM2Index).LocCoord;
                        LCS_Pos_GFBL = [0.13594; -0.06066; -0.00002]*1000;
                        Glob_Pos_GFBL = Glob_T_FUL*HM.Segments(LFootIndex).RamsisLCS_T_ALLCS*[LCS_Pos_GFBL;1];
                        fprintf(fid,['\tGFBL \t\t\tball-joint-l\t\t',num2str(Glob_Pos_GFBL(1:3)','%-11.3f'),'\n']);
                    end
                end
            end
            if RThighIndex ~= 0
                Glob_T_OSR = [Glob_R_OSR,Glob_Pos_GHUR(1:3);0,0,0,1];
                LCS_Pos_GKNR = HM.Segments(RThighIndex).LocalPoints(RKJCIndex).LocCoord;
                Glob_Pos_GKNR = Glob_T_OSR*HM.Segments(RThighIndex).RamsisLCS_T_ALLCS*[LCS_Pos_GKNR;1];
                PointIndex = getCellIndex(PointTranslator,'RKJC');
                fprintf(fid,['\t',PointTranslator{PointIndex,3},'\t\t\t',PointTranslator{PointIndex,2},'\t\t',num2str(Glob_Pos_GKNR(1:3)','%-11.3f'),'\n']);
                if RShankIndex ~= 0
                    Glob_T_USR = [Glob_R_USR,Glob_Pos_GKNR(1:3);0,0,0,1];
                    LCS_Pos_GSPR = HM.Segments(RShankIndex).LocalPoints(RAJCIndex).LocCoord;
                    Glob_Pos_GSPR = Glob_T_USR*HM.Segments(RShankIndex).RamsisLCS_T_ALLCS*[LCS_Pos_GSPR;1];
                    PointIndex = getCellIndex(PointTranslator,'RAJC');
                    fprintf(fid,['\t',PointTranslator{PointIndex,3},'\t\t\t',PointTranslator{PointIndex,2},'\t\t',num2str(Glob_Pos_GSPR(1:3)','%-11.3f'),'\n']);
                    if RFootIndex ~= 0
                        Glob_T_FUR = [Glob_R_FUR,Glob_Pos_GSPR(1:3);0,0,0,1];
                        %             LCS_Pos_GFBR = HM.Segments(FootIndex).LocalALs(RFM2Index).LocCoord;
                        LCS_Pos_GFBR = [0.13594; -0.06066; -0.00002]*1000;
                        Glob_Pos_GFBR = Glob_T_FUR*HM.Segments(RFootIndex).RamsisLCS_T_ALLCS*[LCS_Pos_GFBR;1];
                        fprintf(fid,['\tGFBR \t\t\tball-joint-l\t\t',num2str(Glob_Pos_GFBR(1:3)','%-11.3f'),'\n']);
                    end
                end
            end
            fprintf(fid,'ENDSECTION SKELETON_POINTS\n');
%             % Joints
            fprintf(fid,'SECTION JOINTS\n');
            fprintf(fid,['\tGHZ \t',num2str(Glob_Pos_GHZ(1:3)','%-11.3f'),'\t',num2str(Angle.GHZ_a1*180/pi,'%-8.3f'),'\t',...
                num2str(Angle.GHZ_a2*180/pi,'%-8.3f'),'\t',num2str(Angle.GHZ_a3*180/pi,'%-8.3f'),'\n']);
            if LThighIndex ~= 0
                fprintf(fid,['\tGHUL \t',num2str(Glob_Pos_GHUL(1:3)','%-11.3f'),'\t',num2str((Angle.GHUL_a1*180/pi),'%-8.3f'),'\t',...
                    num2str((Angle.GHUL_a2*180/pi),'%-8.3f'),'\t',num2str((Angle.GHUL_a3*180/pi),'%-8.3f'),'\n']);
                if LShankIndex ~= 0
                    fprintf(fid,['\tGKNL \t',num2str(Glob_Pos_GKNL(1:3)','%-11.3f'),'\t',num2str(Angle.GKNL_a1*180/pi,'%-8.3f'),'\t',...
                        num2str(Angle.GKNL_a2*180/pi,'%-8.3f'),'\t',num2str(Angle.GKNL_a3*180/pi,'%-8.3f'),'\n']);
                    if LFootIndex ~= 0
                        fprintf(fid,['\tGSPL \t',num2str(Glob_Pos_GSPL(1:3)','%-11.3f'),'\t',num2str(Angle.GSPL_a1*180/pi,'%-8.3f'),'\t',...
                            num2str(Angle.GSPL_a2*180/pi,'%-8.3f'),'\t',num2str(Angle.GSPL_a3*180/pi,'%-8.3f'),'\n']);
                    end
                end
            end
            if RThighIndex ~= 0
                fprintf(fid,['\tGHUR \t',num2str(Glob_Pos_GHUR(1:3)','%-11.3f'),'\t',num2str((Angle.GHUR_a1*180/pi),'%-8.3f'),'\t',...
                    num2str((Angle.GHUR_a2*180/pi),'%-8.3f'),'\t',num2str((Angle.GHUR_a3*180/pi),'%-8.3f'),'\n']);
                if RShankIndex ~= 0
                    fprintf(fid,['\tGKNR \t',num2str(Glob_Pos_GKNR(1:3)','%-11.3f'),'\t',num2str(Angle.GKNR_a1*180/pi,'%-8.3f'),'\t',...
                        num2str(Angle.GKNR_a2*180/pi,'%-8.3f'),'\t',num2str(Angle.GKNR_a3*180/pi,'%-8.3f'),'\n']);
                    if RFootIndex ~= 0
                        fprintf(fid,['\tGSPR \t',num2str(Glob_Pos_GSPR(1:3)','%-11.3f'),'\t',num2str(Angle.GSPR_a1*180/pi,'%-8.3f'),'\t',...
                            num2str(Angle.GSPR_a2*180/pi,'%-8.3f'),'\t',num2str(Angle.GSPR_a3*180/pi,'%-8.3f'),'\n']);
                    end
                end
            end
            fprintf(fid,'ENDSECTION JOINTS\n');
%             % Joint System
            fprintf(fid,'SECTION JOINTS_COORDINATE_SYSTEMS\n');
            fprintf(fid,['\tGHZ \t',num2str(Glob_Pos_GHZ(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GHZ(:,1)','%-8.3f'),'\t',...
                num2str(Glob_R_GHZ(:,2)','%-8.3f'),'\t',num2str(Glob_R_GHZ(:,3)','%-8.3f'),'\n']);
            Glob_R_GLK = Glob_R_BEC*BEC_R_GLK;
            fprintf(fid,['\tGLK \t',num2str(Glob_Pos_GLK(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GLK(:,1)','%-8.3f'),'\t',...
                    num2str(Glob_R_GLK(:,2)','%-8.3f'),'\t',num2str(Glob_R_GLK(:,3)','%-8.3f'),'\n']);
            if LThighIndex ~= 0
                Glob_R_GHUL = Glob_R_BEC*BEC_R_GHUL;
                fprintf(fid,['\tGHUL \t',num2str(Glob_Pos_GHUL(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GHUL(:,1)','%-8.3f'),'\t',...
                    num2str(Glob_R_GHUL(:,2)','%-8.3f'),'\t',num2str(Glob_R_GHUL(:,3)','%-8.3f'),'\n']);
                if LShankIndex ~= 0
                    Glob_R_GKNL = Glob_R_OSL*OSL_R_GKNL;
                    fprintf(fid,['\tGKNL \t',num2str(Glob_Pos_GKNL(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GKNL(:,1)','%-8.3f'),'\t',...
                        num2str(Glob_R_GKNL(:,2)','%-8.3f'),'\t',num2str(Glob_R_GKNL(:,3)','%-8.3f'),'\n']);
                    if LFootIndex ~= 0
                        Glob_R_GSPL = Glob_R_USL*USL_R_GSPL;
                        fprintf(fid,['\tGSPL \t',num2str(Glob_Pos_GSPL(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GSPL(:,1)','%-8.3f'),'\t',...
                            num2str(Glob_R_GSPL(:,2)','%-8.3f'),'\t',num2str(Glob_R_GSPL(:,3)','%-8.3f'),'\n']);
                        Glob_R_GFBL = Glob_R_FUL*FUL_R_GFBL;
                        fprintf(fid,['\tGFBL \t',num2str(Glob_Pos_GFBL(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GFBL(:,1)','%-8.3f'),'\t',...
                            num2str(Glob_R_GFBL(:,2)','%-8.3f'),'\t',num2str(Glob_R_GFBL(:,3)','%-8.3f'),'\n']);
                    end
                end
            end
            if RThighIndex ~= 0
                Glob_R_GHUR = Glob_R_BEC*BEC_R_GHUR;
                fprintf(fid,['\tGHUR \t',num2str(Glob_Pos_GHUR(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GHUR(:,1)','%-8.3f'),'\t',...
                    num2str(Glob_R_GHUR(:,2)','%-8.3f'),'\t',num2str(Glob_R_GHUR(:,3)','%-8.3f'),'\n']);
                if RShankIndex ~= 0
                    Glob_R_GKNR = Glob_R_OSR*OSR_R_GKNR;
                    fprintf(fid,['\tGKNR \t',num2str(Glob_Pos_GKNR(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GKNR(:,1)','%-8.3f'),'\t',...
                        num2str(Glob_R_GKNR(:,2)','%-8.3f'),'\t',num2str(Glob_R_GKNR(:,3)','%-8.3f'),'\n']);
                    if RFootIndex ~= 0
                        Glob_R_GSPR = Glob_R_USR*USR_R_GSPR;
                        fprintf(fid,['\tGSPR \t',num2str(Glob_Pos_GSPR(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GSPR(:,1)','%-8.3f'),'\t',...
                            num2str(Glob_R_GSPR(:,2)','%-8.3f'),'\t',num2str(Glob_R_GSPR(:,3)','%-8.3f'),'\n']);
                        Glob_R_GFBR = Glob_R_FUR*FUR_R_GFBR;
                        fprintf(fid,['\tGFBR \t',num2str(Glob_Pos_GFBR(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GFBR(:,1)','%-8.3f'),'\t',...
                            num2str(Glob_R_GFBR(:,2)','%-8.3f'),'\t',num2str(Glob_R_GFBR(:,3)','%-8.3f'),'\n']);
                    end
                end
            end
            fprintf(fid,'ENDSECTION JOINTS_COORDINATE_SYSTEMS\n');
%             % Body System
            fprintf(fid,'SECTION BODY_ELEMENT_COORDINATE_SYSTEMS\n');
            fprintf(fid,['\tGHZ \t',num2str(Glob_Pos_GHZ(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_BEC(:,1)','%-8.3f'),'\t',...
                num2str(Glob_R_BEC(:,2)','%-8.3f'),'\t',num2str(Glob_R_BEC(:,3)','%-8.3f'),'\n']);
            if LThighIndex ~= 0
                fprintf(fid,['\tGHUL \t',num2str(Glob_Pos_GHUL(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_OSL(:,1)','%-8.3f'),'\t',...
                    num2str(Glob_R_OSL(:,2)','%-8.3f'),'\t',num2str(Glob_R_OSL(:,3)','%-8.3f'),'\n']);
                if LShankIndex ~= 0
                    fprintf(fid,['\tGKNL \t',num2str(Glob_Pos_GKNL(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_USL(:,1)','%-8.3f'),'\t',...
                        num2str(Glob_R_USL(:,2)','%-8.3f'),'\t',num2str(Glob_R_USL(:,3)','%-8.3f'),'\n']);
                    if LFootIndex ~= 0
                        fprintf(fid,['\tGSPL \t',num2str(Glob_Pos_GSPL(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_FUL(:,1)','%-8.3f'),'\t',...
                            num2str(Glob_R_FUL(:,2)','%-8.3f'),'\t',num2str(Glob_R_FUL(:,3)','%-8.3f'),'\n']);
                    end
                end
            end
            if RThighIndex ~= 0
                fprintf(fid,['\tGHUR \t',num2str(Glob_Pos_GHUR(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_OSR(:,1)','%-8.3f'),'\t',...
                    num2str(Glob_R_OSR(:,2)','%-8.3f'),'\t',num2str(Glob_R_OSR(:,3)','%-8.3f'),'\n']);
                if RShankIndex ~= 0
                    fprintf(fid,['\tGKNR \t',num2str(Glob_Pos_GKNR(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_USR(:,1)','%-8.3f'),'\t',...
                        num2str(Glob_R_USR(:,2)','%-8.3f'),'\t',num2str(Glob_R_USR(:,3)','%-8.3f'),'\n']);
                    if RFootIndex ~= 0
                        fprintf(fid,['\tGSPR \t',num2str(Glob_Pos_GSPR(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_FUR(:,1)','%-8.3f'),'\t',...
                            num2str(Glob_R_FUR(:,2)','%-8.3f'),'\t',num2str(Glob_R_FUR(:,3)','%-8.3f'),'\n']);
                    end
                end
            end
            fprintf(fid,'ENDSECTION BODY_ELEMENT_COORDINATE_SYSTEMS\n'); 
            % Skin_Points
            fprintf(fid,'SECTION SKIN_POINTS\n');
            fprintf(fid,'   BEC_1_1	  -931.962    764.185    779.167\n');
            fprintf(fid,'ENDSECTION SKIN_POINTS\n');
%             % Marker
            fprintf(fid,'SECTION MARKER\n');
            SegmentTranslator(1,:) = {'Pelvis','BEC',Glob_T_BEC};
            SegmentTranslatorIndex = 2;
            if LThighIndex ~= 0
                SegmentTranslator(SegmentTranslatorIndex,:) = {'LeftThigh','OSL',Glob_T_OSL};
                SegmentTranslatorIndex = SegmentTranslatorIndex +1;
                if LShankIndex ~= 0
                    SegmentTranslator(SegmentTranslatorIndex,:) = {'LeftShank','USL',Glob_T_USL};
                    SegmentTranslatorIndex = SegmentTranslatorIndex +1;
                    if LFootIndex ~= 0
                        SegmentTranslator(SegmentTranslatorIndex,:) = {'LeftFoot','FUL',Glob_T_FUL};
                        SegmentTranslatorIndex = SegmentTranslatorIndex +1;
                    end
                end
            end
            if RThighIndex ~= 0
                SegmentTranslator(SegmentTranslatorIndex,:) = {'RightThigh','OSR',Glob_T_OSR};
                SegmentTranslatorIndex = SegmentTranslatorIndex +1;
                if RShankIndex ~= 0
                    SegmentTranslator(SegmentTranslatorIndex,:) = {'RightShank','USR',Glob_T_USR};
                    SegmentTranslatorIndex = SegmentTranslatorIndex +1;
                    if RFootIndex ~= 0
                        SegmentTranslator(SegmentTranslatorIndex,:) = {'RightFoot','FUR',Glob_T_FUR};
                        SegmentTranslatorIndex = SegmentTranslatorIndex +1;
                    end
                end
            end
            for i=1:NSegments
                NMarkers = size(HM.Segments(i).LocalMarkers,1);
                SegName = HM.Segments(i).Name;
                SegIndex = getCellIndex(SegmentTranslator,SegName);
                for j=1:NMarkers
                    MarkerName = HM.Segments(i).LocalMarkers(j).Point.Name;
                    MLocCoord = HM.Segments(i).LocalMarkers(j).LocCoord;
                    Glob_Pos_Marker = SegmentTranslator{SegIndex,3}*HM.Segments(i).RamsisLCS_T_ALLCS*[MLocCoord;1];
                    fprintf(fid,['\t',MarkerName,'\t\t',SegmentTranslator{SegIndex,2},'\t\t\t',num2str(Glob_Pos_Marker(1:3)','%-11.3f'),'\n']);
                end
            end
            
            fprintf(fid,'ENDSECTION MARKER\n');
%             % Centres of gravity
            fprintf(fid,'SECTION CENTRES_OF_GRAVITY\n');
            for i=1:NSegments
                SegName = HM.Segments(i).Name;
                SegIndex = getCellIndex(SegmentTranslator,SegName);
                Loc_CoM = HM.Segments(i).CoM.LocCoord;
                Glob_CoM = SegmentTranslator{SegIndex,3}*HM.Segments(i).RamsisLCS_T_ALLCS*[Loc_CoM;1];
                fprintf(fid,['\t',SegmentTranslator{SegIndex,2},'\t\t',num2str(Glob_CoM(1:3)','%-11.3f'),'\n']);
            end
            fprintf(fid,'ENDSECTION CENTRES_OF_GRAVITY\n');
%             % Mass
            fprintf(fid,'SECTION MASS\n');
            for i=1:NSegments
                SegName = HM.Segments(i).Name;
                SegIndex = getCellIndex(SegmentTranslator,SegName);
                fprintf(fid,['\t',SegmentTranslator{SegIndex,2},'\t\t',num2str(HM.Segments(i).Mass','%-11.3f'),'\n']);
            end
            fprintf(fid,'ENDSECTION MASS\n');
            fprintf(fid,'SECTION MOMENT_OF_INERTIA\n');
            for i=1:NSegments
                SegName = HM.Segments(i).Name;
                SegIndex = getCellIndex(SegmentTranslator,SegName);
                ALLCS_I_CM = HM.Segments(i).I;
                RamsisLCS_I_CM = HM.Segments(i).RamsisLCS_T_ALLCS(1:3,1:3) * ALLCS_I_CM * HM.Segments(i).RamsisLCS_T_ALLCS(1:3,1:3)';
                fprintf(fid,['\t',SegmentTranslator{SegIndex,2},'\t\t',...
                    num2str(RamsisLCS_I_CM(1,1),'%-11.3f'),'\t\t',num2str(RamsisLCS_I_CM(1,2),'%-11.3f'),'\t\t',num2str(RamsisLCS_I_CM(1,3),'%-11.3f'),'\t\t',...
                    num2str(RamsisLCS_I_CM(2,1),'%-11.3f'),'\t\t',num2str(RamsisLCS_I_CM(2,2),'%-11.3f'),'\t\t',num2str(RamsisLCS_I_CM(2,3),'%-11.3f'),'\t\t',...
                    num2str(RamsisLCS_I_CM(3,1),'%-11.3f'),'\t\t',num2str(RamsisLCS_I_CM(3,2),'%-11.3f'),'\t\t',num2str(RamsisLCS_I_CM(3,3),'%-11.3f'),'\t\t','\n']);
            end
            fprintf(fid,'ENDSECTION MOMENT_OF_INERTIA\n');
            % close file
            fclose(fid);
            
        end
        function writeSubjectParsIntermed(HM,FileName,Path)
%             FileNameExt = [FileName,'.xml'];
%             % open file
%             fid = fopen([Path,FileNameExt],'w');
%             % Header
%             fprintf(fid,'?<?xml version="1.0" encoding="utf-8"?>\n');
%             fprintf(fid,'<?SubjectParameterIntermediate?>\n');
%             fprintf(fid,'<Subject>\n');
%             fprintf(fid,'\t<Segment Name="Pelvis" Mass="33">\n');
%             fprintf(fid,'\t\t<AL Name="RIAS" XYZ_in_AL-LCS="[0;0;200]"/>\n');
%             fprintf(fid,'\t</Segment>\n');
%             fprintf(fid,'</Subject>\n');
            NSegments = size(HM.Segments,1);
            NJoints = size(HM.Joints,1);
            NSubjectPar = size(HM.SubjectParameters,2);
            % Create a  XML document.
            docNode = com.mathworks.xml.XMLUtils.createDocument('Subject');
            docRootNode = docNode.getDocumentElement;
            % Fill all segments elements
            for i=1:NSegments
                SegmentElement = HM.Segments(i).writeSubjectParsIntermed(docNode);
                docRootNode.appendChild(SegmentElement);
            end
            % Fill Joints elements
%             for i=1:NJoints
%                 JointElement = HM.Joints{i}.writeSubjectParsIntermed(docNode);
%                 docRootNode.appendChild(JointElement);
%             end
            % write subject parameters from alm file
            SubPartElement = docNode.createElement('Subject_Parameters');
            docRootNode.appendChild(SubPartElement);
            for i=1:NSubjectPar
                Space = findstr(HM.SubjectParameters(i).Keyword,' ');
                if ~isspace(Space)
                    NewAcronym = regexprep(HM.SubjectParameters(i).Acronym,' ','_');
                    SubElement = SubPartElement.appendChild(docNode.createElement(NewAcronym));                    
                else
                    SubElement = SubPartElement.appendChild(docNode.createElement(HM.SubjectParameters(i).Acronym));
                end
                NameAttribute = docNode.createAttribute('Name');
                ValueAttribute = docNode.createAttribute('Value');
                NameAttribute.setValue(HM.SubjectParameters(i).Keyword);
                ValueAttribute.setValue(HM.SubjectParameters(i).Value);
                SubElement.setAttributeNode(NameAttribute);
                SubElement.setAttributeNode(ValueAttribute);
            end
            % Save the sample XML document.
            xmlFileName = [Path,FileName,'.xml'];
            xmlwrite(xmlFileName,docNode);
%             edit(xmlFileName);
        end
        function writeEXP(HM,FilePath,FileName)
%             Tmp = load([FilePath, FileName]);
%             StructMat = Tmp.dim;
            deg2rad = pi/180;
            
%             Glob_Pos_GHZ = [300;-1080;460]; %mm S01
              Glob_Pos_GHZ = [420;-1028;335]; %mm S02
%             Angle.GHZ_a1  = -42*deg2rad; Angle.GHZ_a2  = 0*deg2rad; Angle.GHZ_a3  = -6*deg2rad;
%             Angle.GHUR_a3 = -14*deg2rad; Angle.GHUR_a2 = 0*deg2rad; Angle.GHUR_a1 = -8*deg2rad;
%             Angle.GKNR_a2 = 0*deg2rad; Angle.GKNR_a1 = 20*deg2rad; Angle.GKNR_a3 = 0*deg2rad;
%             Angle.GSPR_a3 = 90*deg2rad; Angle.GSPR_a1 = 4*deg2rad; Angle.GSPR_a2 = 0*deg2rad;
%             Angle.GFBR_a3 = 0*deg2rad; Angle.GFBR_a2 = 0*deg2rad; Angle.GFBR_a1 = 0*deg2rad;
%             Angle.GHUL_a3 = -15*deg2rad; Angle.GHUL_a2 = 5*deg2rad; Angle.GHUL_a1 = -25*deg2rad;
%             Angle.GKNL_a2 = -2*deg2rad; Angle.GKNL_a1 = 2*deg2rad; Angle.GKNL_a3 = 0*deg2rad;
%             Angle.GSPL_a3 = 90*deg2rad; Angle.GSPL_a1  = 0*deg2rad;Angle.GSPL_a2  = 0*deg2rad;
%             Angle.GFBL_a3 = 0*deg2rad; Angle.GFBL_a2 = 0*deg2rad; Angle.GFBL_a1 = 0*deg2rad;
%             Angle.GLK_a3 = 0*deg2rad; Angle.GLK_a2= 0*deg2rad;  Angle.GLK_a1 = 0*deg2rad;
%             Angle.GKH_a3 = -24*deg2rad; Angle.GKH_a2 = 0*deg2rad; Angle.GKH_a1= 0*deg2rad;
%             Angle.GBRK_a3 = -14*deg2rad; Angle.GBRK_a2 = 0*deg2rad; Angle.GBRK_a1 = 0*deg2rad;
%             Angle.GSBR_a3 = 15*deg2rad; Angle.GSBR_a2 = 20*deg2rad; Angle.GSBR_a1 = 0*deg2rad;
%             Angle.GSBL_a3 = -22*deg2rad; Angle.GSBL_a2 = 7*deg2rad; Angle.GSBL_a1 = 0*deg2rad;
%             Angle.GSR_a3  = -25*deg2rad;  Angle.GSR_a2 = 60*deg2rad; Angle.GSR_a1  = 0*deg2rad;
%             Angle.GSL_a3  = -20*deg2rad;  Angle.GSL_a2 = -60*deg2rad; Angle.GSL_a1  = 0*deg2rad;
%             Angle.GELR_a2 = -40*deg2rad; Angle.GELR_a1 = 5*deg2rad; Angle.GELR_a3 = 0*deg2rad;
%             Angle.GELL_a2 = -40*deg2rad; Angle.GELL_a1 = 10*deg2rad; Angle.GELL_a3 = 0*deg2rad;
%             Angle.GHAR_a3 = -15*deg2rad; Angle.GHAR_a2 = 0*deg2rad; Angle.GHAR_a1 = 0*deg2rad;
%             Angle.GHAL_a3 = -10*deg2rad; Angle.GHAL_a2 = 0*deg2rad; Angle.GHAL_a1 = 0*deg2rad;
            Angle.GHZ_a1  = 0*deg2rad; Angle.GHZ_a2  = 0*deg2rad; Angle.GHZ_a3  = 0*deg2rad;
            Angle.GHUR_a3 = 0*deg2rad; Angle.GHUR_a2 = 0*deg2rad; Angle.GHUR_a1 = 0*deg2rad;
            Angle.GKNR_a2 = 0*deg2rad; Angle.GKNR_a1 = 0*deg2rad; Angle.GKNR_a3 = 0*deg2rad;
            Angle.GSPR_a3 = 0*deg2rad; Angle.GSPR_a1 = 0*deg2rad; Angle.GSPR_a2 = 0*deg2rad;
            Angle.GFBR_a3 = 0*deg2rad; Angle.GFBR_a2 = 0*deg2rad; Angle.GFBR_a1 = 0*deg2rad;
            Angle.GHUL_a3 = 0*deg2rad; Angle.GHUL_a2 = 0*deg2rad; Angle.GHUL_a1 = 0*deg2rad;
            Angle.GKNL_a2 = 0*deg2rad; Angle.GKNL_a1 = 0*deg2rad; Angle.GKNL_a3 = 0*deg2rad;
            Angle.GSPL_a3 = 0*deg2rad; Angle.GSPL_a1  = 0*deg2rad;Angle.GSPL_a2  = 0*deg2rad;
            Angle.GFBL_a3 = 0*deg2rad; Angle.GFBL_a2 = 0*deg2rad; Angle.GFBL_a1 = 0*deg2rad;
            Angle.GLK_a3 = 0*deg2rad; Angle.GLK_a2= 0*deg2rad;  Angle.GLK_a1 = 0*deg2rad;
            Angle.GKH_a3 = 0*deg2rad; Angle.GKH_a2 = 0*deg2rad; Angle.GKH_a1= 0*deg2rad;
            Angle.GBRK_a3 = 0*deg2rad; Angle.GBRK_a2 = 0*deg2rad; Angle.GBRK_a1 = 0*deg2rad;
            Angle.GSBR_a3 = 0*deg2rad; Angle.GSBR_a2 = 0*deg2rad; Angle.GSBR_a1 = 0*deg2rad;
            Angle.GSBL_a3 = 0*deg2rad; Angle.GSBL_a2 = 0*deg2rad; Angle.GSBL_a1 = 0*deg2rad;
            Angle.GSR_a3  = 0*deg2rad;  Angle.GSR_a2 = 0*deg2rad; Angle.GSR_a1  = 0*deg2rad;
            Angle.GSL_a3  = 0*deg2rad;  Angle.GSL_a2 = 0*deg2rad; Angle.GSL_a1  = 0*deg2rad;
            Angle.GELR_a2 = 0*deg2rad; Angle.GELR_a1 = 0*deg2rad; Angle.GELR_a3 = 0*deg2rad;
            Angle.GELL_a2 = 0*deg2rad; Angle.GELL_a1 = 0*deg2rad; Angle.GELL_a3 = 0*deg2rad;
            Angle.GHAR_a3 = 0*deg2rad; Angle.GHAR_a2 = 0*deg2rad; Angle.GHAR_a1 = 0*deg2rad;
            Angle.GHAL_a3 = 0*deg2rad; Angle.GHAL_a2 = 0*deg2rad; Angle.GHAL_a1 = 0*deg2rad;
            
            Angle.PKSP_a1 =0*deg2rad; Angle.PKSP_a2 = 0*deg2rad; Angle.PKSP_a3 = 0*deg2rad;
            Angle.GM1L_a1 =0*deg2rad; Angle.GM1L_a2 = 0*deg2rad; Angle.GM1L_a3 = 0*deg2rad;
            Angle.GM1R_a1 =0*deg2rad; Angle.GM1R_a2 = 0*deg2rad; Angle.GM1R_a3 = 0*deg2rad;
            
            Glob_R_BEC = rot123s(Angle.GHZ_a1,Angle.GHZ_a2,Angle.GHZ_a3);
            BEC_R_OSR = rot321s(Angle.GHUR_a3,Angle.GHUR_a2,Angle.GHUR_a1);
            OSR_R_USR = rot21(Angle.GKNR_a2,Angle.GKNR_a1);
            USR_R_FUR = rot31(Angle.GSPR_a3,Angle.GSPR_a1);
%             FUR_R_FBR = rot_z(Angle.GFBR_a3);
            BEC_R_OSL = rot321s(Angle.GHUL_a3,Angle.GHUL_a2,Angle.GHUL_a1);
            OSL_R_USL = rot21(Angle.GKNL_a2,Angle.GKNL_a1);
            USL_R_FUL = rot31(Angle.GSPL_a3,Angle.GSPL_a1);
%             FUL_R_FBL = rot_z(Angle.GFBL_a3);
%             BEC_R_ULW = rot32(Angle.GLK_a3,Angle.GLK_a2);
%             BEC_R_ULW = rot321s(Angle.GLK_a3,Angle.GLK_a2,Angle.GLK_a1);
%             OLW_R_ULW = rot321s(Angle.GLL_a3,Angle.GLL_a2,Angle.GLL_a1);
%             ULW_R_UBW = rot321s(Angle.GBL_a3,Angle.GBL_a2,Angle.GBL_a1);
%             UBW_R_OBW = rot321s(Angle.GBB_a3,Angle.GBB_a2,Angle.GBB_a1);
%             OBW_R_UHW = rot321s(Angle.GHB_a3,Angle.GHB_a2,Angle.GHB_a1);
%             UHW_R_OHW = rot321s(Angle.GHH_a3,Angle.GHH_a2,Angle.GHH_a1);
            BRK_R_KO  = rot321s(Angle.GKH_a3,Angle.GKH_a2,Angle.GKH_a1);
            BEC_R_BRK = rot321s(Angle.GBRK_a3,Angle.GBRK_a2,Angle.GBRK_a1);
            BRK_R_SBR = rot321s(Angle.GSBR_a3,Angle.GSBR_a2,Angle.GSBR_a1);
            BRK_R_SBL = rot321s(Angle.GSBL_a3,Angle.GSBL_a2,Angle.GSBL_a1);
            SBR_R_OAR = rot321s(Angle.GSR_a3,Angle.GSR_a2,Angle.GSR_a1);
            SBL_R_OAL = rot321s(Angle.GSL_a3,Angle.GSL_a2,Angle.GSL_a1);
            OAR_R_UAR = rot21(Angle.GELR_a2,Angle.GELR_a1);
            OAL_R_UAL = rot21(Angle.GELL_a2,Angle.GELL_a1);
            UAR_R_HAR = rot32(Angle.GHAR_a3,Angle.GHAR_a2);
            UAL_R_HAL = rot32(Angle.GHAL_a3,Angle.GHAL_a2);
            
            Glob_R_BEC = rot_y(-pi/2)*rot_x(pi/2)*Glob_R_BEC;                   RigidTrans(1,:) = {'BEC',Glob_R_BEC};
            Glob_R_OSR = Glob_R_BEC * rot_y(pi) * BEC_R_OSR;                    RigidTrans = [RigidTrans; {'OSR',Glob_R_OSR}];
            Glob_R_OSL = Glob_R_BEC * rot_y(pi) * BEC_R_OSL;                    RigidTrans = [RigidTrans; {'OSL',Glob_R_OSL}];
            Glob_R_USR = Glob_R_OSR * rot_x(-pi/2) * OSR_R_USR * rot_x(pi/2);   RigidTrans = [RigidTrans; {'USR',Glob_R_USR}];
            Glob_R_USL = Glob_R_OSL * rot_x(-pi/2) * OSL_R_USL * rot_x(pi/2);   RigidTrans = [RigidTrans; {'USL',Glob_R_USL}];
            Glob_R_FUR = Glob_R_USR * USR_R_FUR * rot_z (pi/2);                 RigidTrans = [RigidTrans; {'FUR',Glob_R_FUR}];
            Glob_R_FUL = Glob_R_USL * USL_R_FUL * rot_z (pi/2);                 RigidTrans = [RigidTrans; {'FUL',Glob_R_FUL}];
%             Glob_R_FBR = Glob_R_FUR * FUR_R_FBR;                                RigidTrans = [RigidTrans; {'FBR',Glob_R_FBR}];
%             Glob_R_FBL = Glob_R_FUL * FUL_R_FBL;                                RigidTrans = [RigidTrans; {'FBL',Glob_R_FBL}];
%             Glob_R_ULW = Glob_R_BEC * BEC_R_ULW;                                RigidTrans = [RigidTrans; {'ULW',Glob_R_ULW}];
%             Glob_R_OLW = Glob_R_ULW * OLW_R_ULW;                                RigidTrans = [RigidTrans; {'OLW',Glob_R_OLW}];
%             Glob_R_UBW = Glob_R_OLW * ULW_R_UBW;                                RigidTrans = [RigidTrans; {'UBW',Glob_R_UBW}];
%             Glob_R_OBW = Glob_R_UBW * UBW_R_OBW;                                RigidTrans = [RigidTrans; {'OBW',Glob_R_OBW}];
%             Glob_R_UHW = Glob_R_OBW * OBW_R_UHW;                                RigidTrans = [RigidTrans; {'UHW',Glob_R_UHW}];
%             Glob_R_OHW = Glob_R_UHW * UHW_R_OHW;                                RigidTrans = [RigidTrans; {'OHW',Glob_R_OHW}];
            Glob_R_BRK = Glob_R_BEC *  rot_z(-pi/2) * rot_y(pi) * BEC_R_BRK;    RigidTrans = [RigidTrans; {'BRK',Glob_R_BRK}];
            Glob_R_KO  = Glob_R_BRK * rot_z(pi/2) * rot_x(pi)* BRK_R_KO;        RigidTrans = [RigidTrans; {'KO',Glob_R_KO}];
            Glob_R_SBR = Glob_R_BRK * rot_y(-pi/2) * rot_x(-pi/2) * BRK_R_SBR;  RigidTrans = [RigidTrans; {'SBR',Glob_R_SBR}];
            Glob_R_SBL = Glob_R_BRK * rot_y(pi/2) * rot_x(pi/2) *   BRK_R_SBL;  RigidTrans = [RigidTrans; {'SBL',Glob_R_SBL}];
            Glob_R_OAR = Glob_R_SBR * SBR_R_OAR;                                RigidTrans = [RigidTrans; {'OAR',Glob_R_OAR}];
            Glob_R_OAL = Glob_R_SBL * SBL_R_OAL;                                RigidTrans = [RigidTrans; {'OAL',Glob_R_OAL}];
            Glob_R_UAR = Glob_R_OAR * rot_x(-pi/2) * OAR_R_UAR * rot_x(pi/2);   RigidTrans = [RigidTrans; {'UAR',Glob_R_UAR}];
            Glob_R_UAL = Glob_R_OAL * rot_x(-pi/2) * OAL_R_UAL * rot_x(pi/2);   RigidTrans = [RigidTrans; {'UAL',Glob_R_UAL}];
            Glob_R_HAR = Glob_R_UAR * UAR_R_HAR;                                RigidTrans = [RigidTrans; {'HAR',Glob_R_HAR}];
            Glob_R_HAL = Glob_R_UAL * UAL_R_HAL;                                RigidTrans = [RigidTrans; {'HAL',Glob_R_HAL}];
            
            
            % Table of RAMSIS-RPx names
%             BodyTranslator( 1,:) = {'pelvis','BEC'};
%             BodyTranslator( 2,:) = {'lower lumbar spine','ULW'};
%             BodyTranslator( 3,:) = {'upper lumbar spine','OLW'};
%             BodyTranslator( 4,:) = {'lower thoracic spine','UBW'};
%             BodyTranslator( 5,:) = {'upper thoracic spine','OBW'};
%             BodyTranslator( 6,:) = {'lower cervical spine','UHW'};
%             BodyTranslator( 7,:) = {'upper cervical spine','OHW'};
%             BodyTranslator( 8,:) = {'head','KO'};
%             BodyTranslator( 9,:) = {'chest','BRK'};
%             BodyTranslator(10,:) = {'right clavicle','SBR'};
%             BodyTranslator(11,:) = {'left clavicle','SBL'};
%             BodyTranslator(12,:) = {'right upper arm','OAR'};
%             BodyTranslator(13,:) = {'left upper arm','OAL'};
%             BodyTranslator(14,:) = {'right lower arm','UAR'};
%             BodyTranslator(15,:) = {'left lower arm','UAL'};
%             BodyTranslator(16,:) = {'right hand','HAR'};
%             BodyTranslator(17,:) = {'left hand','HAL'};
%             BodyTranslator(18,:) = {'right thigh','OSR'};
%             BodyTranslator(19,:) = {'left thigh','OSL'};
%             BodyTranslator(20,:) = {'right lower leg','USR'};
%             BodyTranslator(21,:) = {'left lower leg','USL'};
%             BodyTranslator(22,:) = {'right foot','FUR'};
%             BodyTranslator(23,:) = {'left foot','FUL'};
%             BodyTranslator(24,:) = {'right ball of foot','FBR'};
%             BodyTranslator(25,:) = {'left ball of foot','FBL'};
%             % PUFO. Fingers are not currently in the model.
%             % Markers on fingers are added to body Hand
%             BodyTranslator(26,:) = {'right thumb finger 3','HAR'};
%             BodyTranslator(27,:) = {'right index finger 3','HAR'};
%             BodyTranslator(28,:) = {'right middle finger 3','HAR'};
%             BodyTranslator(29,:) = {'left index finger 3','HAL'};
            
            HM.Segments(2).LocalPoints(1).Point.GlobalCoord = Glob_Pos_GHZ;
            NSegments = size(HM.Segments,1);
            for i=1:NSegments
                if HM.Segments(i).Fixed ~= 1 % If is not Ground
                %if ~strcmpi(HM.Segments(i).Name,'Ground')
                    NPoints = size(HM.Segments(i).LocalPoints,1);
                    NMarkers = size(HM.Segments(i).LocalMarkers,1);
                    HM.Segments(i).LocalPoints(1).LocCoord = [0;0;0];
                    RIndex = getCellIndex(RigidTrans,HM.Segments(i).Name);
                    for j=2:NPoints
%                         PName = HM.Segments(i).LocalPoints(j).Point.Name;
%                         HM.Segments(i).LocalPoints(j).LocCoord = (StructMat.(PName))';
                        HM.Segments(i).LocalPoints(j).Point.GlobalCoord = RigidTrans{RIndex,2}*1000*HM.Segments(i).LocalPoints(j).LocCoord + HM.Segments(i).LocalPoints(1).Point.GlobalCoord;
                    end
                    for j=1:NMarkers
                        HM.Segments(i).LocalMarkers(j).Point.GlobalCoord = RigidTrans{RIndex,2}*1000*HM.Segments(i).LocalMarkers(j).LocCoord + HM.Segments(i).LocalPoints(1).Point.GlobalCoord;                   
                    end
                end
            end
%             NMarkers = length(StructMat.Marker);
%             for i=1:NMarkers
%                 MarkerName   = upper(StructMat.Marker(i).Name);
%                 SegmentIndex = getCellIndex(BodyTranslator,StructMat.Marker(i).BodyName);
%                 SegmentName  = BodyTranslator{SegmentIndex,2};
%                 SegmentIndex = getVecIndex(SegmentName,HM.Segments);
%                 if SegmentIndex>0
%                     RIndex = getCellIndex(RigidTrans,HM.Segments(SegmentIndex).Name);
%                     MarkerIndex = getVecIndex(MarkerName,HM.Segments(SegmentIndex).LocalMarkers);
%                     HM.Segments(SegmentIndex).LocalMarkers(MarkerIndex).LocCoord = StructMat.Marker(i).LocalPosition';
%                     HM.Segments(SegmentIndex).LocalMarkers(MarkerIndex).Point.GlobalCoord = RigidTrans{RIndex,2}*HM.Segments(SegmentIndex).LocalMarkers(MarkerIndex).LocCoord + HM.Segments(SegmentIndex).LocalPoints(1).Point.GlobalCoord;                   
%                 else
% %                     warning('Body %s is not in Human Model',SegmentName);
%                 end
%             end
            
            FileNameExt = [[FilePath,FileName(1:end-4)],'.EXP'];
            
            
            
            % open file
            fid = fopen(FileNameExt,'w');
            % Header
            fprintf(fid,'RAMSIS_MODEL_DATA_EXPORT_FILE V2.00\n');
            % Points
            fprintf(fid,'SECTION SKELETON_POINTS\n');
            NPoints = size(HM.Points,1);
            for i=2:NPoints
                fprintf(fid,['\t',HM.Points(i).Name,'\t\t\t hip-center\t\t',num2str(HM.Points(i).GlobalCoord','%-11.3f'),'\n']);
            end
            fprintf(fid,'ENDSECTION SKELETON_POINTS\n');
            %             % Joints
            fprintf(fid,'SECTION JOINTS\n');
            for i=2:NPoints
                fprintf(fid,['\t',HM.Points(i).Name,' \t',num2str(HM.Points(i).GlobalCoord','%-11.3f'),'\t\t',num2str(Angle.([HM.Points(i).Name,'_a1'])*180/pi,'%-8.3f'),'\t',...
                num2str(Angle.([HM.Points(i).Name,'_a2'])*180/pi,'%-8.3f'),'\t',num2str(Angle.([HM.Points(i).Name,'_a3'])*180/pi,'%-8.3f'),'\n']);
            end
            fprintf(fid,'ENDSECTION JOINTS\n');
%             % Joint System
            fprintf(fid,'SECTION JOINTS_COORDINATE_SYSTEMS\n');
            fprintf(fid,'    GHZ	     0.000      0.000      0.000	     0.000      0.000      1.000	    -1.000      0.000      0.000	     0.000     -1.000      0.000\n');
%             Glob_R_GLK = Glob_R_BEC*BEC_R_GLK;
%             fprintf(fid,['\tGLK \t',num2str(Glob_Pos_GLK(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GLK(:,1)','%-8.3f'),'\t',...
%                 num2str(Glob_R_GLK(:,2)','%-8.3f'),'\t',num2str(Glob_R_GLK(:,3)','%-8.3f'),'\n']);
%             Glob_R_GHUL = Glob_R_BEC*BEC_R_GHUL;
%             fprintf(fid,['\tGHUL \t',num2str(Glob_Pos_GHUL(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GHUL(:,1)','%-8.3f'),'\t',...
%                 num2str(Glob_R_GHUL(:,2)','%-8.3f'),'\t',num2str(Glob_R_GHUL(:,3)','%-8.3f'),'\n']);
%             Glob_R_GKNL = Glob_R_OSL*OSL_R_GKNL;
%             fprintf(fid,['\tGKNL \t',num2str(Glob_Pos_GKNL(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GKNL(:,1)','%-8.3f'),'\t',...
%                 num2str(Glob_R_GKNL(:,2)','%-8.3f'),'\t',num2str(Glob_R_GKNL(:,3)','%-8.3f'),'\n']);
%             Glob_R_GSPL = Glob_R_USL*USL_R_GSPL;
%             fprintf(fid,['\tGSPL \t',num2str(Glob_Pos_GSPL(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GSPL(:,1)','%-8.3f'),'\t',...
%                 num2str(Glob_R_GSPL(:,2)','%-8.3f'),'\t',num2str(Glob_R_GSPL(:,3)','%-8.3f'),'\n']);
%             Glob_R_GFBL = Glob_R_FUL*FUL_R_GFBL;
%             fprintf(fid,['\tGFBL \t',num2str(Glob_Pos_GFBL(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_GFBL(:,1)','%-8.3f'),'\t',...
%                 num2str(Glob_R_GFBL(:,2)','%-8.3f'),'\t',num2str(Glob_R_GFBL(:,3)','%-8.3f'),'\n']);
            
            fprintf(fid,'ENDSECTION JOINTS_COORDINATE_SYSTEMS\n');
%             % Body System
            fprintf(fid,'SECTION BODY_ELEMENT_COORDINATE_SYSTEMS\n');
            for i=1:NSegments
                if HM.Segments(i).Fixed ~= 1 % If is not Fixed(Ground)                    
                    JointName = HM.Segments(i).LocalPoints(1).Name;
                    Global_Pos_Joint = HM.Segments(i).LocalPoints(1).Point.GlobalCoord;
                    RIndex = getCellIndex(RigidTrans,HM.Segments(i).Name);
                    Glob_Vec_U = RigidTrans{RIndex,2}(:,1);
                    Glob_Vec_V = RigidTrans{RIndex,2}(:,2);
                    Glob_Vec_W = RigidTrans{RIndex,2}(:,3);
                    fprintf(fid,['\t',JointName,' \t',num2str(Global_Pos_Joint','%-11.3f'),'\t\t',num2str(Glob_Vec_U','%-8.3f'),'\t',...
                        num2str(Glob_Vec_V','%-8.3f'),'\t',num2str(Glob_Vec_W','%-8.3f'),'\n']);
                end
            end
%             fprintf(fid,['\tGHZ \t',num2str(Glob_Pos_GHZ(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_BEC(:,1)','%-8.3f'),'\t',...
%                 num2str(Glob_R_BEC(:,2)','%-8.3f'),'\t',num2str(Glob_R_BEC(:,3)','%-8.3f'),'\n']);
%             fprintf(fid,['\tGHUL \t',num2str(Glob_Pos_GHUL(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_OSL(:,1)','%-8.3f'),'\t',...
%                 num2str(Glob_R_OSL(:,2)','%-8.3f'),'\t',num2str(Glob_R_OSL(:,3)','%-8.3f'),'\n']);
%             fprintf(fid,['\tGKNL \t',num2str(Glob_Pos_GKNL(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_USL(:,1)','%-8.3f'),'\t',...
%                 num2str(Glob_R_USL(:,2)','%-8.3f'),'\t',num2str(Glob_R_USL(:,3)','%-8.3f'),'\n']);
%             fprintf(fid,['\tGSPL \t',num2str(Glob_Pos_GSPL(1:3)','%-11.3f'),'\t\t',num2str(Glob_R_FUL(:,1)','%-8.3f'),'\t',...
%                 num2str(Glob_R_FUL(:,2)','%-8.3f'),'\t',num2str(Glob_R_FUL(:,3)','%-8.3f'),'\n']);
            fprintf(fid,'ENDSECTION BODY_ELEMENT_COORDINATE_SYSTEMS\n'); 
            % Skin_Points
            fprintf(fid,'SECTION SKIN_POINTS\n');
            fprintf(fid,'   BEC_1_1	  -931.962    764.185    779.167\n');
            fprintf(fid,'ENDSECTION SKIN_POINTS\n');
%             % Marker
            fprintf(fid,'SECTION MARKER\n');
%             SegmentTranslator(1,:) = {'Pelvis','BEC',Glob_T_BEC};
%             SegmentTranslator(2,:) = {'LeftThigh','OSL',Glob_T_OSL};
%             SegmentTranslator(3,:) = {'LeftShank','USL',Glob_T_USL};
%             SegmentTranslator(4,:) = {'LeftFoot','FUL',Glob_T_FUL};
            for i=1:NSegments
                if HM.Segments(i).Fixed ~= 1 % If is not Fixed(Ground)                    
                    NMarkers = size(HM.Segments(i).LocalMarkers,1);
                    SegName = HM.Segments(i).Name;
                    %SegIndex = getCellIndex(SegmentTranslator,SegName);
                    for j=1:NMarkers
                        MarkerName = HM.Segments(i).LocalMarkers(j).Point.Name;
                        fprintf(fid,['\t',MarkerName,'\t\t',SegName,'\t\t\t',num2str(HM.Segments(i).LocalMarkers(j).Point.GlobalCoord','%-11.3f'),'\n']);
                    end
                end
            end
            
            fprintf(fid,'ENDSECTION MARKER\n');
%             % Centres of gravity
            fprintf(fid,'SECTION CENTRES_OF_GRAVITY\n');
%             for i=1:NSegments
%                 SegName = HM.Segments(i).Name;
%                 SegIndex = getCellIndex(SegmentTranslator,SegName);
%                 Loc_CoM = HM.Segments(i).CoM.LocCoord;
%                 Glob_CoM = SegmentTranslator{SegIndex,3}*HM.Segments(i).RamsisLCS_T_ALLCS*[Loc_CoM;1];
                fprintf(fid,'      STE	     0.000      0.000      0.000\n');
%             end
            fprintf(fid,'ENDSECTION CENTRES_OF_GRAVITY\n');
%             % Mass
            fprintf(fid,'SECTION MASS\n');
%             for i=1:NSegments
%                 SegName = HM.Segments(i).Name;
%                 SegIndex = getCellIndex(SegmentTranslator,SegName);
                fprintf(fid,'      STE	     0.000\n');
%             end
            fprintf(fid,'ENDSECTION MASS\n');
            fprintf(fid,'SECTION MOMENT_OF_INERTIA\n');
%             for i=1:NSegments
%                 SegName = HM.Segments(i).Name;
%                 SegIndex = getCellIndex(SegmentTranslator,SegName);
%                 ALLCS_I_CM = HM.Segments(i).I;
%                 RamsisLCS_I_CM = HM.Segments(i).RamsisLCS_T_ALLCS(1:3,1:3) * ALLCS_I_CM * HM.Segments(i).RamsisLCS_T_ALLCS(1:3,1:3)';
%                 fprintf(fid,['\t',SegmentTranslator{SegIndex,2},'\t\t',...
%                     num2str(RamsisLCS_I_CM(1,1),'%-11.3f'),'\t\t',num2str(RamsisLCS_I_CM(1,2),'%-11.3f'),'\t\t',num2str(RamsisLCS_I_CM(1,3),'%-11.3f'),'\t\t',...
%                     num2str(RamsisLCS_I_CM(2,1),'%-11.3f'),'\t\t',num2str(RamsisLCS_I_CM(2,2),'%-11.3f'),'\t\t',num2str(RamsisLCS_I_CM(2,3),'%-11.3f'),'\t\t',...
%                     num2str(RamsisLCS_I_CM(3,1),'%-11.3f'),'\t\t',num2str(RamsisLCS_I_CM(3,2),'%-11.3f'),'\t\t',num2str(RamsisLCS_I_CM(3,3),'%-11.3f'),'\t\t','\n']);
                fprintf(fid,'	FUL		255.472		-88.399		31.823		-88.399		1145.645		-14.144		31.823		-14.144		1082.882\n');
%             end
            fprintf(fid,'ENDSECTION MOMENT_OF_INERTIA\n');
            % close file
            fclose(fid);
        end
        function writePSP(HM,FileName, SubjectPath, ExpPath, PAM_OutputPath, PAM_DeskPath)
            
            disp([''])
            % Change current working directory
            CurrentDir = pwd;
            cd(PAM_DeskPath); % for default soft installation is PAM_DeskPath='C:\Program Files\ESI Group\PAM-DHErgo\2010.1\Binary\Desk';

            % Define Two environment variables needed by PAM software
            % 1) Path with the experimental data
            setenv('DHERGO_EXPDIR', ExpPath);
            % 2) Path for PAM intermediate files
            setenv('DHERGO_PAMWORKDIR', PAM_OutputPath);
            
            % write subject parameters (.psp) for PAM model
            LogFilePSP = [FileName,'_PSPout','.log'];
            dos(['pamcomfort-client -CreateSubject HM50KR 01_SD > ',LogFilePSP,' 2>&1']);  % funciona
            dos(['move ',LogFilePSP,' "',SubjectPath,'"']);

            % go back to original working directory
            cd(CurrentDir);   
        end

    end
    
end

