classdef HUMAN_PARSER_PSP< handle
    %HUMAN_PARSER_PSP 
    
    properties
        PspDom          % Document Object Model with the definition of psp file
        mm2m = 0.001;    % Constant
        deg2rad=pi/180; % Constant
    end
    
    methods
        function HPpsp = HUMAN_PARSER_PSP() 
        end
        function Sement_R_Joint = calRMatrix(HPpsp,Looc_Vec_x,Looc_Vec_y)
            Looc_Vec_z = cross(Looc_Vec_x, Looc_Vec_y);
            Looc_Vec_y = cross(Looc_Vec_z, Looc_Vec_x);
            Looc_Vec_x = Looc_Vec_x / norm(Looc_Vec_x);
            Looc_Vec_y = Looc_Vec_y / norm(Looc_Vec_y);
            Looc_Vec_z = Looc_Vec_z / norm(Looc_Vec_z);
            Sement_R_Joint = [Looc_Vec_x, Looc_Vec_y, Looc_Vec_z];
        end
        function parseSegmentsInfo(HPpsp,Human)
            SegmentTranslator(1,:) = {'Pelvis','HM50KR_posture'};
            SegmentTranslator(2,:) = {'LeftThigh','Left_leg'};
            SegmentTranslator(3,:) = {'LeftShank','Left_calf'};
            SegmentTranslator(4,:) = {'LeftFoot','Left_foot'};
            NSegments = size(Human.Segments,1);
            % Recorremos todos los segmento del archivo psp
            PspSystem = HPpsp.PspDom.getElementsByTagName('System');
%             PspSegments = PspSystem.item(0).getElementsByTagName('Segment');
            PspSegments = HPpsp.PspDom.getElementsByTagName('Segment');
            for i = 1:PspSegments.getLength
                SegmentIndex = 0;
                SegmentListItem = PspSegments.item(i-1);
                SegmentName = char(SegmentListItem.getAttribute('Name')); 
                SegmentMass = char(SegmentListItem.getAttribute('Mass'));
                
%                 SegmentCoM = char(SegmentListItem.getAttribute('CoM'));
                SegmentI = char(SegmentListItem.getAttribute('I'));
                
                    
                SegmentJoints = SegmentListItem.getElementsByTagName('Joint');
                SegmentPoints = SegmentListItem.getElementsByTagName('Point');
                JointListItem = SegmentJoints.item(0);
                if ~isempty(JointListItem)
                    JointLCS = JointListItem.getElementsByTagName('LocalCoordSystem');
                    JointLCSItem = JointLCS.item(0);
                    Topology = 0;
                    if isempty (JointLCSItem)
                        Topology = 1;
                    else
                       NoNum = str2num(char(JointLCSItem.getAttribute('Origin')));
                       if isempty(NoNum)
                           Topology = 1;
                       end
                    end
                end
                if ~isempty(SegmentMass)
                    SegIndex = getCellIndex(SegmentTranslator,SegmentName);
                    SegHuIndex = getVecIndex(SegmentTranslator(SegIndex,2),Human.Segments);
                    Human.Segments(SegHuIndex).Mass = str2num(SegmentMass);
                    % Get CoM
                    SegmentCoM = SegmentListItem.getElementsByTagName('CoM');
                    SCoM = SegmentCoM.item(0);
                    SNCoM = char(SCoM.getAttribute('XYZ_in_AL-LCS'));
                    % Add CoM in list of points in Human
                    AddedPoint = Human.addCoM([SegmentTranslator(SegIndex,2),'_CoM']);
                    NModelCoM = size(Human.CoMs,1);
                    % Add CoM in Segment
                    Human.Segments(SegHuIndex).CoM = LOCAL_POINT(AddedPoint,[SegmentTranslator(SegIndex,2),'_CoM']);
                    Human.Segments(SegHuIndex).CoM.LocCoord = HPpsp.mm2m * str2num(SNCoM);
                    % Get MoI
                    SegmentMoI = SegmentListItem.getElementsByTagName('MoI');
                    SMoI = SegmentMoI.item(0);
                    SIxx = str2num(char(SMoI.getAttribute('Ixx')));
                    SIxy = str2num(char(SMoI.getAttribute('Ixy')));
                    SIxz = str2num(char(SMoI.getAttribute('Ixz')));
                    SIyx = str2num(char(SMoI.getAttribute('Iyx')));
                    SIyy = str2num(char(SMoI.getAttribute('Iyy')));
                    SIyz = str2num(char(SMoI.getAttribute('Iyz')));
                    SIzx = str2num(char(SMoI.getAttribute('Izx')));
                    SIzy = str2num(char(SMoI.getAttribute('Izy')));
                    SIzz = str2num(char(SMoI.getAttribute('Izz')));
                    SegmentI = [SIxx,SIxy,SIxz; SIyx,SIyy,SIyz; SIzx,SIzy,SIzz];
                    Human.Segments(SegHuIndex).I = SegmentI*0.0001;% Pass from kg·cm2 to kg·m2;
                end
                    
                % Calculamos la posición del segmento psp correspondiente en nuestra clase
                for j=1:NSegments
                    ClassSName = regexprep(Human.Segments(j).Name,'_',' ');
                    if strcmpi(ClassSName,SegmentName) && Topology==0
                        SegmentIndex = getVecIndex(Human.Segments(j).Name,Human.Segments);
                        break;
                    end
                end
                if SegmentIndex ~=0
%                     if ~isempty(SegmentMass)
%                         Human.Segments(SegmentIndex).Mass = str2num(SegmentMass);
%                     end
%                     if ~isempty(SegmentCoM)
%                         % Add CoM in list of points in Human
%                         AddedPoint = Human.addCoM([SegmentName,'_CoM']);
%                         NModelCoM = size(Human.CoMs,1);
%                         % Add CoM in Segment
%                         Human.Segments(SegmentIndex).CoM = LOCAL_POINT(AddedPoint,[SegmentName,'_CoM']);
%                         Human.Segments(SegmentIndex).CoM.LocCoord = HPpsp.mm2m * str2num(SegmentCoM);
%                     end
%                     if ~isempty(SegmentI)
%                         Human.Segments(SegmentIndex).I = str2num(SegmentI)*0.00001;% Pass from kg·cm2 to kg·m2;
%                     end
                    % Recorremos los joints de cada segmento
                    for j=1:SegmentJoints.getLength
                        JointListItem = SegmentJoints.item(j-1);
                        JointName = char(JointListItem.getAttribute('Mnemonic'));
                        JointLCS = JointListItem.getElementsByTagName('LocalCoordSystem');
                        JointLCSItem = JointLCS.item(0);
                        NLPoints = size(Human.Segments(SegmentIndex).LocalPoints,1);
                        NLVectors = size(Human.Segments(SegmentIndex).LocalVectors,1);
                        % Obtenemos la posición local del punto mirando en el origen del joint
                        for k=1:NLPoints
                            if strcmpi(Human.Segments(SegmentIndex).LocalPoints(k).Point.Name,JointName)
                                PointLoocCoord = str2num(char(JointLCSItem.getAttribute('Origin')));
                                Human.Segments(SegmentIndex).LocalPoints(k).LocCoord = HPpsp.mm2m * PointLoocCoord';
                                break;
                            end
                        end
                        % Obtenemos las coordenadas locales de los vectores de los joints que tengamos
                        for k=1:NLVectors
                            if strcmpi(Human.Segments(SegmentIndex).LocalVectors(k).Vector.Name(1:(end-3)),JointName)
                                if strcmpi(Human.Segments(SegmentIndex).LocalVectors(k).Vector.Name((end-1)),'U')
                                    if isempty(Human.Segments(SegmentIndex).LocalVectors(k).LocCoord)
                                        VLoocCoord = str2num(char(JointLCSItem.getAttribute('U')));
                                        Human.Segments(SegmentIndex).LocalVectors(k).LocCoord = VLoocCoord';
                                    end
                                elseif strcmpi(Human.Segments(SegmentIndex).LocalVectors(k).Vector.Name((end-1)),'V')
                                    if isempty(Human.Segments(SegmentIndex).LocalVectors(k).LocCoord)
                                        VLoocCoord = str2num(char(JointLCSItem.getAttribute('V')));
                                        Human.Segments(SegmentIndex).LocalVectors(k).LocCoord = VLoocCoord';
                                    end
                                elseif strcmpi(Human.Segments(SegmentIndex).LocalVectors(k).Vector.Name((end-1)),'W')
                                    if isempty(Human.Segments(SegmentIndex).LocalVectors(k).LocCoord)
                                        VLoocCoord = str2num(char(JointLCSItem.getAttribute('W')));
                                        Human.Segments(SegmentIndex).LocalVectors(k).LocCoord = VLoocCoord';
                                    end
                                end
                            end
                        end
                    end
                    % Recorremos los points de cada segmento
                    for j=1:SegmentPoints.getLength
                        PointListItem = SegmentPoints.item(j-1);
                        PointName = char(PointListItem.getAttribute('Mnemonic'));
                        PointType = char(PointListItem.getAttribute('Type'));
                        % Si el segmento es sólo distal no joint obtenemos sus coordenadas locales
                        if strcmpi(PointType,'Distal')
                            for k=1:NLPoints
                                if strcmpi(Human.Segments(SegmentIndex).LocalPoints(k).Point.Name,PointName)
                                    if isempty(Human.Segments(SegmentIndex).LocalPoints(k).LocCoord)
                                        PointLoocCoord = str2num(char(PointListItem.getAttribute('LocalPosition')));
                                        Human.Segments(SegmentIndex).LocalPoints(k).LocCoord = HPpsp.mm2m *PointLoocCoord';
                                        break;
                                    end
                                end
                            end
                            % Obtenemos los markers y sus coordenadas locales en el segmento
                        elseif strcmpi(PointType,'External Marker')
                            % add marker in list of markers in human
                            AddedMarker = Human.addMarker(PointName);
                            
                            %  Add Local Point to Segment and its local coordinates
                            Human.Segments(SegmentIndex).addMarker(AddedMarker);
                            NLocalMarker = size(Human.Segments(SegmentIndex).LocalMarkers,1);
                            MarkerLoocCoord = str2num(char(PointListItem.getAttribute('LocalPosition')));
                            Human.Segments(SegmentIndex).LocalMarkers(NLocalMarker).LocCoord = HPpsp.mm2m * MarkerLoocCoord';
                            
                        end
                    end
                end
                    
            end 
             % add marker in q
            Human.addMarkerInq();
        end
        function parseSensor(HPpsp,Human)
            NSensors = size(Human.Sensors,1);
            PspSegments = HPpsp.PspDom.getElementsByTagName('Segment');
            for i = 1:PspSegments.getLength
                SegmentIndex = 0;
                SegmentListItem = PspSegments.item(i-1);
                SegmentName = char(SegmentListItem.getAttribute('Name'));
                SegmentJoints = SegmentListItem.getElementsByTagName('Joint');
                JointListItem = SegmentJoints.item(0);
                if ~isempty(JointListItem)
                    JointLCS = JointListItem.getElementsByTagName('LocalCoordSystem');
                    JointLCSItem = JointLCS.item(0);
                    Topology = 0;
                    if isempty (JointLCSItem)
                        Topology = 1;
                    end
                end
                for j=1:NSensors
                    if ~strcmpi(Human.Sensors{j}.Type,'TRS')
                        SenSeg1Name = regexprep(Human.Sensors{j}.Segment1.Name,'_',' ');
                        SenSeg2Name = regexprep(Human.Sensors{j}.Segment2.Name,'_',' ');
                        if strcmpi(SenSeg1Name,SegmentName)&& Topology==0
                            SegmentIndex = getVecIndex(Human.Sensors{j}.Segment1.Name,Human.Segments);
                            Sensor1Pos = j; % Hacer algo cuando sean más de uno ejem pelvis
                            %                         break;
                        elseif strcmpi(SenSeg2Name,SegmentName)&& Topology==0
                            SegmentIndex = getVecIndex(Human.Sensors{j}.Segment2.Name,Human.Segments);
                            Sensor2Pos = j;
                        end
                    end
                end
                if SegmentIndex ~=0
                    SegmentJoints = SegmentListItem.getElementsByTagName('Joint');
                    for j=1:SegmentJoints.getLength
                        JointListItem = SegmentJoints.item(j-1);
                        JointName = char(JointListItem.getAttribute('Mnemonic'));
                        JointLCS = JointListItem.getElementsByTagName('LocalCoordSystem');
                        JointLCSItem = JointLCS.item(0);
                        ULoocCoord = str2num(char(JointLCSItem.getAttribute('U')));
                        VLoocCoord = str2num(char(JointLCSItem.getAttribute('V')));
%                         WLoocCoord = str2num(char(JointLCSItem.getAttribute('W')));
                        if strcmpi(JointName,Human.Sensors{Sensor1Pos}.Name)&& strcmpi(SegmentName,regexprep(Human.Sensors{Sensor1Pos}.Segment1.Name,'_',' '))
                            if ischar(Human.Sensors{Sensor1Pos}.Perm1x)
                                if ~isempty(ULoocCoord)
                                    Human.Sensors{Sensor1Pos}.Perm1x = ULoocCoord';
                                end
                            end
                            if ischar(Human.Sensors{Sensor1Pos}.Perm1y)
                                 if ~isempty(VLoocCoord)
                                    Human.Sensors{Sensor1Pos}.Perm1y = VLoocCoord';
                                 end
                            end
                        elseif strcmpi(JointName,Human.Sensors{Sensor2Pos}.Name)&& strcmpi(SegmentName,regexprep(Human.Sensors{Sensor2Pos}.Segment2.Name,'_',' '))
                            if ischar(Human.Sensors{Sensor2Pos}.Perm2x)
                                if ~isempty(ULoocCoord)
                                    Human.Sensors{Sensor2Pos}.Perm2x = ULoocCoord';
                                end
                            end
                            if ischar(Human.Sensors{Sensor2Pos}.Perm2y)
                                if ~isempty(VLoocCoord)
                                    Human.Sensors{Sensor2Pos}.Perm2y = VLoocCoord';
                                end
                            end
                        end
                    end
                end
            end
        end
        function parseXMat(HPpsp,Human)
            
%             Angle.HPT_a1  =  0*HPpsp.deg2rad; 
%             Angle.HPT_a2  =  0*HPpsp.deg2rad; 
%             Angle.HPT_a3  =  0*HPpsp.deg2rad;
%             Angle.LHJ_a3 =   0*HPpsp.deg2rad; 
%             Angle.LHJ_a2 =  90*HPpsp.deg2rad; 
%             Angle.LHJ_a1 =   0*HPpsp.deg2rad;
%             Angle.LKJ_a1 =   0*HPpsp.deg2rad;
%             Angle.LKJ_a2 = -90*HPpsp.deg2rad;
%             Angle.LKJ_a3 = 180*HPpsp.deg2rad;
%             Angle.LAJ_a1 = -90*HPpsp.deg2rad;
%             Angle.LAJ_a2 =   0*HPpsp.deg2rad;
%             Angle.LAJ_a3 = -90*HPpsp.deg2rad;
%             
%             Glob_R_Pelvis = rot321s(Angle.HPT_a3,Angle.HPT_a2,Angle.HPT_a1);
%             LHJ_R_LThihg = rot321s(Angle.LHJ_a3,Angle.LHJ_a2,Angle.LHJ_a1);
%             LKJ_R_LShank = rot321s(Angle.LKJ_a3,Angle.LKJ_a2,Angle.LKJ_a1);
%             LAJ_R_LFoot = rot321s(Angle.LAJ_a3,Angle.LAJ_a2,Angle.LAJ_a1);
            
%             Angle.HPT_a1  =  0*HPpsp.deg2rad; 
%             Angle.HPT_a2  =  21*HPpsp.deg2rad; 
%             Angle.HPT_a3  =  6*HPpsp.deg2rad;
%             Angle.LHJ_a3 =   0*HPpsp.deg2rad; 
%             Angle.LHJ_a2 =  -15*HPpsp.deg2rad; 
%             Angle.LHJ_a1 =   45*HPpsp.deg2rad;
%             Angle.LKJ_a1 =   0*HPpsp.deg2rad;
%             Angle.LKJ_a2 = -70*HPpsp.deg2rad;
%             Angle.LKJ_a3 = 20*HPpsp.deg2rad;
%             Angle.LAJ_a1 = -16*HPpsp.deg2rad;
%             Angle.LAJ_a2 =   0*HPpsp.deg2rad;
%             Angle.LAJ_a3 = 0*HPpsp.deg2rad;
            
            Angle.HPT_a1  = 0*HPpsp.deg2rad; 
            Angle.HPT_a2  = 0*HPpsp.deg2rad; 
            Angle.HPT_a3  = 90*HPpsp.deg2rad;
            Angle.LHJ_a1 = 60*HPpsp.deg2rad; 
            Angle.LHJ_a2 = 0*HPpsp.deg2rad; 
            Angle.LHJ_a3 = 0*HPpsp.deg2rad;
            Angle.LKJ_a1 = -80*HPpsp.deg2rad;
            Angle.LKJ_a2 = 0*HPpsp.deg2rad;
            Angle.LKJ_a3 = -20*HPpsp.deg2rad;
            Angle.LAJ_a1 = 0*HPpsp.deg2rad;
            Angle.LAJ_a2 = 0*HPpsp.deg2rad;
            Angle.LAJ_a3 = 0*HPpsp.deg2rad;
            
            Glob_R_HPT = rot123s(Angle.HPT_a1,Angle.HPT_a2,Angle.HPT_a3);
            LHJ_R_LHJ = rot123s(Angle.LHJ_a1,Angle.LHJ_a2,Angle.LHJ_a3);
            LKJ_R_LKJ = rot123s(Angle.LKJ_a1,Angle.LKJ_a2,Angle.LKJ_a3);
            LAJ_R_LAJ = rot123s(Angle.LAJ_a1,Angle.LAJ_a2,Angle.LAJ_a3);            
            % Recorremos todos los segmento del archivo psp
            PspSegments = HPpsp.PspDom.getElementsByTagName('Segment');
            for i = 1:PspSegments.getLength
                SegmentListItem = PspSegments.item(i-1);
                SegmentName = char(SegmentListItem.getAttribute('Name'));
                SegmentJoints = SegmentListItem.getElementsByTagName('Joint');
                for j=1:SegmentJoints.getLength
                    JointListItem = SegmentJoints.item(j-1);
                    JointName = char(JointListItem.getAttribute('Mnemonic'));
                    JointLCS = JointListItem.getElementsByTagName('LocalCoordSystem');
                    JointLCSItem = JointLCS.item(0);
                    if ~isempty(JointLCSItem)
                        ULoocCoord = str2num(char(JointLCSItem.getAttribute('U')));
                        VLoocCoord = str2num(char(JointLCSItem.getAttribute('V')));
                        WLoocCoord = str2num(char(JointLCSItem.getAttribute('W')));
                        if strcmpi(SegmentName,'HM50KR posture')&& strcmpi(JointName,'HPT')
                            Pelvis_R_HPT = HPpsp.calRMatrix(ULoocCoord',VLoocCoord');
                        elseif strcmpi(SegmentName,'HM50KR posture')&& strcmpi(JointName,'LHJ')
                            if ~isempty(ULoocCoord)
                                Pelvis_R_LHJ = HPpsp.calRMatrix(ULoocCoord',VLoocCoord');
                            end
                        elseif strcmpi(SegmentName,'Left leg')&& strcmpi(JointName,'LKJ')
                             if ~isempty(ULoocCoord)
                                 LThihg_R_LKJ = HPpsp.calRMatrix(ULoocCoord',VLoocCoord');
                             end
                        elseif strcmpi(SegmentName,'Left calf')&& strcmpi(JointName,'LAJ')
                             if ~isempty(ULoocCoord)
                                 LShank_R_LAJ = HPpsp.calRMatrix(ULoocCoord',VLoocCoord');
                             end
                        end
                    end
                end
            end
            P1 =[0,-1,0;0,0,-1;1,0,0];
            P2 =[0,0,1;-1,0,0;0,-1,0];
%             Glob_R_Pelvis = Glob_R_Pelvis* Pelvis_R_HPT;                                    RigidTrans(1,:) = {'HM50KR_posture',Glob_R_Pelvis};
%             Glob_R_LThihg = Glob_R_Pelvis * Pelvis_R_LHJ* Pelvis_R_LThihg *LThihg_R_LHJ;    RigidTrans = [RigidTrans; {'Left_leg',Glob_R_LThihg}];
%             Glob_R_LShank = Glob_R_LThihg * LThihg_R_LKJ * LThihg_R_LShank * LShank_R_LKJ;  RigidTrans = [RigidTrans; {'Left_calf',Glob_R_LShank}];
%             Glob_R_LFoot  = Glob_R_LShank * LShank_R_LAJ * LShank_R_LFoot * LFoot_R_LAJ;    RigidTrans = [RigidTrans; {'Left_foot',Glob_R_LFoot}];
            
%             Glob_R_Pelvis = Glob_R_Pelvis* Pelvis_R_HPT;                     RigidTrans(1,:) = {'HM50KR_posture',Glob_R_Pelvis};
%             Glob_R_Pelvis = Glob_R_HPT* Pelvis_R_HPT;                       RigidTrans(1,:) = {'HM50KR_posture',Glob_R_Pelvis};
%             Glob_R_LThihg = Glob_R_Pelvis * Pelvis_R_LHJ * LHJ_R_LThihg;   RigidTrans = [RigidTrans; {'Left_leg',Glob_R_LThihg}];
% %             Glob_R_LShank = Glob_R_LThihg * LThihg_R_LKJ * LKJ_R_LShank * LThihg_R_LKJ';  RigidTrans = [RigidTrans; {'Left_calf',Glob_R_LShank}];
%             Glob_R_LShank = Glob_R_LThihg * LThihg_R_LKJ * LKJ_R_LShank;  RigidTrans = [RigidTrans; {'Left_calf',Glob_R_LShank}];
%             Glob_R_LFoot  = Glob_R_LShank * LShank_R_LAJ * LAJ_R_LFoot;   RigidTrans = [RigidTrans; {'Left_foot',Glob_R_LFoot}];
            
            Glob_R_Pelvis = Glob_R_HPT* Pelvis_R_HPT;                       RigidTrans(1,:) = {'HM50KR_posture',Glob_R_Pelvis};
            Glob_R_LThihg = Glob_R_Pelvis * P1*LHJ_R_LHJ * P2;   RigidTrans = [RigidTrans; {'Left_leg',Glob_R_LThihg}];
%             Glob_R_LShank = Glob_R_LThihg * LThihg_R_LKJ * LKJ_R_LShank * LThihg_R_LKJ';  RigidTrans = [RigidTrans; {'Left_calf',Glob_R_LShank}];
            Glob_R_LShank = Glob_R_LThihg * P1 *LKJ_R_LKJ * P2;  RigidTrans = [RigidTrans; {'Left_calf',Glob_R_LShank}];
            Glob_R_LFoot  = Glob_R_LShank * P1 * LAJ_R_LAJ * P2;   RigidTrans = [RigidTrans; {'Left_foot',Glob_R_LFoot}];NTrans = size(RigidTrans);
            for i=1:NTrans
                SegmentIndex = getVecIndex(RigidTrans{i,1},Human.Segments);
                if SegmentIndex ~= 0
                    Human.Segments(SegmentIndex).Glob_R_Seg =  RigidTrans{i,2};
                end
            end
        end
        function readxml(HPpsp,PathXML,FileXML)
            % READXML Read the psp file in xml format
            HPpsp.PspDom = xmlread([PathXML,FileXML]);
        end
        
    end
    
end

