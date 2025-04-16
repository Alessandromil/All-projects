classdef SEGMENT < handle
    %Segment parts of model
    %   
    
    properties
        Name              % Name of Segment                                           char
        Distals           % The distal of the Segment                                 [NDistal x 1]
                          % Distals(i).Child                                          handle Seg Child
                          % Distals(i).Joint                                          handle to the Joint Seg-Child
        Proximal          % The proximal of the Segment
                          % Proximal.Parent                                           handle Seg Parent
                          % Proximal.Joint                                            handle to the Joint Seg-Parent
        CoM               % Center of mass of Seg                                     LOCAL_POINT
        LocalVectors      % Vector of Vectors in Seg                                  LOCAL_VECTOR[]
        LocalVectorsAdded % Vector of Vectors to added constraints                    LOCAL_VECTOR[]
        LocalMarkers      % Vector of markers in Seg                                  LOCAL_POINT[]
        LocalPoints       % Vector of points in Seg                                   LOCAL_POINT[]
        LocalPointsRel    % LocalRelPoints% Vector of points rel in Seg                               LOCAL_POINT[]
        LocalALs          % Vector of anatomical landamarks in Seg                    LOCAL_POINT[]
        Fixed = 0;        % 0 free segment, 1 fixed segment,
        Root = 0;         % 0 free segment, 1 Root (More than one child)
        Glob_R_Seg        % The R matrix calculated from subject parameters files       double(3x3)
        Glob_Pos_OrSeg    % The position of the origing of the Seg calculated from subject parameter files double(3x1)
        Graphic           % Data of graphics is a cell(i,4); i= Number of graphics
                          % cell(i,1)= Name char/ cell(i,2)= Trans cell(3x1)/cell(i,3)= Rot cell(3x1)/ cell(i,4)= Scale cell(3x1)
        WireFrame         % struct with Wire frame definition
                          %   WireFrame.DrawSeq  Cell with names of point connected by the wire
                          %   WireFrame.Color    string with the color name for Compamm
                          %   WireFrame.Radius   string with the wire radius
        BasisType         % Type of Basis = 0 Basis formed by 3 vectors
                          %               = 1 Basis formed by 2 vectors and 2 points
        SegRef            % is a struct (1x1) that contain:
                            % Diri: strruc (1x1) that contain:
                            %      Global : char value (3x1)
        InvR              % is the inverse of [SegRef.Dir1.Local, SegRef.Dir2.Local, SegRef.Dir3.Local]   char{3,3}
        M                 % Mass Matrix double(12x12)
        I                 % Inercia moment Matrix in CoM??????                        double(3x3)
        Gender            % The gender of the subject(MALE or FEMALE)                                 char{}
        Length            % Length of the segment depending table 6.2 of D8           double
        Mass              % Mass of the segment                                       double
        TCS_Pos_Landmark  % The position of Landmark in TCS
        TCS_Pos_Markers   % The position of Markers in TCS
        TCS_Pos_TOALs     % The position of Landmark of the Trasfered origin segment in TCS
        TOLocalALs        % Vector of the AL of the Transfered origin segment         LOCAL_POINT[]
        Glob_R_LCS
        Glob_Pos_OrLCS
        PreSegMarkPos     % Global coord of the marker in previus segment             double[3xnMarkers]
        Glob_Pos_ParentOr % Glob_Pos_OrParent
        RamsisLCS_T_ALLCS % Transformation matrix between Ramsis and LCS              double[4x4]
        w                 % angular velocity of the segment                           double[3x1]
        wdot              % angular acceleration of the segment                       dounle[3x1]
        vdotOrBCS         % Glob_vdot_OrBCS % acceleration of the origin                                double[3x1]
        Seg_F_CoM         % The force in the CoM (Inertia Force)                      double[3x1]
        Seg_M_CoM         % The moment in the CoM (Inertia Moment)                    double[3x1]
        F_Ext             % External forces in the segment
                          % F_Ext.Sys                                                 double
                          %      0=LCS or 1=GlobalCS
                          % F_Ext.Value                                               double[NFrame x 3]
                          % F_Ext.Pos                                                 double[NFrame x 3]
        M_Ext             % Extenal moments in the segment                            double[NFrame x 3]
        SubAdditPar       % Subject additional Parameters not palpated                struct
                          % SubAdditPar.(ParameterName)
    end
    
    methods
        function S = SEGMENT(Name)
            % Constructor of Class
            S.Name = Name;
        end
        function addMarker(S,Marker)
            S.LocalMarkers =[S.LocalMarkers;LOCAL_POINT(Marker,S.Name)];
        end
        function addAL(S,AL)
            S.LocalALs =[S.LocalALs;LOCAL_POINT(AL,S.Name)];
        end
        function addDistal(S,Child,Joint)
            Distal.Child = Child;
            Distal.Joint = Joint;
            S.Distals = [S.Distals;Distal];
        end
        function addProximal(S,Parent,Joint)
            S.Proximal.Parent = Parent;
            S.Proximal.Joint = Joint;
        end
        function addPoint(S,Point)
            S.LocalPoints =[S.LocalPoints;LOCAL_POINT(Point,S.Name)];
        end
        function addRelPoint(S,PointRel)
            S.LocalPointsRel = [S.LocalPointsRel;LOCAL_POINT(PointRel,S.Name)];
        end
        function addVector(S,Vector)
            S.LocalVectors =[S.LocalVectors;LOCAL_VECTOR(Vector,S.Name)];
        end
        function addVector2Segment(S,Vector,LocCoord)
            LocalVector = LOCAL_VECTOR(Vector,S.Name);
            LocalVector.LocCoord = LocCoord;
            S.LocalVectorsAdded = [S.LocalVectorsAdded;LocalVector];
        end
        function calcRd(S,Glob_Pos_OrBody, Glob_Vec_x, Glob_Vec_y)
            % Glob_Pos_OrBody is an array (3x1)
            % Glob_Vec_x is an array (3x1)
            % Glob_Vec_y is an array (3x1)
            % Glob_Vec_z is an array (3x1)
            
            % calculate rotation matrix
            Glob_Vec_z = cross(Glob_Vec_x, Glob_Vec_y);
            Glob_Vec_y = cross(Glob_Vec_z, Glob_Vec_x);
            Glob_Vec_x = Glob_Vec_x / norm(Glob_Vec_x);
            Glob_Vec_y = Glob_Vec_y / norm(Glob_Vec_y);
            Glob_Vec_z = Glob_Vec_z / norm(Glob_Vec_z);
            Glob_R_Loc = [Glob_Vec_x, Glob_Vec_y, Glob_Vec_z];
            
            % build Body_T
            S.Glob_Pos_OrSeg = Glob_Pos_OrBody;
            S.Glob_R_Seg     = Glob_R_Loc;
        end
        function calcInertiaForces(S,InvDyn,q_i,q_2dot_i,Frame)
            NChilds = size(S.Distals,1);
            for i=1:NChilds
                ChildSeg = S.Distals(i).Child;
                Glob_R_Seg = S.getRd(q_i);
                Glob_R_Child  = ChildSeg.getRd(q_i);
                Child_R_Seg = Glob_R_Child'* Glob_R_Seg;
%                 Seg_R_Child = S.Distals(i).Joint.getR(Frame);
%                 Child_R_Seg = Seg_R_Child';
                Seg_w = S.w;
                Seg_wdot = S.wdot;
                [Child_Seg2Child_Anglesdot,Child_Seg2Child_Angles2dot]= S.calcSeg2ChildAngDers(S.Distals(i).Joint,Child_R_Seg,Frame);
%                 Seg2Child_Anglesdot = S.Distals(i).Joint.Sensor.Child_R_Joint * S.Distals(i).Joint.Anglesdot(Frame,:)';
%                 Seg2Child_Angles2dot = S.Distals(i).Joint.Sensor.Child_R_Joint * S.Distals(i).Joint.Angles2dot(Frame,:)';
                if ~strcmp(S.Distals(i).Joint.Type,'Float')
                    Seg_PtJnt2Child = S.Distals(i).Joint.Point1.LocCoord;
                    Child_PtJnt2Seg = S.Distals(i).Joint.Point2.LocCoord;
                else
                    PosInq = S.Distals(i).Joint.TrsSensor.EndPoint.PosInq;
                    Seg_PtJnt2Child = q_i(PosInq:PosInq+2)';
                    PointName = S.Distals(i).Joint.TrsSensor.EndPoint.Name;
                    PointIndex = getVecIndex(PointName,S.Distals(i).Joint.Sensor.Segment2.LocalPoints);
                    Child_PtJnt2Seg = S.Distals(i).Joint.Sensor.Segment2.LocalPoints(PointIndex).LocCoord;
                end
                Seg_vdotOrBCS = S.vdotOrBCS;
                if strcmp(S.Distals(i).Joint.Type,'Float')
                    PosInq = S.Distals(i).Joint.TrsSensor.EndPoint.PosInq;
                    Seg_vdotOrBCS = S.vdotOrBCS + Glob_R_Seg'*q_2dot_i(PosInq:PosInq + 2)';
                end
                
                Child_PtCoM = ChildSeg.CoM.LocCoord;
                Child_m = ChildSeg.Mass;
                ChildCoM_I = ChildSeg.I;% ChildSeg.Jeye(3);
                [Child_w,Child_wdot,Child_vdotOrBCS,Child_F_CoM,Child_M_CoM]= InvDyn.calcInertiaForces(Child_R_Seg,...
                    Seg_w,Seg_wdot,Child_Seg2Child_Anglesdot,Child_Seg2Child_Angles2dot,Seg_PtJnt2Child,...
                    Seg_vdotOrBCS,-Child_PtJnt2Seg,Child_PtCoM,Child_m,ChildCoM_I);
                ChildSeg.w = Child_w;
                ChildSeg.wdot = Child_wdot;
                ChildSeg.vdotOrBCS = Child_vdotOrBCS;
                ChildSeg.Seg_F_CoM = Child_F_CoM;
                ChildSeg.Seg_M_CoM = Child_M_CoM;
                if ~isempty(ChildSeg.Distals)
                    ChildSeg.calcInertiaForces(InvDyn,q_i,q_2dot_i,Frame);
                end
            end
            
        end
        function calcJointForces(S,InvDyn,q_i,Frame)
            % Comentario de que se calcula el joint que une al padre
            NFExt = size(S.F_Ext,2);
            NMExt = size(S.M_Ext,2);
            NChilds = size(S.Distals,1);
            Seg_F_Ext  = [0;0;0];
            Seg_Fd_Ext = [0;0;0];
            Seg_M_Ext  = [0;0;0];
            Seg_M_Jnt2Child = [0;0;0];
            Seg_F_Jnt2Child  = [0;0;0];
            Seg_Fd_Jnt2Child = [0;0;0];
            [Glob_R_Seg, Glob_Pos_OrSeg] = S.getRd(q_i);
            for i=1:NFExt
                if S.F_Ext(i).Sys == 0
                    Seg_F_Ext  = Seg_F_Ext + S.F_Ext(i).Value(Frame,:)';
                    Seg_Fd_Ext = Seg_Fd_Ext + cross(S.F_Ext(i).Pos(Frame,:),S.F_Ext(i).Value(Frame,:))';
                elseif S.F_Ext(i).Sys == 1
                    Seg_i_F_Ext = Glob_R_Seg'*S.F_Ext(i).Value(Frame,:)';
                    Seg_i_PosFExt =  Glob_R_Seg'* (S.F_Ext(i).Pos(Frame,:)'-Glob_Pos_OrSeg);
                    Seg_F_Ext = Seg_F_Ext + Seg_i_F_Ext;
                    Seg_Fd_Ext = Seg_Fd_Ext + cross(Seg_i_PosFExt,Seg_i_F_Ext);                    
                end
            end
            for i=1:NMExt
                Seg_M_Ext = Seg_M_Ext + S.M_Ext(i).Value(Frame);
            end
            for i=1:NChilds
                Glob_R_Child  = S.Distals(i).Child.getRd(q_i);
                ChildProxJointU = S.Distals(i).Child.Proximal.Joint.Sensor.Perm2x;
                ChildProxJointV = S.Distals(i).Child.Proximal.Joint.Sensor.Perm2y;
                ChildProxJointW = cross(ChildProxJointU,ChildProxJointV);
                ChildProxJoint_R_Child = [ChildProxJointU,ChildProxJointV,ChildProxJointW];
                Seg_R_ChildProxJoint = Glob_R_Seg'*Glob_R_Child*ChildProxJoint_R_Child';
                Seg_F_Jnt2Child  = Seg_F_Jnt2Child - (Seg_R_ChildProxJoint*S.Distals(i).Joint.F(Frame,:)');
                Seg_Pos_Jnt2Child = S.Distals(i).Joint.Point1.LocCoord;
                
                Seg_Fd_Jnt2Child = Seg_Fd_Jnt2Child - cross(Seg_Pos_Jnt2Child,(Seg_R_ChildProxJoint*S.Distals(i).Joint.F(Frame,:)'));
                Seg_M_Jnt2Child = Seg_M_Jnt2Child - (Seg_R_ChildProxJoint*S.Distals(i).Joint.M(Frame,:)');

%                 Antes                
%                 Seg_R_Child = Glob_R_Seg'*Glob_R_Child;
%                 Seg_F_Jnt2Child  = Seg_F_Jnt2Child - (Seg_R_Child*S.Distals(i).Joint.F(Frame,:)');
%                 Seg_Pos_Jnt2Child = S.Distals(i).Joint.Point1.LocCoord;
%                 
%                 Seg_Fd_Jnt2Child = Seg_Fd_Jnt2Child - cross(Seg_Pos_Jnt2Child,(Seg_R_Child*S.Distals(i).Joint.F(Frame,:)'));
%                 Seg_M_Jnt2Child = Seg_M_Jnt2Child - (Seg_R_Child*S.Distals(i).Joint.M(Frame,:)');
            end
            Seg_Fd_CoM = cross(S.CoM.LocCoord,S.Seg_F_CoM);
            if ~strcmp(S.Proximal.Joint.Type,'Float')
                Seg_Pos_Jnt2Parent = S.Proximal.Joint.Point2.LocCoord;
            else
                PointName = S.Proximal.Joint.TrsSensor.EndPoint.Name;
                PointIndex = getVecIndex(PointName,S.Proximal.Joint.Sensor.Segment2.LocalPoints);
                Seg_Pos_Jnt2Parent = S.Proximal.Joint.Sensor.Segment2.LocalPoints(PointIndex).LocCoord;
            end
            
%             JointIndex = getLocVecIndex(JointName,S.LocalPoints);
%             Seg_Pos_Jnt2Parent = S.LocalPoints(JointIndex).LocCoord;
            [Seg_F_Jnt2Parent,Seg_M_Jnt2Parent] = InvDyn.calcJointForces(S.Seg_F_CoM,Seg_F_Ext,Seg_F_Jnt2Child,...
                Seg_M_Jnt2Child,Seg_M_Ext,S.Seg_M_CoM,Seg_Fd_Jnt2Child,Seg_Fd_Ext,Seg_Fd_CoM,Seg_Pos_Jnt2Parent);
            ProxJointU = S.Proximal.Joint.Sensor.Perm2x;
            ProxJointV = S.Proximal.Joint.Sensor.Perm2y;
            PorxJointW = cross(ProxJointU,ProxJointV);
            ProxJoint_R_Seg = [ProxJointU,ProxJointV,PorxJointW];
            S.Proximal.Joint.F(Frame,:) = ProxJoint_R_Seg * Seg_F_Jnt2Parent;
            S.Proximal.Joint.M(Frame,:) = ProxJoint_R_Seg * Seg_M_Jnt2Parent;
            % Para tener las fuerzas en los joints en globales
            S.Proximal.Joint.FGlob(Frame,:) = Glob_R_Seg * Seg_F_Jnt2Parent;
            S.Proximal.Joint.MGlob(Frame,:) = Glob_R_Seg * Seg_M_Jnt2Parent;
            
            if (S.Proximal.Parent.Root == 0) && ~isempty(S.Proximal.Parent.Proximal) % Free segment && No Ground
                S.Proximal.Parent.calcJointForces(InvDyn,q_i,Frame);
            elseif (S.Proximal.Parent.Root == 1) % Root segment
                NChilds = size(S.Proximal.Parent.Distals,1);
                AllChild = 1;
                for i=1:NChilds
                    if isempty(S.Proximal.Parent.Distals(i).Joint.F)
                        AllChild = 0;
                    end
                end
                if AllChild == 1
                   S.Proximal.Parent.calcJointForces(InvDyn,q_i,Frame); 
                end
            end
        end
        function [Child_Seg2Child_Anglesdot,Child_Seg2Child_Angles2dot]=calcSeg2ChildAngDers(S,Joint,Child_R_Seg,Frame)
            % This function calc the value of velocities and acelerations between the segment and his child in the child BCS.
            
            % Calc velocities and acelerations in JCS of the segment
            SensorType = Joint.Sensor.Type;
            if strcmp(SensorType,'SPH')
                RotSeq = Joint.Sensor.RotSeq;
                Alpha1 = Joint.Angles(Frame,1);
                Alpha2 = Joint.Angles(Frame,2);
                Alpha1_dot = Joint.Anglesdot(Frame,1);
                Alpha2_dot = Joint.Anglesdot(Frame,2);
                Alpha3_dot = Joint.Anglesdot(Frame,3);
                Alpha1_2dot = Joint.Angles2dot(Frame,1);
                Alpha2_2dot = Joint.Angles2dot(Frame,2);
                Alpha3_2dot = Joint.Angles2dot(Frame,3);
                if strcmp(RotSeq,'123')
                    Axis1 = [1;0;0]; Axis2 = [0;1;0]; Axis3 = [0;0;1];
                    Rot1 = rot_x(Alpha1);
                    Rot2 = rot_y(Alpha2);
                    Rot1dot = [0,0,0;0,0,-Alpha1_dot;0,Alpha1_dot,0]*Rot1;
                    Rot2dot = [0,0,Alpha2_dot;0,0,0;-Alpha2_dot,0,0]*Rot2;   
                elseif strcmp(RotSeq,'321')
                    Axis1 = [0;0;1]; Axis2 = [0;1;0]; Axis3 = [1;0;0];
                    Rot1 = rot_z(Alpha1);
                    Rot2 = rot_y(Alpha2);
                    Rot1dot = [0,-Alpha1_dot,0;Alpha1_dot,0,0;0,0,0]*Rot1;
                    Rot2dot = [0,0,Alpha2_dot;0,0,0;-Alpha2_dot,0,0]*Rot2; 
                end
                JCSSeg_Seg2Child_Anglesdot  = Alpha1_dot*Axis1 + Alpha2_dot*Rot1*Axis2 + Alpha3_dot*Rot1*Rot2*Axis3;
                JCSSeg_Seg2Child_Angles2dot = Alpha1_2dot*Axis1 +  Alpha2_2dot*Rot1*Axis2 + Alpha2_dot*Rot1dot*Axis2 + ...
                                              Alpha3_2dot*Rot1*Rot2*Axis3 + Alpha3_dot*Rot1dot*Rot2*Axis3 + Alpha3_dot*Rot1*Rot2dot*Axis3; 
            else
                error(['Sensor of type ',SensorType,' does not exist. Only SPH sensors are valid'])
            end
            % Put velocities and acelerations in BCS of the segment (Fixed pemutation between BCS and JCS)
            BCS_R_JCS = Joint.Sensor.Parent_R_Joint;
            Seg_Seg2Child_Anglesdot  = BCS_R_JCS * JCSSeg_Seg2Child_Anglesdot;
            Seg_Seg2Child_Angles2dot = BCS_R_JCS * JCSSeg_Seg2Child_Angles2dot;
            % Put velocities and acelerations in BCS of the child segment (rotation between child and segment) 
            Child_Seg2Child_Anglesdot  = Child_R_Seg * Seg_Seg2Child_Anglesdot;
            Child_Seg2Child_Angles2dot = Child_R_Seg * Seg_Seg2Child_Angles2dot;
        end
        function calcALLCS_T_RamsisLCS(S,Glob_Vec_MeasVAxis)
            if strcmp (S.Name, 'Pelvis')
                NMarkers = size(S.LocalMarkers,1);
                Glob_Pos_MarkersStanding = [];
                AL_LCS_Pos_Markers = [];
                PostureIndex = getVecIndex('Standing',S.LocalMarkers(1).Point.Postures);
                for i=1:NMarkers
                    Glob_Pos_MarkersStanding = [Glob_Pos_MarkersStanding,S.LocalMarkers(i).Point.Postures(PostureIndex).Glob];
                    AL_LCS_Pos_Markers = [AL_LCS_Pos_Markers,S.LocalMarkers(i).LocCoord];
                end
                
                [Glob_R_AL_LCS,Or_AL_LCS] = calcOptPose(Glob_Pos_MarkersStanding,AL_LCS_Pos_Markers);
                RHJCIndex = getLocVecIndex('RHJC',S.LocalPoints);
                LHJCIndex = getLocVecIndex('LHJC',S.LocalPoints);
                LCS_Pos_RHJC = S.LocalPoints(RHJCIndex).LocCoord;
                LCS_Pos_LHJC = S.LocalPoints(LHJCIndex).LocCoord;
                LCS_Pos_midHJC = 0.5*(LCS_Pos_RHJC+LCS_Pos_LHJC);
                RamOr = LCS_Pos_midHJC;
                % Method using the PelvisTransversalRotation  Angle: Alpha. Rotation only in z
%                 Alpha = 15.55*pi/180;
%                 Alpha = acos(dot(Glob_Vec_MeasVAxis,Glob_R_AL_LCS(:,2))/(norm(Glob_Vec_MeasVAxis)*norm(Glob_R_AL_LCS(:,2))));
%                 PelX = [1;0;0];
%                 PelY = [0;1;0];
%                 PelZ = [0;0;1];
%                 RamB = -PelZ;
%                 Ram = [cos(Alpha),sin(Alpha),0; sin(Alpha), cos(Alpha) 0;0,0,1]*[PelY,PelX];
%                 RamT = Ram(:,1); RamN = Ram(:,2);
                % Method rotating the 3 angles
                Glob_Vec_RamT = Glob_Vec_MeasVAxis;
                Glob_Vec_RamN = cross(-Glob_R_AL_LCS(:,3),Glob_Vec_MeasVAxis);
                Glob_Vec_RamB = cross(Glob_Vec_RamT,Glob_Vec_RamN);
                Glob_R_Ram = [Glob_Vec_RamT,Glob_Vec_RamN,Glob_Vec_RamB];
                AL_LCS_R_Ramsis = Glob_R_AL_LCS'*Glob_R_Ram;
                RamT = AL_LCS_R_Ramsis(:,1);
                RamN = AL_LCS_R_Ramsis(:,2);
                RamB = AL_LCS_R_Ramsis(:,3);
                
            elseif strcmp (S.Name, 'LeftThigh') || strcmp (S.Name, 'RightThigh')
                if strcmp (S.Name, 'LeftThigh')
                    HJCIndex = getLocVecIndex('LHJC',S.LocalPoints);
                elseif  strcmp (S.Name, 'RightThigh')
                    HJCIndex = getLocVecIndex('RHJC',S.LocalPoints);
                end
                LCS_Pos_HJC = S.LocalPoints(HJCIndex).LocCoord;
                ThiX = [1;0;0];
                ThiY = [0;1;0];
                ThiZ = [0;0;1];
                RamB = ThiZ;
                RamT = -ThiY;
                RamN = ThiX;
                RamOr = LCS_Pos_HJC;
            elseif strcmp (S.Name, 'LeftShank') || strcmp (S.Name, 'RightShank')
                if strcmp (S.Name, 'LeftShank')
                    KJCIndex = getLocVecIndex('LKJC',S.LocalPoints);
                elseif strcmp (S.Name, 'RightShank')
                    KJCIndex = getLocVecIndex('RKJC',S.LocalPoints);
                end
                LCS_Pos_KJC = S.LocalPoints(KJCIndex).LocCoord;
                Alpha = 0*pi/180;
                ShaX = [1;0;0];
                ShaY = [0;1;0];
                ShaZ = [0;0;1];
                Ram = [cos(Alpha),0,sin(Alpha);0,1,0 ;sin(Alpha),0, cos(Alpha)]*[ShaX,ShaZ];
                RamT = -ShaY;
                RamB = Ram(:,2); RamN = Ram(:,1);
                RamOr = LCS_Pos_KJC;
            elseif strcmp (S.Name, 'LeftFoot') || strcmp (S.Name, 'RightFoot')
                if strcmp (S.Name, 'LeftFoot')
                    AJCIndex = getLocVecIndex('LAJC',S.LocalPoints);
                elseif strcmp (S.Name, 'RightFoot')
                    AJCIndex = getLocVecIndex('RAJC',S.LocalPoints);
                end
                LCS_Pos_AJC = S.LocalPoints(AJCIndex).LocCoord;
                FootX = [1;0;0];
                FootY = [0;1;0];
                FootZ = [0;0;1];
                RamT = FootX;
                RamB = FootZ; RamN = FootY;
%                 RamN = [0;0;1];% Vertical measurement
%                 RamT = [1;0;0];
%                 RamB = cross(RamT,RamN);
                RamOr = LCS_Pos_AJC;
            end
            LCS_R_Ramsis = [RamT,RamN,RamB];
            S.RamsisLCS_T_ALLCS = [LCS_R_Ramsis',-LCS_R_Ramsis'*RamOr;0,0,0,1];
        end
        function AJCIndex = calcAJCinFootLCS(S)
            if strcmp (S.Name, 'RightFoot')
                %                 FTCIndex = getLocVecIndex('RFTC',S.LocalALs);
                %                 FMEIndex = getLocVecIndex('RFME',S.LocalALs);
                %                 FLEIndex = getLocVecIndex('RFLE',S.LocalALs);
                AJCIndex = getLocVecIndex('RAJC',S.LocalPoints);
            elseif strcmp (S.Name, 'LeftFoot')||strcmp (S.Name, 'Foot')
                %                 FTCIndex = getLocVecIndex('LFTC',S.LocalALs);
                %                 FMEIndex = getLocVecIndex('LFME',S.LocalALs);
                %                 FLEIndex = getLocVecIndex('LFLE',S.LocalALs);
                AJCIndex = getLocVecIndex('LAJC',S.LocalPoints);
            end
            S.LocalPoints(AJCIndex).LocCoord = [0;0;0];
        end
        function [RHJCIndex,LHJCIndex]= calcHJCinPelvis(S)
            LIASIndex = getLocVecIndex('LIAS',S.LocalALs);
            RIASIndex = getLocVecIndex('RIAS',S.LocalALs);
            LIPSIndex = getLocVecIndex('LIPS',S.LocalALs);
            RIPSIndex = getLocVecIndex('RIPS',S.LocalALs);
            RHJCIndex = getLocVecIndex('RHJC',S.LocalPoints);
            LHJCIndex = getLocVecIndex('LHJC',S.LocalPoints);
            LIAS = S.LocalALs(LIASIndex).LocCoord';
            RIAS = S.LocalALs(RIASIndex).LocCoord';
            LIPS = S.LocalALs(LIPSIndex).LocCoord';
            RIPS = S.LocalALs(RIPSIndex).LocCoord';
            % To calculate Right HJC 
            S032=[LIAS;LIPS;RIAS;RIPS];
            % The second value must be 1 the select method and second 0 because is right
            [Glb_Pos_Head, Loc_Pos_Head, Rad_Head, Rot_Lc_2_Gl_CS_Ori] = S.calcHJCsInPelvisLCS_ULB(S032, 1, 0);
            S.LocalPoints(RHJCIndex).LocCoord = Glb_Pos_Head;
            % To calculate Left HJC
            S032=[RIAS;RIPS;LIAS;LIPS];
            % The second value must be 1 the select method and second 1 because is left
            [Glb_Pos_Head, Loc_Pos_Head, Rad_Head, Rot_Lc_2_Gl_CS_Ori] = S.calcHJCsInPelvisLCS_ULB(S032, 1, 1);
            S.LocalPoints(LHJCIndex).LocCoord = Glb_Pos_Head;
        end
        function LJCIndex = calcLJCinPelvis(S)
            % get the values of ALs in LCS
            LIASIndex = getLocVecIndex('LIAS',S.LocalALs);
            RIASIndex = getLocVecIndex('RIAS',S.LocalALs);
            LIPSIndex = getLocVecIndex('LIPS',S.LocalALs);
            RIPSIndex = getLocVecIndex('RIPS',S.LocalALs);
            RHJCIndex = getLocVecIndex('RHJC',S.LocalPoints);
            LHJCIndex = getLocVecIndex('LHJC',S.LocalPoints);
            LJCIndex = getLocVecIndex('MLJC',S.LocalPoints);
            LCS_Pos_LIAS = S.LocalALs(LIASIndex).LocCoord;
            LCS_Pos_RIAS = S.LocalALs(RIASIndex).LocCoord;
            LCS_Pos_LIPS = S.LocalALs(LIPSIndex).LocCoord;
            LCS_Pos_RIPS = S.LocalALs(RIPSIndex).LocCoord;
            LCS_Pos_RHJC = S.LocalPoints(RHJCIndex).LocCoord;
            LCS_Pos_LHJC = S.LocalPoints(LHJCIndex).LocCoord;

            % get the values for a average pelvis of ALs in Dumas CS
            if strcmpi (S.Gender,'male')
                DuCs_Pos_RIAS = [78;7;112];
				DuCs_Pos_LIAS = [78;7;-112];
				DuCs_Pos_midIPS = [-102;7;0];
				DuCs_Pos_RHJC = [56;-75;81];	 
				DuCs_Pos_LHJC = [56;-75;-81];
            elseif strcmpi (S.Gender,'female')
				DuCs_Pos_RIAS = [87;-13;119];
				DuCs_Pos_LIAS = [87;-13;-119];
				DuCs_Pos_midIPS = [-108;13;0];
				DuCs_Pos_RHJC = [54;-93;88];	 
				DuCs_Pos_LHJC = [54;-93;-88];
            end
            % Calculate the scalate matrix
            DuCs_Pos_midIAS = 0.5*(DuCs_Pos_RIAS + DuCs_Pos_LIAS);
            DuCs_Pos_midHJC = 0.5*(DuCs_Pos_RHJC + DuCs_Pos_LHJC);
            DuCs_Vec_midIASmidIPS = DuCs_Pos_midIPS - DuCs_Pos_midIAS;
            DuCs_Vec_midIASmidHJC = DuCs_Pos_midHJC - DuCs_Pos_midIAS;
            
            LCS_Pos_midIAS = 0.5*(LCS_Pos_RIAS+LCS_Pos_LIAS);
			LCS_Pos_midHJC = 0.5*(LCS_Pos_RHJC+LCS_Pos_LHJC);
			LCS_Pos_midIPS = 0.5*(LCS_Pos_RIPS+LCS_Pos_LIPS);
			LCS_Vec_midIASmidIPS = LCS_Pos_midIPS - LCS_Pos_midIAS;
			LCS_Vec_midIASmidHJC = LCS_Pos_midHJC - LCS_Pos_midIAS;
            
            Sx = abs(LCS_Vec_midIASmidIPS(1)/DuCs_Vec_midIASmidIPS(1));
            Sy = abs(LCS_Vec_midIASmidHJC(2)/DuCs_Vec_midIASmidHJC(2));
            Sz = 0;
            Scal = [Sx 0 0; 0 Sy 0; 0 0 Sz];
            % The LCS and DumasCS are equal but with different centre
            % The LCS centre is in midIAS and DumasCS in LJC
            LCS_Pos_AvgLJC = -DuCs_Pos_midIAS; 
            LCS_Pos_LJC = Scal * LCS_Pos_AvgLJC;

            S.LocalPoints(LJCIndex).LocCoord = LCS_Pos_LJC;
        end
        function AJCIndex = calcAJCinShank(S)
            if strcmp (S.Name, 'RightShank')
                %                 FTCIndex = getLocVecIndex('RFTC',S.LocalALs);
                %                 FMEIndex = getLocVecIndex('RFME',S.LocalALs);
                %                 FLEIndex = getLocVecIndex('RFLE',S.LocalALs);
                AJCIndex = getLocVecIndex('RAJC',S.LocalPoints);
            elseif strcmp (S.Name, 'LeftShank')||strcmp (S.Name, 'Shank')
                %                 FTCIndex = getLocVecIndex('LFTC',S.LocalALs);
                %                 FMEIndex = getLocVecIndex('LFME',S.LocalALs);
                %                 FLEIndex = getLocVecIndex('LFLE',S.LocalALs);
                AJCIndex = getLocVecIndex('LAJC',S.LocalPoints);
            end
            S.LocalPoints(AJCIndex).LocCoord = [0;0;0];
        end
        function KJCIndex = calcKJCinShank(S,TCS_Pos_OrThigh,TCS_Pos_OrShank,LCS_R_TCS)
            if strcmp (S.Name, 'RightShank')
                %                 FTCIndex = getLocVecIndex('RFTC',S.LocalALs);
                %                 FMEIndex = getLocVecIndex('RFME',S.LocalALs);
                %                 FLEIndex = getLocVecIndex('RFLE',S.LocalALs);
                KJCIndex = getLocVecIndex('RKJC',S.LocalPoints);
            elseif strcmp (S.Name, 'LeftShank')||strcmp (S.Name, 'Shank')
                %                 FTCIndex = getLocVecIndex('LFTC',S.LocalALs);
                %                 FMEIndex = getLocVecIndex('LFME',S.LocalALs);
                %                 FLEIndex = getLocVecIndex('LFLE',S.LocalALs);
                KJCIndex = getLocVecIndex('LKJC',S.LocalPoints);
            end
%             LCS_Pos_OrThigh = changeCoordSys(Glob_Pos_OrThigh,S.Glob_Pos_OrLCS,S.Glob_R_LCS');
            LCS_Pos_OrThigh = changeCoordSys(TCS_Pos_OrThigh,TCS_Pos_OrShank,LCS_R_TCS);
            S.LocalPoints(KJCIndex).LocCoord = LCS_Pos_OrThigh;
        end
        function [TCS_Pos_HJC,HJCIndex] = calcHJCinThigh(S)
            if strcmp (S.Name, 'RightThigh')
                FTCIndex = getLocVecIndex('RFTC',S.LocalALs);
                FMEIndex = getLocVecIndex('RFME',S.LocalALs);
                FLEIndex = getLocVecIndex('RFLE',S.LocalALs);
                HJCIndex = getLocVecIndex('RHJC',S.LocalPoints);
                RL = 0; % 0 In the Right
            elseif strcmp (S.Name, 'LeftThigh')||strcmp (S.Name, 'Thigh')
                FTCIndex = getLocVecIndex('LFTC',S.LocalALs);
                FMEIndex = getLocVecIndex('LFME',S.LocalALs);
                FLEIndex = getLocVecIndex('LFLE',S.LocalALs);
                HJCIndex = getLocVecIndex('LHJC',S.LocalPoints);
                RL = 1; % 1 In the Left
            end
%             FTC = S.LocalALs(FTCIndex).LocCoord';
%             FME = S.LocalALs(FMEIndex).LocCoord';
%             FLE = S.LocalALs(FLEIndex).LocCoord';
            FTC = S.TCS_Pos_Landmark(:,FTCIndex)';
            FME = S.TCS_Pos_Landmark(:,FMEIndex)';
            FLE = S.TCS_Pos_Landmark(:,FLEIndex)';
            AL = [FTC;FME;FLE];
            [Glb_Pos_Head, Loc_Pos_Head, Rad_Head, Rot_Lc_2_Gl_CS_Ori] = S.calcKJCinThighLCS_ULB(AL,RL);
            TCS_Pos_HJC = Glb_Pos_Head(1,:)';
            
        end
        function [TCS_Pos_HJC,HJCIndex] = calcHJCinThighWithPelvisALs(S)
            % get the Index of Pelvis ALs 
            LIASIndex = getLocVecIndex('LIAS',S.TOLocalALs);
            RIASIndex = getLocVecIndex('RIAS',S.TOLocalALs);
            LIPSIndex = getLocVecIndex('LIPS',S.TOLocalALs);
            RIPSIndex = getLocVecIndex('RIPS',S.TOLocalALs);
            % get the values of Pelvis ALs in Thigh TCS
            LIAS = S.TCS_Pos_TOALs(:,LIASIndex)';
            RIAS = S.TCS_Pos_TOALs(:,RIASIndex)';
            LIPS = S.TCS_Pos_TOALs(:,LIPSIndex)';
            RIPS = S.TCS_Pos_TOALs(:,RIPSIndex)';
            if strcmpi(S.Name,'RightThigh')
                HJCIndex = getLocVecIndex('RHJC',S.LocalPoints);
                % To calculate Right HJC
                S032=[LIAS;LIPS;RIAS;RIPS];
                % The second value must be 1 the select method and second 0 because is right
                [Glb_Pos_Head, Loc_Pos_Head, Rad_Head, Rot_Lc_2_Gl_CS_Ori] = S.calcHJCsInPelvisLCS_ULB(S032, 1, 0);
                TCS_Pos_HJC = Glb_Pos_Head;
            elseif strcmpi(S.Name,'LeftThigh')
                HJCIndex = getLocVecIndex('LHJC',S.LocalPoints);
                % To calculate Left HJC
                S032=[RIAS;RIPS;LIAS;LIPS];
                % The second value must be 1 the select method and second 1 because is left
                [Glb_Pos_Head, Loc_Pos_Head, Rad_Head, Rot_Lc_2_Gl_CS_Ori] = S.calcHJCsInPelvisLCS_ULB(S032, 1, 1);
                TCS_Pos_HJC = Glb_Pos_Head;
            end
        end
        function KJCIndex = calcKJCinThigh(S)
            if strcmp (S.Name, 'RightThigh')
                %                 FTCIndex = getLocVecIndex('RFTC',S.LocalALs);
                %                 FMEIndex = getLocVecIndex('RFME',S.LocalALs);
                %                 FLEIndex = getLocVecIndex('RFLE',S.LocalALs);
                KJCIndex = getLocVecIndex('RKJC',S.LocalPoints);
            elseif strcmp (S.Name, 'LeftThigh')||strcmp (S.Name, 'Thigh')
                %                 FTCIndex = getLocVecIndex('LFTC',S.LocalALs);
                %                 FMEIndex = getLocVecIndex('LFME',S.LocalALs);
                %                 FLEIndex = getLocVecIndex('LFLE',S.LocalALs);
                KJCIndex = getLocVecIndex('LKJC',S.LocalPoints);
            end
            S.LocalPoints(KJCIndex).LocCoord = [0;0;0];
        end
        function calcGlobWithOptPose(S,Glob_R_Seg,Glob_Pos_Or)
            % CALCGLOBWITHOPTPOSE calcs the value of the Global coordinates of Vectors, Points and Markers with
            % optimal posture obtained by motion file
            NPoints  = size(S.LocalPoints,1);
            NVectors = size(S.LocalVectors,1);
            NMarkers = size(S.LocalMarkers,1);
            for i=1:NVectors
                S.LocalVectors(i).Vector.GlobalCoord = Glob_R_Seg*S.LocalVectors(i).LocCoord;
            end
            for i=1:NPoints
                S.LocalPoints(i).Point.GlobalCoord = Glob_Pos_Or + Glob_R_Seg*S.LocalPoints(i).LocCoord;
                % The propertie of the point OptPose is uptdated for future checking
                S.LocalPoints(i).Point.OptPose = 1;
            end
            for i=1:NMarkers
                S.LocalMarkers(i).Point.GlobalCoord = Glob_Pos_Or + Glob_R_Seg*S.LocalMarkers(i).LocCoord;
            end
        end
        function q_t = calcq_tWithOptPose(S,Glob_R_Seg,Glob_Pos_Or,Frame,q_t)
            NPoints  = size(S.LocalPoints,1);
            NVectors = size(S.LocalVectors,1);
            NMarkers = size(S.LocalMarkers,1);
            for i=1:NVectors
                if S.LocalVectors(i).Vector.Fixed == 0
                    PosInq = S.LocalVectors(i).Vector.PosInq;
                    q_t(Frame,PosInq:PosInq+2) = Glob_R_Seg*S.LocalVectors(i).LocCoord;
                end
            end
            for i=1:NPoints
                if S.LocalPoints(i).Point.Fixed == 0
                    PosInq = S.LocalPoints(i).Point.PosInq;
                    q_t(Frame,PosInq:PosInq+2) = Glob_Pos_Or + Glob_R_Seg*S.LocalPoints(i).LocCoord;
                end
            end
            for i=1:NMarkers
                PosInq = S.LocalMarkers(i).Point.PosInq;
                q_t(Frame,PosInq:PosInq+2) = Glob_Pos_Or + Glob_R_Seg*S.LocalMarkers(i).LocCoord;
            end
        end
        function calcLCS(S)
            if strcmp(S.Name, 'Pelvis')
                [TCS_Pos_OrPelvis,LCS_R_TCS] = S.calcPelvisLCS;
                S.setMarkersAndALsInLCS(TCS_Pos_OrPelvis,LCS_R_TCS);
                [RHJCIndex,LHJCIndex]= S.calcHJCinPelvis();
                LJCIndex = S.calcLJCinPelvis();
                Length = S.LocalPoints(LJCIndex).LocCoord - ((S.LocalPoints(RHJCIndex).LocCoord+S.LocalPoints(LHJCIndex).LocCoord)/2);
                S.Length = norm(Length);
                              
            elseif strcmp (S.Name, 'RightThigh')||strcmp (S.Name, 'LeftThigh')||strcmp (S.Name, 'Thigh')
%                 [TCS_Pos_HJC,HJCIndex] = S.calcHJCinThigh();
                [TCS_Pos_HJC,HJCIndex] = S.calcHJCinThighWithPelvisALs();
                KJCIndex = S.calcKJCinThigh();
                [TCS_Pos_OrThigh,LCS_R_TCS] = S.calcThighLCS(TCS_Pos_HJC);
                S.setMarkersAndALsInLCS(TCS_Pos_OrThigh,LCS_R_TCS);
                LCS_Pos_HJC  = changeCoordSys(TCS_Pos_HJC,TCS_Pos_OrThigh,LCS_R_TCS);
                S.LocalPoints(HJCIndex).LocCoord = LCS_Pos_HJC;
                Length = S.LocalPoints(HJCIndex).LocCoord - S.LocalPoints(KJCIndex).LocCoord;
                S.Length = norm(Length);
                                
            elseif strcmp (S.Name, 'RightShank')||strcmp (S.Name, 'LeftShank')||strcmp (S.Name, 'Shank')
                [TCS_Pos_OrThigh,TCS_Pos_OrShank,LCS_R_TCS] = S.calcShankLCS();
                S.setMarkersAndALsInLCS(TCS_Pos_OrShank,LCS_R_TCS);
                KJCIndex = S.calcKJCinShank(TCS_Pos_OrThigh,TCS_Pos_OrShank,LCS_R_TCS);
                AJCIndex = S.calcAJCinShank();
                Length = S.LocalPoints(KJCIndex).LocCoord - S.LocalPoints(AJCIndex).LocCoord;
                S.Length = norm(Length);
            elseif strcmp (S.Name, 'RightFoot')||strcmp (S.Name, 'LeftFoot')||strcmp (S.Name, 'Foot')
                [TCS_Pos_OrFoot,LCS_R_TCS] = S.calcFootLCS();
                if strcmp (S.Name, 'RightFoot')
                    FM1Index = getLocVecIndex('RFM1',S.LocalALs);
                    FM5Index = getLocVecIndex('RFM5',S.LocalALs);
                elseif strcmp (S.Name, 'LeftFoot')||strcmp (S.Name, 'Foot')
                    FM1Index = getLocVecIndex('LFM1',S.LocalALs);
                    FM5Index = getLocVecIndex('LFM5',S.LocalALs);
                end
                S.setMarkersAndALsInLCS(TCS_Pos_OrFoot,LCS_R_TCS);
                AJCIndex = S.calcAJCinFootLCS(); 
                Length = S.LocalPoints(AJCIndex).LocCoord - ((S.LocalALs(FM1Index).LocCoord + S.LocalALs(FM5Index).LocCoord)/2);
                S.Length = norm(Length);
            elseif strcmp (S.Name,'Thorax')
                [TCS_Pos_OrThorax,LCS_R_TCS] = S.calcThoraxLCS();
                S.setMarkersAndALsInLCS(TCS_Pos_OrThorax,LCS_R_TCS);
            elseif strcmp (S.Name,'RightClavicle') || strcmp (S.Name,'LeftClavicle')
                [TCS_Pos_OrClavicle,LCS_R_TCS] = S.calcClavicleLCS();
                S.setMarkersAndALsInLCS(TCS_Pos_OrClavicle,LCS_R_TCS);
            elseif strcmp (S.Name,'RightScapula') || strcmp (S.Name,'LeftScapula')
                [TCS_Pos_OrScapula,LCS_R_TCS] = S.calcScapulaLCS();
                S.setMarkersAndALsInLCS(TCS_Pos_OrScapula,LCS_R_TCS);
            elseif strcmp (S.Name,'RightHumerus') || strcmp (S.Name,'LeftHumerus')
                S.calcGHJCinHumerus;
                [TCS_Pos_OrHumerus,LCS_R_TCS] = S.calcHumerusLCS();
                S.setMarkersAndALsInLCS(TCS_Pos_OrHumerus,LCS_R_TCS);
            elseif strcmp (S.Name,'RightForearm') || strcmp (S.Name,'LeftForearm')
                [TCS_Pos_OrForearm,LCS_R_TCS] = S.calcForearmLCS();
                S.setMarkersAndALsInLCS(TCS_Pos_OrForearm,LCS_R_TCS);
            elseif strcmp (S.Name,'RightUlna') || strcmp (S.Name,'LeftUlna')
                [TCS_Pos_OrUlna,LCS_R_TCS] = S.calcUlnaLCS();
                S.setMarkersAndALsInLCS(TCS_Pos_OrUlna,LCS_R_TCS);
            elseif strcmp (S.Name,'RightRadius') || strcmp (S.Name,'LeftRadius')
                [TCS_Pos_OrRadius,LCS_R_TCS] = S.calcRadiusLCS();
                S.setMarkersAndALsInLCS(TCS_Pos_OrRadius,LCS_R_TCS);
            elseif strcmp (S.Name,'RightHand') || strcmp (S.Name,'LeftHand')
                [TCS_Pos_OrHand,LCS_R_TCS] = S.calcHandLCS();
                S.setMarkersAndALsInLCS(TCS_Pos_OrHand,LCS_R_TCS);    
            else
                error('%s is not a valid name for the segment.',S.Name);
            end

        end
        function [TCS_Pos_OrLCS,LCS_R_TCS] = calcClavicleLCS(S)
            % Check if is Right or Left Clavicle
            if strcmpi(S.Name,'RightClavicle')
                PreFix = 'R';
            elseif strcmpi(S.Name,'LeftClavicle')
                PreFix = 'L';
            end
            % get ALs Index
            CSJIndex = getLocVecIndex([PreFix,'CSJ'],S.LocalALs);
            CAJIndex = getLocVecIndex([PreFix,'CAJ'],S.LocalALs);
            SJNIndex = getLocVecIndex('SJN',S.TOLocalALs);
            SXSIndex = getLocVecIndex('SXS',S.TOLocalALs);
            TV3Index = getLocVecIndex('TV3',S.TOLocalALs);
            TV8Index = getLocVecIndex('TV8',S.TOLocalALs);
            % check if ALs have been measured
            if CSJIndex==0, error(['Clavicle landmark ',PreFix,'CSJ not found in subjec parameter file']); end
            if CAJIndex==0, error(['Clavicle landmark ',PreFix,'CAJ not found in subjec parameter file']); end
            % get ALs coordinates in TCS
            CSJ = S.TCS_Pos_Landmark(:,CSJIndex);
            CAJ = S.TCS_Pos_Landmark(:,CAJIndex);
            SJN = S.TCS_Pos_TOALs(:,SJNIndex);
            SXS = S.TCS_Pos_TOALs(:,SXSIndex);
            TV3 = S.TCS_Pos_TOALs(:,TV3Index);
            TV8 = S.TCS_Pos_TOALs(:,TV8Index);
            YTho = (((SJN+TV3)/2)-((SXS-TV8)/2));
            % APPLY LCS DEFINITION to TCS coordinates
            % Originin CSJ
            TCS_Pos_OrLCS = CSJ;
            % The axes of LCS
            TCS_R_LCS(:,3) = (CAJ-CSJ)/norm(CAJ-CSJ);
            if strcmp(PreFix,'R')
                TCS_R_LCS(:,1) = cross(YTho,TCS_R_LCS(:,3))/norm(cross(YTho,TCS_R_LCS(:,3)));
            elseif strcmp(PreFix,'L')
                TCS_R_LCS(:,1) = cross(TCS_R_LCS(:,3),YTho)/norm(cross(TCS_R_LCS(:,3),YTho));
            end
            TCS_R_LCS(:,2) = cross(TCS_R_LCS(:,3),TCS_R_LCS(:,1));
            LCS_R_TCS = TCS_R_LCS';
        end
        function [TCS_Pos_OrLCS,LCS_R_TCS] = calcFootLCS(S)
            
            if strcmpi(S.Name,'RightFoot')
                PreFix = 'R';
            elseif strcmpi(S.Name,'LeftFoot')
                PreFix = 'L';
            end
            FCCIndex = getLocVecIndex([PreFix,'FCC'],S.LocalALs);
            FM1Index = getLocVecIndex([PreFix,'FM1'],S.LocalALs);
            FM5Index = getLocVecIndex([PreFix,'FM5'],S.LocalALs);
            FM2Index = getLocVecIndex([PreFix,'FM2'],S.LocalALs);
            TAMIndex = getLocVecIndex([PreFix,'TAM'],S.TOLocalALs);
            FALIndex = getLocVecIndex([PreFix,'FAL'],S.TOLocalALs);
            
            if FCCIndex==0, error(['Foot landmark ',PreFix,'FCC not found in subjec parameter file']); end
            if FM1Index==0, error(['Foot landmark ',PreFix,'FM1 not found in subjec parameter file']); end
            if FM5Index==0, error(['Foot landmark ',PreFix,'FM5 not found in subjec parameter file']); end
            if FM2Index==0
                disp(['    WARNING: Foot landmark ',PreFix,'FM2 not found in subjec parameter file']);
                disp( '             TUM Foot LCS definition will be used based on FCC, FM1 & FM5');
            end
                  
            % get AL values in TCS
            FCC = S.TCS_Pos_Landmark(:,FCCIndex);
            FM1 = S.TCS_Pos_Landmark(:,FM1Index);
            FM5 = S.TCS_Pos_Landmark(:,FM5Index);
            TAM = S.TCS_Pos_TOALs(:,TAMIndex); 
            FAL = S.TCS_Pos_TOALs(:,FALIndex);
            TCS_Pos_OrShank = (TAM+FAL)/2;
            
            if FM2Index ~= 0 % Foot LCS definition based on 4 ALs: RFCC, RFM1, RFM5 & RFM2
                
                FM2 = S.TCS_Pos_Landmark(:,FM2Index);
                % Origin of foot: origin of Shank
                TCS_Pos_OrLCS = TCS_Pos_OrShank;
                % The axes of LCS
                if strcmp(PreFix,'L') % Left
                    TCS_R_LCS(:,2) = cross(FM1-FCC,FM5-FCC)/norm(cross(FM1-FCC,FM5-FCC));
                elseif  strcmp(PreFix,'R')% Right
                    TCS_R_LCS(:,2) = cross(FM5-FCC,FM1-FCC)/norm(cross(FM5-FCC,FM1-FCC));
                end
                C = dot(TCS_R_LCS(:,2),(FM2-FCC))*TCS_R_LCS(:,2);
                TCS_R_LCS(:,1) = ((FM2-FCC)-C)/norm(((FM2-FCC)-C));
                TCS_R_LCS(:,3) = cross(TCS_R_LCS(:,1),TCS_R_LCS(:,2));
                LCS_R_TCS = TCS_R_LCS';
                % S.drawLCSinTCS(TCS_Pos_OrLCS,TCS_R_LCS);
                
            else % FM2Index == 0 % Foot LCS definition based on 3 ALs: FCC, FM1 & FM5 (TUM definition)
                % Origin of foot: origin of Shank
                TCS_Pos_OrLCS = TCS_Pos_OrShank;
                % The axes of LCS
                if strcmp(PreFix,'L') % Left
                    % Y axis
                    TCS_R_LCS(:,2) = cross(FM1-FCC,FM5-FCC)/norm(cross(FM1-FCC,FM5-FCC));
                    if isfield(S.SubAdditPar,'Dist_LFM1_LFM2')
                        Dist_FM1_FM2 = S.SubAdditPar.Dist_LFM1_LFM2;
                    else
                        error('Parameter Dist_LFM1_LFM2 does not exist in Subject interface file')
                    end
                    
                elseif strcmp(PreFix,'R') % Right
                    % Y axis
                    TCS_R_LCS(:,2) = cross(FM5-FCC,FM1-FCC)/norm(cross(FM5-FCC,FM1-FCC));
                    if isfield(S.SubAdditPar,'Dist_RFM1_RFM2')
                        Dist_FM1_FM2 = S.SubAdditPar.Dist_RFM1_RFM2;
                    else
                        error('Parameter Dist_RFM1_RFM2 does not exist in Subject interface file')
                    end
                end
                % X axis: line conecting FCC and a point between FM1 and FM5
                uVec_FM1_FM5 = (FM5-FM1)/norm(FM5-FM1);
                FM2 = FM1 + Dist_FM1_FM2 * uVec_FM1_FM5;
                TCS_R_LCS(:,1) = (FM2-FCC)/norm(FM2-FCC);                
                TCS_R_LCS(:,3) = cross(TCS_R_LCS(:,1),TCS_R_LCS(:,2));
                LCS_R_TCS = TCS_R_LCS';
%                 S.drawLCSinTCS(TCS_Pos_OrLCS,TCS_R_LCS);
                
            end
        end
        function [TCS_Pos_OrLCS,LCS_R_TCS] = calcForearmLCS(S)
            % Check if is Right or Left Clavicle
            if strcmpi(S.Name,'RightForearm')
                PreFix = 'R';
            elseif strcmpi(S.Name,'LeftForearm')
                PreFix = 'L';
            end
            % get ALs Index
            USPIndex = getLocVecIndex([PreFix,'USP'],S.LocalALs);
            RSPIndex = getLocVecIndex([PreFix,'RSP'],S.LocalALs);
            % check if ALs have been measured
            if USPIndex==0, error(['Forearm landmark ',PreFix,'USP not found in subjec parameter file']); end
            if RSPIndex==0, error(['Forearm landmark ',PreFix,'RSP not found in subjec parameter file']); end
            % get ALs coordinates in TCS
            RSP = S.TCS_Pos_Landmark(:,RSPIndex);
            USP = S.TCS_Pos_Landmark(:,USPIndex);
            HME = S.TCS_Pos_TransLandmark(:,HMEIndex);
            HLE = S.TCS_Pos_TransLandmark(:,HLEIndex);
            % APPLY LCS DEFINITION to TCS coordinates
            % Originin USP
            TCS_Pos_OrLCS = USP;
            % The axes of LCS
            TCS_R_LCS(:,2) = ((HME+HLE/2)-USP)/norm(((HME+HLE/2)-USP));
            if strcmp(PreFix,'R')
                TCS_R_LCS(:,1) = cross(((HME+HLE/2)-USP),RSP-USP)/norm(cross(((HME+HLE/2)-USP),RSP-USP));
            elseif strcmp(PreFix,'L')
                TCS_R_LCS(:,1) = cross(RSP-USP,((HME+HLE/2)-USP))/norm(cross(RSP-USP,((HME+HLE/2)-USP)));
            end
            TCS_R_LCS(:,3) = cross(TCS_R_LCS(:,1),TCS_R_LCS(:,2));
            LCS_R_TCS = TCS_R_LCS';
        end
        function [TCS_Pos_OrLCS,LCS_R_TCS] = calcHandLCS(S)
             % Check if is Right or Left Clavicle
            if strcmpi(S.Name,'RightHand')
                PreFix = 'R';
            elseif strcmpi(S.Name,'LeftHand')
                PreFix = 'L';
            end
            % get ALs Index
            HM2Index = getLocVecIndex([PreFix,'HM2'],S.LocalALs);
            HM5Index = getLocVecIndex([PreFix,'HM5'],S.LocalALs);
            RSPIndex = getLocVecIndex([PreFix,'RSP'],S.LocalALs);
            % check if ALs have been measured
            if HM2Index==0, error(['Hand landmark ',PreFix,'HM2 not found in subjec parameter file']); end
            if HM5Index==0, error(['Hand landmark ',PreFix,'HM5 not found in subjec parameter file']); end
            % get ALs coordinates in TCS
            HM2 = S.TCS_Pos_Landmark(:,HM2Index);
            HM5 = S.TCS_Pos_Landmark(:,HM5Index);
            RSP = S.TCS_Pos_TOALs(:,RSPIndex);
            % APPLY LCS DEFINITION to TCS coordinates
            % Originin USP
            TCS_Pos_OrLCS = RSP;
            % The axes of LCS
            TCS_R_LCS(:,3) = (HM5+HM2)/norm(HM5+HM2);
            if strcmp(PreFix,'R')
                TCS_R_LCS(:,1) = cross(HM2-RSP,HM5-RSP)/norm(cross(HM2-RSP,HM5-RSP));
            elseif strcmp(PreFix,'L')
                TCS_R_LCS(:,1) = cross(HM5-RSP,HM2-RSP)/norm(cross(HM5-RSP,HM2-RSP));
            end
            TCS_R_LCS(:,2) = cross(TCS_R_LCS(:,3),TCS_R_LCS(:,1));
            LCS_R_TCS = TCS_R_LCS';
            
        end
        function [TCS_Pos_OrLCS,LCS_R_TCS] = calcHumerusLCS(S)
            % Check if is Right or Left Humerus
            if strcmpi(S.Name,'RightHumerus')
                PreFix = 'R';
            elseif strcmpi(S.Name,'LeftHumerus')
                PreFix = 'L';
            end
            % get ALs Index
            HMEIndex = getLocVecIndex([PreFix,'HME'],S.LocalALs);
            HLEIndex = getLocVecIndex([PreFix,'HLE'],S.LocalALs);
            USPIndex = getLocVecIndex([PreFix,'USP'],S.TOLocalALs);
            GHJCIndex = getLocVecIndex([PreFix,'GHJC'],S.LocalPoints);
            % check if ALs have been measured
            if HMEIndex==0, error(['Humerus landmark ',PreFix,'HME not found in subjec parameter file']); end
            if HLEIndex==0, error(['Humerus landmark ',PreFix,'HLE not found in subjec parameter file']); end
            % get ALs coordinates in TCS
            HME = S.TCS_Pos_Landmark(:,HMEIndex);
            HLE = S.TCS_Pos_Landmark(:,HLEIndex);
            GHJC = S.TCS_Pos_Point(:,GHJCIndex);
            USP = S.TCS_Pos_TOALs(:,USPIndex);
            YFa = ((HME+HLE/2)-USP);
            % APPLY LCS DEFINITION to TCS coordinates
            % Originin GHJC
            TCS_Pos_OrLCS = GHJC;
            % The axes of LCS
            TCS_R_LCS(:,2) = (GHJC-((HME+HLE)/2))/norm(GHJC-((HME+HLE)/2));
            TCS_R_LCS(:,3) = cross(TCS_R_LCS(:,2),YFa)/norm(cross(TCS_R_LCS(:,2),YFa));
            TCS_R_LCS(:,1) = cross(TCS_R_LCS(:,2),TCS_R_LCS(:,3));
            LCS_R_TCS = TCS_R_LCS';
        end
        function [TCS_Pos_OrLCS,LCS_R_TCS] = calcPelvisLCS(S)
            % get ALs indices
            LIASIndex = getLocVecIndex('LIAS',S.LocalALs);
            RIASIndex = getLocVecIndex('RIAS',S.LocalALs);
            LIPSIndex = getLocVecIndex('LIPS',S.LocalALs);
            RIPSIndex = getLocVecIndex('RIPS',S.LocalALs);
            % check if ALs have been measured
            if LIASIndex==0, error('Pelvis landmark LIAS not found in subjec parameter file'); end
            if RIASIndex==0, error('Pelvis landmark RIAS not found in subjec parameter file'); end
            if LIPSIndex==0, error('Pelvis landmark LIPS not found in subjec parameter file'); end
            if RIPSIndex==0, error('Pelvis landmark RIPS not found in subjec parameter file'); end
            % get ALs coordinates in TCS            
            LIAS = S.TCS_Pos_Landmark(:,LIASIndex);
            RIAS = S.TCS_Pos_Landmark(:,RIASIndex);
            LIPS = S.TCS_Pos_Landmark(:,LIPSIndex);
            RIPS = S.TCS_Pos_Landmark(:,RIPSIndex);
            % APPLY LCS DEFINITION to TCS coordinates
            % Origin mid-point between LIAS and RIAS
            TCS_Pos_OrLCS = (LIAS + RIAS)/2;
            % The axes of LCS
            TCS_R_LCS(:,3) = (RIAS - LIAS) / norm(RIAS - LIAS);
            TCS_R_LCS(:,2) = cross(((LIPS+RIPS)/2)-LIAS,TCS_R_LCS(:,3))/norm(cross(((LIPS+RIPS)/2)-LIAS,TCS_R_LCS(:,3)));
            TCS_R_LCS(:,1) = cross(TCS_R_LCS(:,2),TCS_R_LCS(:,3));
            LCS_R_TCS = TCS_R_LCS';
%             S.drawLCSinTCS(TCS_Pos_OrLCS,TCS_R_LCS);
        end
        function [TCS_Pos_OrLCS,LCS_R_TCS] = calcRadiusLCS(S)
           % Check if is Right or Left Clavicle
            if strcmpi(S.Name,'RightRadius')
                PreFix = 'R';
            elseif strcmpi(S.Name,'LeftRadius')
                PreFix = 'L';
            end
            % get ALs Index
            RSPIndex = getLocVecIndex([PreFix,'RSP'],S.LocalALs);
            HLEIndex = getLocVecIndex([PreFix,'HME'],S.TOLocalALs);
            USPIndex = getLocVecIndex([PreFix,'USP'],S.TOLocalALs);
            % check if ALs have been measured
            if RSPIndex==0, error(['Radius landmark ',PreFix,'RSP not found in subjec parameter file']); end
            % get ALs coordinates in TCS
            RSP = S.TCS_Pos_Landmark(:,RSPIndex);
            HLE = S.TCS_Pos_TOALs(:,HLEIndex);
            USP = S.TCS_Pos_TOALs(:,USPIndex);
            
            % APPLY LCS DEFINITION to TCS coordinates
            % Originin RSP
            TCS_Pos_OrLCS = RSP;
            % The axes of LCS
            TCS_R_LCS(:,2) = (HLE-RSP)/norm(HLE-RSP);
            if strcmp(PreFix,'R')
                TCS_R_LCS(:,1) = cross(USP-RSP,HLE-RSP)/norm(cross(USP-RSP,HLE-RSP));
            elseif strcmp(PreFix,'L')
                TCS_R_LCS(:,1) = cross(HLE-RSP,USP-RSP)/norm(cross(HLE-RSP,USP-RSP));
            end
            TCS_R_LCS(:,3) = cross(TCS_R_LCS(:,1),TCS_R_LCS(:,2));
            LCS_R_TCS = TCS_R_LCS'; 
        end
        function [TCS_Pos_OrLCS,LCS_R_TCS] = calcScapulaLCS(S)
            % Check if is Right or Left Scapula
            if strcmpi(S.Name,'RightScapula')
                PreFix = 'R';
            elseif strcmpi(S.Name,'LeftScapula')
                PreFix = 'L';
            end
            % get ALs Index
            SRSIndex = getLocVecIndex([PreFix,'SRS'],S.LocalALs);
            SIAIndex = getLocVecIndex([PreFix,'SIA'],S.LocalALs);
            SAAIndex = getLocVecIndex([PreFix,'SIA'],S.LocalALs);
            % check if ALs have been measured
            if SRSIndex==0, error(['Scapula landmark ',PreFix,'SRS not found in subjec parameter file']); end
            if SIAIndex==0, error(['Scapula landmark ',PreFix,'SIA not found in subjec parameter file']); end
            if SAAIndex==0, error(['Scapula landmark ',PreFix,'SAA not found in subjec parameter file']); end
            % get ALs coordinates in TCS
            SRS = S.TCS_Pos_Landmark(:,SRSIndex);
            SIA = S.TCS_Pos_Landmark(:,SIAIndex);
            SAA = S.TCS_Pos_Landmark(:,SAAIndex);
            % APPLY LCS DEFINITION to TCS coordinates
            % Originin CSJ
            TCS_Pos_OrLCS = SAA;
            % The axes of LCS
            TCS_R_LCS(:,3) = (SAA-SRS)/norm(SAA-SRS);
            if strcmp(PreFix,'R')
                TCS_R_LCS(:,1) = cross(SIA-SAA,SRS-SAA)/norm(cross(SIA-SAA,SRS-SAA));
            elseif strcmp(PreFix,'L')
                TCS_R_LCS(:,1) = cross(SRS-SAA,SIA-SAA)/norm(cross(SRS-SAA,SIA-SAA));
            end
            TCS_R_LCS(:,2) =  cross(TCS_R_LCS(:,3),TCS_R_LCS(:,1));
            LCS_R_TCS = TCS_R_LCS';
        end
        function [TCS_Pos_OrThigh,TCS_Pos_OrLCS,LCS_R_TCS] = calcShankLCS(S)
            
            NMarkers = size (S.LocalMarkers,1);
            if strcmpi(S.Name,'RightShank')
                PreFix = 'R';
            elseif strcmpi(S.Name,'LeftShank')
                PreFix = 'L';
            end
            TAMIndex = getLocVecIndex([PreFix,'TAM'],S.LocalALs);
            FALIndex = getLocVecIndex([PreFix,'FAL'],S.LocalALs);
            FNEIndex = getLocVecIndex([PreFix,'FNE'],S.LocalALs);
            FMEIndex = getLocVecIndex([PreFix,'FME'],S.TOLocalALs);
            FLEIndex = getLocVecIndex([PreFix,'FLE'],S.TOLocalALs);
            
            % check if ALs have been measured
            if TAMIndex==0, error(['Shank landmark ',PreFix,'TAM not found in subjec parameter file']); end
            if FALIndex==0, error(['Shank landmark ',PreFix,'FAL not found in subjec parameter file']); end
            if FNEIndex==0, error(['Shank landmark ',PreFix,'FNE not found in subjec parameter file']); end
            % get AL values in TCS
            TAM = S.TCS_Pos_Landmark(:,TAMIndex);
            FAL = S.TCS_Pos_Landmark(:,FALIndex);
            FNE = S.TCS_Pos_Landmark(:,FNEIndex);
            FME = S.TCS_Pos_TOALs(:,FMEIndex);
            FLE = S.TCS_Pos_TOALs(:,FLEIndex);
            TCS_Pos_OrThigh = (FME+FLE)/2;
            % Calc SHANK LCS --------------------------------
            % Origin mid-point between TAM and FAL
            TCS_Pos_OrLCS = (TAM + FAL)/2;
            % LCS axes
            TCS_R_LCS(:,2) = (TCS_Pos_OrThigh - TCS_Pos_OrLCS)/norm(TCS_Pos_OrThigh - TCS_Pos_OrLCS);
            if strcmp(PreFix,'L') % Left
               TCS_R_LCS(:,1) = cross(FNE-TCS_Pos_OrLCS,TCS_R_LCS(:,2))/norm(cross(FNE-TCS_Pos_OrLCS,TCS_R_LCS(:,2)));
            elseif strcmp(PreFix,'R')
               TCS_R_LCS(:,1) = cross(TCS_R_LCS(:,2),FNE-TCS_Pos_OrLCS)/norm(cross(TCS_R_LCS(:,2),FNE-TCS_Pos_OrLCS)); 
            end
            TCS_R_LCS(:,3) = cross(TCS_R_LCS(:,1),TCS_R_LCS(:,2));
            LCS_R_TCS = TCS_R_LCS';
%             S.drawLCSinTCS(TCS_Pos_OrLCS,TCS_R_LCS);

        end
        function [TCS_Pos_OrLCS,LCS_R_TCS] = calcThighLCS(S,TCS_Pos_HJC)
            if strcmpi(S.Name,'RightThigh')
                PreFix = 'R';
            elseif strcmpi(S.Name,'LeftThigh')
                PreFix = 'L';
            end
            FMEIndex = getLocVecIndex([PreFix,'FME'],S.LocalALs);
            FLEIndex = getLocVecIndex([PreFix,'FLE'],S.LocalALs);
                        
            % check if ALs have been measured
            if FMEIndex==0, error(['Thigh landmark ',PreFix,'FME not found in subjec parameter file']); end
            if FLEIndex==0, error(['Thigh landmark ',PreFix,'FLE not found in subjec parameter file']); end
           
            % Data from TCS in Posture from ALM ------------
            FME = S.TCS_Pos_Landmark(:,FMEIndex);
            FLE = S.TCS_Pos_Landmark(:,FLEIndex);
            
            % Calc SHANK LCS --------------------------------
            % Origin mid-point between FME and FLE
            TCS_Pos_OrLCS = (FME + FLE)/2;
            % The axes of LCS
            TCS_R_LCS(:,2) = (TCS_Pos_HJC - TCS_Pos_OrLCS)/norm(TCS_Pos_HJC - TCS_Pos_OrLCS);
            if strcmp(PreFix,'L') % Left
                TCS_R_LCS(:,1) = cross(FLE-FME,TCS_Pos_HJC-FME)/norm(cross(FLE-FME,TCS_Pos_HJC-FME));
            elseif strcmp(PreFix,'R')% Right
                TCS_R_LCS(:,1) = cross(TCS_Pos_HJC-FME,FLE-FME)/norm(cross(TCS_Pos_HJC-FME,FLE-FME));
            end
            TCS_R_LCS(:,3) = cross(TCS_R_LCS(:,1),TCS_R_LCS(:,2));
            LCS_R_TCS = TCS_R_LCS';
%             S.drawLCSinTCS(TCS_Pos_OrLCS,TCS_R_LCS);
        end
        function [TCS_Pos_OrLCS,LCS_R_TCS] = calcThoraxLCS(S)
            % get ALs index
            SJNIndex = getLocVecIndex('SJN',S.LocalALs);
            SXSIndex = getLocVecIndex('SXS',S.LocalALs);
            TV3Index = getLocVecIndex('TV3',S.LocalALs);
            TV8Index = getLocVecIndex('TV8',S.LocalALs);
            % check if ALs have been measured
            if SJNIndex==0, error('Thorax landmark SJN not found in subjec parameter file'); end
            if SXSIndex==0, error('Thorax landmark SXS not found in subjec parameter file'); end
            if TV3Index==0, error('Thorax landmark TV3 not found in subjec parameter file'); end
            if TV8Index==0, error('Thorax landmark TV8 not found in subjec parameter file'); end
            % get ALs coordinates in TCS            
            SJN = S.TCS_Pos_Landmark(:,SJNIndex);
            SXS = S.TCS_Pos_Landmark(:,SXSIndex);
            TV3 = S.TCS_Pos_Landmark(:,TV3Index);
            TV8 = S.TCS_Pos_Landmark(:,TV8Index);
            % APPLY LCS DEFINITION to TCS coordinates
            % Origin SJN
            TCS_Pos_OrLCS = SJN;
            % The axes of LCS
            TCS_R_LCS(:,2) = (((SJN+TV3)/2)-((SXS-TV8)/2))/norm (((SJN+TV3)/2)-((SXS-TV8)/2));
            TCS_R_LCS(:,3) = cross((TV3-SJN),(((SXS-TV8)/2)-SJN))/ norm (cross((TV3-SJN),(((SXS-TV8)/2)-SJN)));
            TCS_R_LCS(:,1) = cross(TCS_R_LCS(:,2),TCS_R_LCS(:,3));
            LCS_R_TCS = TCS_R_LCS';            
        end
        function [TCS_Pos_OrLCS,LCS_R_TCS] = calcUlnaLCS(S)
            % Check if is Right or Left Clavicle
            if strcmpi(S.Name,'RightUlna')
                PreFix = 'R';
            elseif strcmpi(S.Name,'LeftUlna')
                PreFix = 'L';
            end
            % get ALs Index
            USPIndex = getLocVecIndex([PreFix,'USP'],S.LocalALs);
            HMEIndex = getLocVecIndex([PreFix,'HME'],S.TOLocalALs);
            HLEIndex = getLocVecIndex([PreFix,'HME'],S.TOLocalALs);
            % check if ALs have been measured
            if USPIndex==0, error(['Ulna landmark ',PreFix,'USP not found in subjec parameter file']); end
            % get ALs coordinates in TCS
            USP = S.TCS_Pos_Landmark(:,USPIndex);
            HME = S.TCS_Pos_TOALs(:,HMEIndex);
            HLE = S.TCS_Pos_TOALs(:,HLEIndex);
            % APPLY LCS DEFINITION to TCS coordinates
            % Originin USP
            TCS_Pos_OrLCS = USP;
            % The axes of LCS
            TCS_R_LCS(:,2) = ((HME+HLE/2)-USP)/norm(((HME+HLE/2)-USP));
            if strcmp(PreFix,'R')
                TCS_R_LCS(:,1) = cross(HME-USP,HLE-USP)/norm(cross(HME-USP,HLE-USP));
            elseif strcmp(PreFix,'L')
                TCS_R_LCS(:,1) = cross(HLE-USP,HME-USP)/norm(cross(HLE-USP,HME-USP));
            end
            TCS_R_LCS(:,3) = cross(TCS_R_LCS(:,1),TCS_R_LCS(:,2));
            LCS_R_TCS = TCS_R_LCS';
        end
        function [R, d] = calOptPose(S,Glob_Pos_Markers, TCS_Pos_Markers)
            % CALOPTPOSE Algorithm based on SODERKVIST and WEDIN (1993) that calculates
            % the optimal pose of a body given the global and local coordinates of n points (markers) of the body.
            %
            %   [R, d] = calOptPose(Glob_Pos_Markers, TCS_Pos_Markers)
            %   Inputs:
            %     + Glob_Pos_Markers are the global coordinates of the markers located in
            %       the body - double array (3 x NMarkers)
            %     + TCS_Pos_Markers are the local coordinates of the markers located in
            %       the body referred to TCS - double array(3 x NMarkers)
            %   Outputs:
            %     + R is the rotation matrix between the global coordinate system and
            %       TCS - double array(3x3)
            %     + d is the position vector of the TCS origin referred to global coordinate
            %       system - double array(3x1)

            globMark = Glob_Pos_Markers;
            locMark  = TCS_Pos_Markers;
            
            NMarkers = size(globMark,2);
            
            % 1) -----------------------
            if isempty(locMark)
                meanGlobMark = zeros(3, 1);
                meanLocMark = zeros(3, 1);
            else
                meanGlobMark = (1/NMarkers) * sum(globMark, 2); % sum by rows
                meanLocMark  = (1/NMarkers) * sum(locMark, 2); % sum by rows
            end
            
            % 2) -----------------------
            A = zeros(3, NMarkers);
            for i=1 : NMarkers
                A(:, i) =  locMark(:, i) - meanLocMark;
            end
            B = zeros(3, NMarkers);
            for i=1 : NMarkers
                B(:, i) =  globMark(:, i) - meanGlobMark;
            end
            
            % 3) -----------------------
            C = B * A';
            
            % 4) -----------------------
            [P, Gamma, Q] = svd(C);
            
            % 5) -----------------------
            R = P * diag([1; 1; det(P*Q')]) * Q';
            
            % 6) -----------------------
            d = meanGlobMark - R * meanLocMark;
        end
        function calcInertiaInCoM(S,Ixx,Iyy,Izz,Ixy,Ixz,Iyz,CoM)
            IxxCm = Ixx - S.Mass*(CoM(2)^2+CoM(3)^2);
            IyyCm = Iyy - S.Mass*(CoM(1)^2+CoM(3)^2);
            IzzCm = Izz - S.Mass*(CoM(2)^2+CoM(1)^2);
            IxyCm = Ixy + S.Mass*CoM(1)*CoM(2);
            IxzCm = Ixz + S.Mass*CoM(1)*CoM(3);
            IyzCm = Iyz + S.Mass*CoM(3)*CoM(2);
            S.I =[IxxCm,IxyCm,IxzCm;IxyCm,IyyCm,IyzCm;IxzCm,IyzCm,IzzCm]; 
            
        end
        function calcInertiaPars(S,HumanMass,CoM)
            if strcmp(S.Name, 'Pelvis')
                S.calcInertiaParsPelvis(HumanMass,CoM);
            elseif strcmp (S.Name, 'RightThigh')||strcmp (S.Name, 'LeftThigh')||strcmp (S.Name, 'Thigh')
                S.calcInertiaParsThigh(HumanMass,CoM);
            elseif strcmp (S.Name, 'RightShank')||strcmp (S.Name, 'LeftShank')||strcmp (S.Name, 'Shank')
                S.calcInertiaParsShank(HumanMass,CoM);
            elseif  strcmp (S.Name, 'RightFoot')||strcmp (S.Name, 'LeftFoot')||strcmp (S.Name, 'Foot')
                S.calcInertiaParsFoot(HumanMass,CoM);
            else
                error('%s is not a valid name segment.',S.Name);
            end
        end
        function calcInertiaParsFoot(S,HumanMass,CoM)
            S.CoM=LOCAL_POINT(CoM,S.Name);
            if strcmpi (S.Gender,'male')
                S.Mass = 1.2*HumanMass/100;
                S.CoM.LocCoord = S.Length*[38.2;-15.1;2.6]/100;
                Ixx = 0.17*0.17*S.Length*S.Length*S.Mass;
                Ixy = 0.13*0.13*S.Length*S.Length*S.Mass;
                Ixz = -0.08*0.08*S.Length*S.Length*S.Mass;
                Iyy = 0.37*0.37*S.Length*S.Length*S.Mass;
                Iyz = 0;
                Izz = 0.36*0.36*S.Length*S.Length*S.Mass;
            elseif strcmpi (S.Gender,'female')
                S.Mass = 1.0*HumanMass/100;
                S.CoM.LocCoord = S.Length*[27.0;-21.8;3.9]/100;
                Ixx = 0.17*0.17*S.Length*S.Length*S.Mass;
                Ixy = -0.10*0.10*S.Length*S.Length*S.Mass;
                Ixz = 0.06*0.06*S.Length*S.Length*S.Mass;
                Iyy = 0.36*0.36*S.Length*S.Length*S.Mass;
                Iyz = -0.04*0.04*S.Length*S.Length*S.Mass;
                Izz = 0.35*0.35*S.Length*S.Length*S.Mass;
            else
                error('%s is not a valid gender.',S.Gender);
            end
            S.I =[Ixx,Ixy,Ixz;Ixy,Iyy,Iyz;Ixz,Iyz,Izz];     
        end
        function calcInertiaParsPelvis(S,HumanMass,CoM)
            S.CoM=LOCAL_POINT(CoM,S.Name);
            if strcmpi (S.Gender,'male')
                S.Mass = 14.2*HumanMass/100;
                S.CoM.LocCoord = S.Length*[-97.2;-28.0;-0.6]/100;
                Ixx = 1.01*1.01*S.Length*S.Length*S.Mass;
                Ixy = -0.25*0.25*S.Length*S.Length*S.Mass;
                Ixz = -0.12*0.12*S.Length*S.Length*S.Mass;
                Iyy = 1.06*1.06*S.Length*S.Length*S.Mass;
                Iyz = -0.08*0.08*S.Length*S.Length*S.Mass;
                Izz = 0.95*0.95*S.Length*S.Length*S.Mass;
            elseif strcmpi (S.Gender,'female')
                S.Mass = 14.6*HumanMass/100;
                S.CoM.LocCoord = S.Length*[-100.9;-23.2;0.2]/100;
                Ixx = 0.91*0.91*S.Length*S.Length*S.Mass;
                Ixy = -0.34*0.34*S.Length*S.Length*S.Mass;
                Ixz = -0.01*0.01*S.Length*S.Length*S.Mass;
                Iyy = 1.00*1.00*S.Length*S.Length*S.Mass;
                Iyz = -0.01*0.01*S.Length*S.Length*S.Mass;
                Izz = 0.79*0.79*S.Length*S.Length*S.Mass;
            else
                error('%s is not a valid gender.',S.Gender);
            end
            S.I =[Ixx,Ixy,Ixz;Ixy,Iyy,Iyz;Ixz,Iyz,Izz];
        end
        function calcInertiaParsShank(S,HumanMass,CoM)
            S.CoM=LOCAL_POINT(CoM,S.Name);
            if strcmpi (S.Gender,'male')
                S.Mass = 4.8*HumanMass/100;
                S.CoM.LocCoord = S.Length*[-4.8;59;0.7]/100;
                Ixx = 0.28*0.28*S.Length*S.Length*S.Mass;
                Ixy = -0.04*0.04*S.Length*S.Length*S.Mass;
                Ixz = -0.02*0.02*S.Length*S.Length*S.Mass;
                Iyy = 0.10*0.10*S.Length*S.Length*S.Mass;
                Iyz = 0.05*0.05*S.Length*S.Length*S.Mass;
                Izz = 0.28*0.28*S.Length*S.Length*S.Mass;
            elseif strcmpi (S.Gender,'female')
                S.Mass = 4.5*HumanMass/100;
                S.CoM.LocCoord = S.Length*[-4.9;59.6;3.1]/100;
                Ixx = 0.28*0.28*S.Length*S.Length*S.Mass;
                Ixy = 0.02*0.02*S.Length*S.Length*S.Mass;
                Ixz = 0.01*0.01*S.Length*S.Length*S.Mass;
                Iyy = 0.10*0.10*S.Length*S.Length*S.Mass;
                Iyz = 0.06*0.06*S.Length*S.Length*S.Mass;
                Izz = 0.28*0.28*S.Length*S.Length*S.Mass;
            else
                error('%s is not a valid gender.',S.Gender);
            end
            S.I =[Ixx,Ixy,Ixz;Ixy,Iyy,Iyz;Ixz,Iyz,Izz];                       

        end
        function calcInertiaParsThigh(S,HumanMass,CoM)
            S.CoM=LOCAL_POINT(CoM,S.Name);
            if strcmpi (S.Gender,'male')
                S.Mass = 12.3*HumanMass/100;
                S.CoM.LocCoord = S.Length*[-4.1;57.1;3.3]/100;
                Ixx = 0.29*0.29*S.Length*S.Length*S.Mass;
                Ixy = 0.07*0.07*S.Length*S.Length*S.Mass;
                Ixz = -0.02*0.02*S.Length*S.Length*S.Mass;
                Iyy = 0.15*0.15*S.Length*S.Length*S.Mass;
                Iyz = -0.07*0.07*S.Length*S.Length*S.Mass;
                Izz = 0.30*0.3*S.Length*S.Length*S.Mass;
            elseif strcmpi (S.Gender,'female')
                S.Mass = 14.6*HumanMass/100;
                S.CoM.LocCoord = S.Length*[-7.7;62.3;0.9]/100;
                Ixx = 0.31*0.31*S.Length*S.Length*S.Mass;
                Ixy = 0.07*0.07*S.Length*S.Length*S.Mass;
                Ixz = -0.02*0.02*S.Length*S.Length*S.Mass;
                Iyy = 0.19*0.19*S.Length*S.Length*S.Mass;
                Iyz = -0.07*0.07*S.Length*S.Length*S.Mass;
                Izz = 0.32*0.32*S.Length*S.Length*S.Mass;
            else
                error('%s is not a valid gender.',S.Gender);
            end
            S.I =[Ixx,Ixy,Ixz;Ixy,Iyy,Iyz;Ixz,Iyz,Izz];  
        end 
        function calcTCS1(S)
            % old Option to calculate TCS
            NAnaLandmarks = size(S.LocalALs,1);
            NMarkers = size (S.LocalMarkers,1);
            TCS_Pos_Markers = zeros(3,NMarkers);
            for i=1:NAnaLandmarks
                Glob_Pos_Markers =[];
                for j=1:NMarkers
                    Glob_Pos_Markers =  [Glob_Pos_Markers,S.LocalMarkers(j).Point.MeasuredCoord(:,i)];
                end
                [Glob_Pos_OrCoordSist{i}, Glob_R_BodyCoordSist{i}, BodyCoordSist_Pos_Marker{i}]= calcMarkersTCS(Glob_Pos_Markers);
%                 S.drawGlobTCS(Glob_Pos_OrCoordSist,Glob_R_BodyCoordSist,Glob_Pos_Markers,S.LocalALs(i));
                TCS_Pos_Markers = TCS_Pos_Markers + BodyCoordSist_Pos_Marker{i};
            end
            % % Average local coordinates of the markers in TCS
            S.TCS_Pos_Markers = (1 / NAnaLandmarks) * TCS_Pos_Markers; 
            disp([S.Name]);
%             DispTCS_Pos_Markers = S.TCS_Pos_Markers
%             DispBodyCoordSist_Pos_Marker = BodyCoordSist_Pos_Marker
            % Optimal pose of the "average TCS" using Soderkvist and Wedin algorithm.
            % From this, TCS_Pos_Landmark can be obtained.
            S.TCS_Pos_Landmark =zeros(3,NAnaLandmarks);
%             disp([S.Name]);
            for i=1:NAnaLandmarks
                Glob_Pos_Markers =[];
                for j=1:NMarkers
                    Glob_Pos_Markers =  [Glob_Pos_Markers,S.LocalMarkers(j).Point.MeasuredCoord(:,i)];
                end
                [Glob_R_TCS{i}, Glob_Pos_OrTCS{i}] = S.calOptPose(Glob_Pos_Markers, S.TCS_Pos_Markers);
%                 for k=1:NMarkers
%                     Loc_Pos_Markers(:,k) = Glob_R_TCS{i}' * (Glob_Pos_Markers(:,k) - Glob_Pos_OrTCS{i});
%                 end
                S.TCS_Pos_Landmark(:,i) = Glob_R_TCS{i}' * (S.LocalALs(i).Point.GlobalCoord - Glob_Pos_OrTCS{i});
%                 S.drawGlobTCS1(Glob_Pos_OrTCS{i},Glob_R_TCS{i},S.LocalALs(i),S.TCS_Pos_Markers,Glob_Pos_Markers)
            end
%             DispOpt_Pos_Markers = Loc_Pos_Markers
%             DispTCS_Pos_ALs = S.TCS_Pos_Landmark
            
        end
        function calcTCS(S)
            NAnaLandmarks = size(S.LocalALs,1);
            NMarkers = size (S.LocalMarkers,1);
%             TCS_Pos_Markers = zeros(3,NMarkers);
%             Solo utilizo el frame del primer landmark palpado cambiar 1 por i            
%             for i=1:NAnaLandmarks
                Glob_Pos_Markers =[];
                for j=1:NMarkers
                    Glob_Pos_Markers =  [Glob_Pos_Markers,S.LocalMarkers(j).Point.MeasuredCoord(:,1)];
                end
                
                [Glob_Pos_OrCoordSist,Glob_R_BodyCoordSist,BodyCoordSist_Pos_Marker]= calcMarkersTCS(Glob_Pos_Markers);
                %S.drawGlobTCS(Glob_Pos_OrCoordSist,Glob_R_BodyCoordSist,Glob_Pos_Markers,S.LocalALs(1));
                % local coordinates of the markers in TCS
                
%             disp([S.Name]);
            S.TCS_Pos_Markers = BodyCoordSist_Pos_Marker;
%             DispTCS_Pos_Markers = BodyCoordSist_Pos_Marker
%             end
%             % % Average local coordinates of the markers in TCS
%             S.TCS_Pos_Markers = (1 / NAnaLandmarks) * TCS_Pos_Markers; 
            % Optimal pose of the "average TCS" using Soderkvist and Wedin algorithm.
            % From this, TCS_Pos_Landmark can be obtained.
            S.TCS_Pos_Landmark =zeros(3,NAnaLandmarks);
%             disp([S.Name]);
            for i=1:NAnaLandmarks
                Glob_Pos_Markers =[];
                for j=1:NMarkers
                    Glob_Pos_Markers =  [Glob_Pos_Markers,S.LocalMarkers(j).Point.MeasuredCoord(:,i)];
                end
                %                 [Glob_R_TCS{i}, Glob_Pos_OrTCS{i}] = S.calOptPose(Glob_Pos_Markers, S.TCS_Pos_Markers);
                %                 S.TCS_Pos_Landmark(:,i) = Glob_R_TCS{i}' * (S.LocalALs(i).Point.GlobalCoord - Glob_Pos_OrTCS{i});
                [Glob_R_TCS{i}, Glob_Pos_OrTCS{i}] = S.calOptPose(Glob_Pos_Markers, S.TCS_Pos_Markers);
%                 for k=1:NMarkers
%                     Loc_Pos_Markers(:,k) = Glob_R_TCS{i}' * (Glob_Pos_Markers(:,k) - Glob_Pos_OrTCS{i});
%                 end
                S.TCS_Pos_Landmark(:,i) = Glob_R_TCS{i}' * (S.LocalALs(i).Point.GlobalCoord - Glob_Pos_OrTCS{i});
%                 S.drawGlobTCS1(Glob_Pos_OrTCS{i},Glob_R_TCS{i},S.LocalALs(i),Loc_Pos_Markers,Glob_Pos_OrCoordSist,Glob_R_BodyCoordSist,BodyCoordSist_Pos_Marker,Glob_Pos_Markers)
%                 S.drawGlobTCS1(Glob_Pos_OrTCS{i},Glob_R_TCS{i},S.LocalALs(i),BodyCoordSist_Pos_Marker,Glob_Pos_Markers)
%                 DispGlob_R_TCS = Glob_R_TCS{i}
            end
%             DispOpt_Pos_Markers = Loc_Pos_Markers
%             DispTCS_Pos_ALs = S.TCS_Pos_Landmark
            
        end
        function calcTransferedALs(S,Post_Pos_TOSegMarkers,TCSTO_Pos_TOSegMarkers,TCSTO_Pos_TOSegALs,TOLocalALs,PostureName)
          % Post_Pos_TOSegMarkers  = Position of the markers of transfer origin segment in the Posture (Global) coordinate system
          % TCSTO_Pos_TOSegMarkers = Position of the markers of transfer origin segment in the transfer origin segment TCS
          % TCSTO_Pos_TOSegALs     = Position of the AL of transfer origin segment in the transfer origin segment TCS
          % TOLocalALs             = Vector with all the ALs in the transfer origin segment
          % PostureName            = The name of the Posture when the translation is measure
            %calcTransferredALs
            NTOSegments = size(Post_Pos_TOSegMarkers,2);
            for i=1:NTOSegments
                if ~isempty(Post_Pos_TOSegMarkers{i})
                    % get ALs coordinates in TO-TCS
                    NTOLocalALs = size(TOLocalALs{i},1);
                    TCSTO_Pos_TOALs = [];
                    for j=1:NTOLocalALs
                        TCSTO_Pos_TOALs = [TCSTO_Pos_TOALs,TCSTO_Pos_TOSegALs{i}(:,j)];
                    end
                    % Pass ALs to the Posture(Global) CS
                    [Glob_R_TCSTO,Glob_Pos_OrTCSTO] = S.calOptPose(Post_Pos_TOSegMarkers{i},TCSTO_Pos_TOSegMarkers{i});
                    TCSTO_Pos_OrGlob = - Glob_R_TCSTO' * Glob_Pos_OrTCSTO;
                    Glob_Pos_TOALs = changeCoordSys(TCSTO_Pos_TOALs,TCSTO_Pos_OrGlob,Glob_R_TCSTO);
                    % Pass Markers of the Segment to the Posture (Global) CS
                    PostureIndex = getVecIndex(PostureName,S.LocalMarkers(1).Point.Postures);
                    NMarkers = size(S.LocalMarkers,1);
                    for j=1:NMarkers
                        Glob_Pos_MarkersPosture(:,j)=S.LocalMarkers(j).Point.Postures(PostureIndex).Glob;
                    end
                    % Share the transfered origin segment AL in the segment TCS
                    [Glob_R_TCS,Glob_Pos_OrTCS]=S.calOptPose(Glob_Pos_MarkersPosture,S.TCS_Pos_Markers);
                    S.TOLocalALs =[S.TOLocalALs;TOLocalALs{i}];
                    k = size(S.TCS_Pos_TOALs,2);
                    for j=1:NTOLocalALs
                        S.TCS_Pos_TOALs(:,j+k) = Glob_R_TCS' * (Glob_Pos_TOALs(:,j) - Glob_Pos_OrTCS);
                    end
                end
            end
                        % %                 S.drawPosture(R,d,Glob_R_Shank,Glob_Pos_OrShank,TCSTO_Pos_TOSegMarkers,S.TCS_Pos_Markers,...
%                               TCSTO_Pos_KJC,TCS_Pos_KJC,'LeftThigh',Post_Pos_TOSegMarkers,Glob_Pos_MarkersPosture,Glob_Pos_
%                               KJC,ParentMarkers)
        end
        function checkSubjectParsData(S,File)
            % CHECKSUBJECTPARSDATA check the if exist the points, vectors and markers of the segment in subject parameter file
            NPoints  = size(S.LocalPoints,1);
            NVectors = size(S.LocalVectors,1);
            NMarkers = size(S.LocalMarkers,1);
            for i=1:NPoints
                if isempty(S.LocalPoints(i).LocCoord)
                    error(['Coords of point "',S.LocalPoints(i).Name,' "not defined in file: ',File])
                end
            end
            for i=1:NVectors
                if isempty(S.LocalVectors(i).LocCoord)
                    error(['Coords of vector "',S.LocalVectors(i).Name,' "not defined in file: ',File])
                end
            end
            for i=1:NMarkers
                if isempty(S.LocalMarkers(i).LocCoord)
                    error(['The coord of marker "',S.LocalMarkers(i).Name,' "not defined in file: ',File])
                end
            end
        end
        function drawGlobTCS(S,Glob_Pos_OrCoordSist,Glob_R_BodyCoordSist,Glob_Pos_Markers,AL)
            NMarkers = size(S.LocalMarkers,1);
            NLandMarks = size(S.LocalALs,1);

            figure('Name',['Palpation ',S.Name,': ',AL.Point.Name],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
           %xlabel(Var_x,'Interpreter','none','FontSize',FontSize);
           %ylabel(Var_y,'Interpreter','none','FontSize',FontSize);

            %figure;
            % figure settings
            hold on, axis equal, grid on
            % calculate figure scale
            MaxSmark = max(max(abs(Glob_Pos_Markers)));
            MaxLmark = max(abs(AL.Point.GlobalCoord));
            MaxMarks = [MaxSmark MaxLmark];
            Scale = 0.5 * max(MaxMarks);
             % plot Global axes 
            quiver3(0,0,0,Scale,0,0,0,'k');
            quiver3(0,0,0,0,Scale,0,0,'k');
            quiver3(0,0,0,0,0,Scale,0,'k');
            
            % plot axes names
            text(Scale,0,0,'X   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % X-axis
            text(0,Scale,0,'Y   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Y-axis
            text(0,0,Scale,'Z   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Z-axis
           
            % plot TCS in Global
            quiver3(Glob_Pos_OrCoordSist(1),Glob_Pos_OrCoordSist(2),Glob_Pos_OrCoordSist(3),...
                    Glob_R_BodyCoordSist(1,1),Glob_R_BodyCoordSist(2,1),Glob_R_BodyCoordSist(3,1),Scale,'r');
            quiver3(Glob_Pos_OrCoordSist(1),Glob_Pos_OrCoordSist(2),Glob_Pos_OrCoordSist(3),...
                    Glob_R_BodyCoordSist(1,2),Glob_R_BodyCoordSist(2,2),Glob_R_BodyCoordSist(3,2),Scale,'g');
            quiver3(Glob_Pos_OrCoordSist(1),Glob_Pos_OrCoordSist(2),Glob_Pos_OrCoordSist(3),...
                    Glob_R_BodyCoordSist(1,3),Glob_R_BodyCoordSist(2,3),Glob_R_BodyCoordSist(3,3),Scale,'b');
                
            % plot axes names
%             text(Scale*Glob_R_BodyCoordSist(1,1),Glob_R_BodyCoordSist(2,1),Glob_R_BodyCoordSist(3,1),'X_T_C_S   ',...
%                 'FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % X-axis
%             text(Glob_R_BodyCoordSist(1,2),Scale*Glob_R_BodyCoordSist(2,2),Glob_R_BodyCoordSist(3,2),'Y_T_C_S   ',...
%                 'FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Y-axis
%             text(Glob_R_BodyCoordSist(1,3),Glob_R_BodyCoordSist(2,3),Scale*Glob_R_BodyCoordSist(3,3),'Z_T_C_S   ',...
%                 'FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Z-axis
%             text(Glob_Pos_OrCoordSist(1)+Glob_R_BodyCoordSist(1,1)*Scale,Glob_Pos_OrCoordSist(2)+Glob_R_BodyCoordSist(2,1)*Scale,Glob_Pos_OrCoordSist(3)+Glob_R_BodyCoordSist(1,3)*Scale,'X_T_C_S   ',...
%                 'FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % X-axis
%             text(Glob_R_BodyCoordSist(1,2),Scale*Glob_R_BodyCoordSist(2,2),Glob_R_BodyCoordSist(3,2),'Y_T_C_S   ',...
%                 'FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Y-axis
%             text(Glob_R_BodyCoordSist(1,3),Glob_R_BodyCoordSist(2,3),Scale*Glob_R_BodyCoordSist(3,3),'Z_T_C_S   ',...
%                 'FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Z-axis

            % plot markers in globals
            for i=1:NMarkers
                hold on;
                plot3(Glob_Pos_Markers(1,i), Glob_Pos_Markers(2,i), Glob_Pos_Markers(3,i), 'bo');
                text(Glob_Pos_Markers(1,i), Glob_Pos_Markers(2,i), Glob_Pos_Markers(3,i), ...
                    ['  ', S.LocalMarkers(i).Point.Name], 'FontSize', 10, 'Color', 'b');
                line([Glob_Pos_Markers(1,i),0],[Glob_Pos_Markers(2,i),0],[Glob_Pos_Markers(3,i),0],'Color','b');
            end
            % plot landmark in global
            hold on;
            plot3(AL.Point.GlobalCoord(1),AL.Point.GlobalCoord(2),AL.Point.GlobalCoord(3),'rx');
            line([AL.Point.GlobalCoord(1),0],[AL.Point.GlobalCoord(2),0],[AL.Point.GlobalCoord(3),0],'Color','r');
            text(AL.Point.GlobalCoord(1),AL.Point.GlobalCoord(2),AL.Point.GlobalCoord(3), ...
                    ['  ', AL.Point.Name], 'FontSize', 10, 'Color', 'r');
            % calc distance between markers
            cont=1;
            for i=1:NMarkers-1
                for j=NMarkers:-1:i+1
                    Dist= norm(Glob_Pos_Markers(:,j)-Glob_Pos_Markers(:,i));
                    str(cont)={['Dist  ', S.LocalMarkers(j).Point.Name,'-',S.LocalMarkers(i).Point.Name,' = ',num2str(Dist)]};
                    cont = cont+1;
                end
            end
            text(Glob_Pos_OrCoordSist(1),Glob_Pos_OrCoordSist(2)+Scale,Glob_Pos_OrCoordSist(3)+Scale,str,'FontSize', 10, 'Color', 'k')
            % title
            title([AL.Point.Name,' & ',S.Name,' TCS. Markers in blue; Landmarks in red'], 'FontWeight', 'Bold', 'FontSize', 10);
%             title(str, 'FontWeight', 'Bold', 'FontSize', 12);
        view([43 46]);
        end
        function drawGlobTCS1(S,Glob_Pos_OrOptimal,Glob_R_Optimal,AL,TCS_Pos_Markers,Glob_Pos_Markers)
            NMarkers = size(S.LocalMarkers,1);
            NLandMarks = size(S.LocalALs,1);

            figure('Name',['Palpation ',S.Name,': ',AL.Point.Name],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
           %xlabel(Var_x,'Interpreter','none','FontSize',FontSize);
           %ylabel(Var_y,'Interpreter','none','FontSize',FontSize);

            %figure;
            % figure settings
            hold on, axis equal, grid on
            % calculate figure scale
            MaxSmark = max(max(abs(Glob_Pos_Markers)));
            MaxLmark = max(abs(AL.Point.GlobalCoord));
            MaxMarks = [MaxSmark MaxLmark];
            Scale = 0.5 * max(MaxMarks);
             % plot Global axes 
            quiver3(0,0,0,Scale,0,0,0,'k');
            quiver3(0,0,0,0,Scale,0,0,'k');
            quiver3(0,0,0,0,0,Scale,0,'k');
            
            % plot axes names
            text(Scale,0,0,'X   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % X-axis
            text(0,Scale,0,'Y   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Y-axis
            text(0,0,Scale,'Z   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Z-axis
           
            % plot TCS in Global
%             quiver3(Glob_Pos_OrTCS(1),Glob_Pos_OrTCS(2),Glob_Pos_OrTCS(3),...
%                     Glob_R_TCS(1,1),Glob_R_TCS(2,1),Glob_R_TCS(3,1),Scale,'r');
%             quiver3(Glob_Pos_OrTCS(1),Glob_Pos_OrTCS(2),Glob_Pos_OrTCS(3),...
%                     Glob_R_TCS(1,2),Glob_R_TCS(2,2),Glob_R_TCS(3,2),Scale,'g');
%             quiver3(Glob_Pos_OrTCS(1),Glob_Pos_OrTCS(2),Glob_Pos_OrTCS(3),...
%                     Glob_R_TCS(1,3),Glob_R_TCS(2,3),Glob_R_TCS(3,3),Scale,'b');
            
            % plot Optimal in Global    
            quiver3(Glob_Pos_OrOptimal(1),Glob_Pos_OrOptimal(2),Glob_Pos_OrOptimal(3),...
                    Glob_R_Optimal(1,1),Glob_R_Optimal(2,1),Glob_R_Optimal(3,1),Scale,'r');
            quiver3(Glob_Pos_OrOptimal(1),Glob_Pos_OrOptimal(2),Glob_Pos_OrOptimal(3),...
                    Glob_R_Optimal(1,2),Glob_R_Optimal(2,2),Glob_R_Optimal(3,2),Scale,'g');
            quiver3(Glob_Pos_OrOptimal(1),Glob_Pos_OrOptimal(2),Glob_Pos_OrOptimal(3),...
                    Glob_R_Optimal(1,3),Glob_R_Optimal(2,3),Glob_R_Optimal(3,3),Scale,'b');
                
            % plot axes names
%             text(Scale*Glob_R_BodyCoordSist(1,1),Glob_R_BodyCoordSist(2,1),Glob_R_BodyCoordSist(3,1),'X_T_C_S   ',...
%                 'FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % X-axis
%             text(Glob_R_BodyCoordSist(1,2),Scale*Glob_R_BodyCoordSist(2,2),Glob_R_BodyCoordSist(3,2),'Y_T_C_S   ',...
%                 'FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Y-axis
%             text(Glob_R_BodyCoordSist(1,3),Glob_R_BodyCoordSist(2,3),Scale*Glob_R_BodyCoordSist(3,3),'Z_T_C_S   ',...
%                 'FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Z-axis
%             text(Glob_Pos_OrCoordSist(1)+Glob_R_BodyCoordSist(1,1)*Scale,Glob_Pos_OrCoordSist(2)+Glob_R_BodyCoordSist(2,1)*Scale,Glob_Pos_OrCoordSist(3)+Glob_R_BodyCoordSist(1,3)*Scale,'X_T_C_S   ',...
%                 'FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % X-axis
%             text(Glob_R_BodyCoordSist(1,2),Scale*Glob_R_BodyCoordSist(2,2),Glob_R_BodyCoordSist(3,2),'Y_T_C_S   ',...
%                 'FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Y-axis
%             text(Glob_R_BodyCoordSist(1,3),Glob_R_BodyCoordSist(2,3),Scale*Glob_R_BodyCoordSist(3,3),'Z_T_C_S   ',...
%                 'FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Z-axis

            % plot markers in globals
            for i=1:NMarkers
                hold on;
                plot3(Glob_Pos_Markers(1,i), Glob_Pos_Markers(2,i), Glob_Pos_Markers(3,i), 'kx');
                text(Glob_Pos_Markers(1,i), Glob_Pos_Markers(2,i), Glob_Pos_Markers(3,i), ...
                    ['  ', S.LocalMarkers(i).Point.Name], 'FontSize', 10, 'Color', 'k');
                line([Glob_Pos_Markers(1,i),0],[Glob_Pos_Markers(2,i),0],[Glob_Pos_Markers(3,i),0],'Color','k');
            end
       
            % plot markers in TCS
%             for i=1:NMarkers
%                 hold on;
%                 plot3(Glob_Pos_OrTCS(1)+Glob_R_TCS(1,:)*TCS_Pos_Markers(:,i),...
%                       Glob_Pos_OrTCS(2)+Glob_R_TCS(2,:)*TCS_Pos_Markers(:,i),...
%                       Glob_Pos_OrTCS(3)+Glob_R_TCS(3,:)*TCS_Pos_Markers(:,i),'mo')
%                 text( Glob_Pos_OrTCS(1)+Glob_R_TCS(1,:)*TCS_Pos_Markers(:,i),...
%                       Glob_Pos_OrTCS(2)+Glob_R_TCS(2,:)*TCS_Pos_Markers(:,i),...
%                       Glob_Pos_OrTCS(3)+Glob_R_TCS(3,:)*TCS_Pos_Markers(:,i),...
%                       ['  ', S.LocalMarkers(i).Point.Name], 'FontSize', 10, 'Color', 'm');
%                 line([Glob_Pos_OrTCS(1)+Glob_R_TCS(1,:)*TCS_Pos_Markers(:,i),Glob_Pos_OrTCS(1)],...
%                      [Glob_Pos_OrTCS(2)+Glob_R_TCS(2,:)*TCS_Pos_Markers(:,i),Glob_Pos_OrTCS(2)],...
%                      [Glob_Pos_OrTCS(3)+Glob_R_TCS(3,:)*TCS_Pos_Markers(:,i),Glob_Pos_OrTCS(3)],...
%                      'Color','m');  
%             end
             % plot markers in Optimal
            for i=1:NMarkers
                hold on;
                plot3(Glob_Pos_OrOptimal(1)+Glob_R_Optimal(1,:)*TCS_Pos_Markers(:,i),...
                      Glob_Pos_OrOptimal(2)+Glob_R_Optimal(2,:)*TCS_Pos_Markers(:,i),...
                      Glob_Pos_OrOptimal(3)+Glob_R_Optimal(3,:)*TCS_Pos_Markers(:,i),'co')
                text( Glob_Pos_OrOptimal(1)+Glob_R_Optimal(1,:)*TCS_Pos_Markers(:,i),...
                      Glob_Pos_OrOptimal(2)+Glob_R_Optimal(2,:)*TCS_Pos_Markers(:,i),...
                      Glob_Pos_OrOptimal(3)+Glob_R_Optimal(3,:)*TCS_Pos_Markers(:,i),...
                      ['  ', S.LocalMarkers(i).Point.Name], 'FontSize', 10, 'Color', 'c');
                line([Glob_Pos_OrOptimal(1)+Glob_R_Optimal(1,:)*TCS_Pos_Markers(:,i),Glob_Pos_OrOptimal(1)],...
                     [Glob_Pos_OrOptimal(2)+Glob_R_Optimal(2,:)*TCS_Pos_Markers(:,i),Glob_Pos_OrOptimal(2)],...
                     [Glob_Pos_OrOptimal(3)+Glob_R_Optimal(3,:)*TCS_Pos_Markers(:,i),Glob_Pos_OrOptimal(3)],...
                     'Color','c');  
            end
            % plot landmark in global
            hold on;
            plot3(AL.Point.GlobalCoord(1),AL.Point.GlobalCoord(2),AL.Point.GlobalCoord(3),'rx');
            line([AL.Point.GlobalCoord(1),0],[AL.Point.GlobalCoord(2),0],[AL.Point.GlobalCoord(3),0],'Color','r');
            text(AL.Point.GlobalCoord(1),AL.Point.GlobalCoord(2),AL.Point.GlobalCoord(3), ...
                    ['  ', AL.Point.Name], 'FontSize', 10, 'Color', 'r');
            % calc distance between markers
            cont=1;
            for i=1:NMarkers-1
                for j=NMarkers:-1:i+1
                    Dist= norm(Glob_Pos_Markers(:,j)-Glob_Pos_Markers(:,i));
                    str(cont)={['Dist  ', S.LocalMarkers(j).Point.Name,'-',S.LocalMarkers(i).Point.Name,' = ',num2str(Dist)]};
                    cont = cont+1;
                end
            end
            text(Glob_Pos_OrOptimal(1),Glob_Pos_OrOptimal(2)+Scale,Glob_Pos_OrOptimal(3)+Scale,str,'FontSize', 10, 'Color', 'k')
            % title
            title([AL.Point.Name,' & ',S.Name,' TCS. Markers in blue; Landmarks in red'], 'FontWeight', 'Bold', 'FontSize', 10);
%             title(str, 'FontWeight', 'Bold', 'FontSize', 12);
        view([43 46]);
        end
        function drawLCS(S)
            % extract input data
            BodyName = S.Name;
            NMarkers = size(S.LocalMarkers,1);
            NLandMarks = size(S.LocalALs,1);
            NPoints = size(S.LocalPoints,1);
            BodyBAF_Pos_Smark = [];
            BodyBAF_Pos_Lmark = [];
            for i=1:NMarkers
                BodyBAF_Pos_Smark = [BodyBAF_Pos_Smark,S.LocalMarkers(i).LocCoord];
            end
            for i=1:NLandMarks
                BodyBAF_Pos_Lmark = [BodyBAF_Pos_Lmark,S.LocalALs(i).LocCoord];
            end
            for i=1:NPoints
                BodyBAF_Pos_Lmark = [BodyBAF_Pos_Lmark,S.LocalPoints(i).LocCoord];
            end
            BodyCM            = S.CoM.LocCoord;
%             JointAxes         = [1,0,0;0,1,0;0,0,1];
%             JointAxesName{1,1}     = [S.Name,'X'];
%             JointAxesName{2,1}     = [S.Name,'Y'];
%             JointAxesName{3,1}     = [S.Name,'Z'];Local axes of the ', BodyName,
            figure('Name',['LCS of ',S.Name],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%             figure;
            % figure settings
            hold on, axis equal, grid on
            % Initialize variable
%             BodyProps = initBodyProps;
            
                      
            % get index body
%             IndexBody = getBodyIndex(BodyProps, BodyName);
            
            % calculate figure scale
            MaxSmark = max(max(abs(BodyBAF_Pos_Smark)));
            MaxLmark = max(max(abs(BodyBAF_Pos_Lmark)));
            MaxMarks = [MaxSmark MaxLmark];
            ScaleBAF = 1.2 * max(MaxMarks);
            
            % plot axes vectors
            quiver3(0,0,0,ScaleBAF,0,0,0,'k');
            quiver3(0,0,0,0,ScaleBAF,0,0,'k');
            quiver3(0,0,0,0,0,ScaleBAF,0,'k');
            
            % plot axes names
            text(ScaleBAF,0,0,'X_L_C_S   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % X-axis
            text(0,ScaleBAF,0,'Y_L_C_S   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Y-axis
            text(0,0,ScaleBAF,'Z_L_C_S   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Z-axis
            
            % Plot calculated skin-markers referred to LCS
            for i = 1 : NMarkers
                hold on;
                plot3(BodyBAF_Pos_Smark(1,i), BodyBAF_Pos_Smark(2,i), BodyBAF_Pos_Smark(3,i), 'bo');
                text(BodyBAF_Pos_Smark(1,i), BodyBAF_Pos_Smark(2,i), BodyBAF_Pos_Smark(3,i), ...
                    ['  ', S.LocalMarkers(i).Point.Name], 'FontSize', 10, 'Color', 'b');
                line([BodyBAF_Pos_Smark(1,i),0],[BodyBAF_Pos_Smark(2,i),0],[BodyBAF_Pos_Smark(3,i),0],'Color','b');
                disp([S.LocalMarkers(i).Point.Name,'=[ ',num2str(BodyBAF_Pos_Smark(1:3,i)','%-11.5f'),']'])
            end
            
            % Plot landmarks referred to LCS
            for i = 1 : NLandMarks
                hold on;
                plot3(BodyBAF_Pos_Lmark(1,i), BodyBAF_Pos_Lmark(2,i), BodyBAF_Pos_Lmark(3,i), 'rx');
                text(BodyBAF_Pos_Lmark(1,i), BodyBAF_Pos_Lmark(2,i), BodyBAF_Pos_Lmark(3,i), ...
                    [' ', S.LocalALs(i).Point.Name, '  '], 'FontSize', 10, 'FontWeight', 'Bold', 'Color', 'r', ...
                    'HorizontalAlignment', 'Right');
                line([BodyBAF_Pos_Lmark(1,i),0],[BodyBAF_Pos_Lmark(2,i),0],[BodyBAF_Pos_Lmark(3,i),0],'Color','r');
            end
            for i = NLandMarks+1 : (NPoints+NLandMarks)
                hold on;
                plot3(BodyBAF_Pos_Lmark(1,i), BodyBAF_Pos_Lmark(2,i), BodyBAF_Pos_Lmark(3,i), 'rx');
                text(BodyBAF_Pos_Lmark(1,i), BodyBAF_Pos_Lmark(2,i), BodyBAF_Pos_Lmark(3,i), ...
                    [' ', S.LocalPoints(i-NLandMarks).Point.Name, '  '], 'FontSize', 10, 'FontWeight', 'Bold', 'Color', 'r', ...
                    'HorizontalAlignment', 'Right');
                line([BodyBAF_Pos_Lmark(1,i),0],[BodyBAF_Pos_Lmark(2,i),0],[BodyBAF_Pos_Lmark(3,i),0],'Color','r');
                disp([S.LocalPoints(i-NLandMarks).Point.Name,'=[ ',num2str(BodyBAF_Pos_Lmark(1:3,i)','%-11.5f'),']'])
            end
            
            % Plot Center of Mass
            
            plot3(BodyCM(1),BodyCM(2),BodyCM(3), 'g^');
            text(BodyCM(1),BodyCM(2),BodyCM(3), '  CM', 'FontSize', 10, 'FontWeight', 'Bold', 'Color', 'g');
            line([BodyCM(1),0],[BodyCM(2),0],[BodyCM(3),0],'Color','g');
            disp([S.Name,'CoM=[ ',num2str(BodyCM','%-11.5f'),']'])
            
            
            % Plot Joint axes
%             quiver3(0,0,0,JointAxes(1,1),JointAxes(2,1),JointAxes(3,1));
%             text(JointAxes(1,1),JointAxes(2,1),JointAxes(3,1),JointAxesName{1},'FontSize',10,'FontWeight','Bold','Color','k','HorizontalAlignment','center');
%             quiver3(0,0,0,JointAxes(1,2),JointAxes(2,2),JointAxes(3,2));
%             text(JointAxes(1,2),JointAxes(2,2),JointAxes(3,2),JointAxesName{2},'FontSize',10,'FontWeight','Bold','Color','k','HorizontalAlignment','center');
%             quiver3(0,0,0,JointAxes(1,1),JointAxes(2,1),JointAxes(3,1));
%             text(JointAxes(1,3),JointAxes(2,3),JointAxes(3,3),JointAxesName{3},'FontSize',10,'FontWeight','Bold','Color','k','HorizontalAlignment','center');
%             if isempty(JointAxes)
%                 
%             elseif strcmp(BodyName,'RThigh')
%                 quiver3(BodyBAF_Pos_Lmark(1,4),BodyBAF_Pos_Lmark(2,4),BodyBAF_Pos_Lmark(3,4),JointAxes(1)*ScaleBAF,JointAxes(2)*ScaleBAF,JointAxes(3)*ScaleBAF);
%                 text((BodyBAF_Pos_Lmark(1,4) + JointAxes(1)*ScaleBAF),(BodyBAF_Pos_Lmark(2,4) + JointAxes(2)*ScaleBAF),(BodyBAF_Pos_Lmark(3,4) + JointAxes(3)*ScaleBAF),JointAxesName,...
%                     'FontSize',10,'FontWeight','Bold','Color','k','HorizontalAlignment','center');
%             elseif strcmp(BodyName,'RShank') | strcmp(BodyName,'LFarm') | strcmp(BodyName,'RFarm')
%                 quiver3(BodyBAF_Pos_Lmark(1,1),BodyBAF_Pos_Lmark(2,1),BodyBAF_Pos_Lmark(3,1),JointAxes(1)*ScaleBAF,JointAxes(2)*ScaleBAF,JointAxes(3)*ScaleBAF);
%                 text((BodyBAF_Pos_Lmark(1,1) + JointAxes(1)*ScaleBAF),(BodyBAF_Pos_Lmark(2,1) + JointAxes(2)*ScaleBAF),(BodyBAF_Pos_Lmark(3,1) + JointAxes(3)*ScaleBAF),JointAxesName,...
%                     'FontSize',10,'FontWeight','Bold','Color','k','HorizontalAlignment','center');
%             elseif strcmp(BodyName,'RUarm') | strcmp(BodyName,'LUarm')
%                 quiver3(BodyBAF_Pos_Lmark(1,4),BodyBAF_Pos_Lmark(2,4),BodyBAF_Pos_Lmark(3,4),JointAxes(1,1)*ScaleBAF,JointAxes(2,1)*ScaleBAF,JointAxes(3,1)*ScaleBAF);
%                 text((BodyBAF_Pos_Lmark(1,4) + JointAxes(1,1)*ScaleBAF),(BodyBAF_Pos_Lmark(2,4) + JointAxes(2,1)*ScaleBAF),(BodyBAF_Pos_Lmark(3,4) + JointAxes(3,1)*ScaleBAF),JointAxesName{1},...
%                     'FontSize',10,'FontWeight','Bold','Color','k','HorizontalAlignment','center');
%                 quiver3(BodyBAF_Pos_Lmark(1,4),BodyBAF_Pos_Lmark(2,4),BodyBAF_Pos_Lmark(3,4),JointAxes(1,2)*ScaleBAF,JointAxes(2,2)*ScaleBAF,JointAxes(3,2)*ScaleBAF);
%                 text((BodyBAF_Pos_Lmark(1,4) + JointAxes(1,2)*ScaleBAF),(BodyBAF_Pos_Lmark(2,4) + JointAxes(2,2)*ScaleBAF),(BodyBAF_Pos_Lmark(1,4) + JointAxes(3,2)*ScaleBAF),JointAxesName{2},...
%                     'FontSize',10,'FontWeight','Bold','Color','k','HorizontalAlignment','center');
%             end
            
            title(['Local axes of the ', BodyName, '. Skin-markers in blue; Landmarks in red'], 'FontWeight', 'Bold', 'FontSize', 12);

            
        end
        function drawLCSinTCS(S,TCS_Pos_OrLCS,TCS_R_LCS)
            NMarkers = size(S.LocalMarkers,1);
            NLandMarks = size(S.LocalALs,1);
            
           %figure('Name',[Var_x,' vs ',Var_y],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',FontSize), hold on
           %xlabel(Var_x,'Interpreter','none','FontSize',FontSize);
           %ylabel(Var_y,'Interpreter','none','FontSize',FontSize);

            figure('Name',[S.Name ' LCS in TCS'],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%             figure;
            % figure settings
            hold on, axis equal, grid on
            % calculate figure scale
            MaxSmark = max(max(abs(S.TCS_Pos_Markers)));
            MaxLmark = max(max(abs(S.TCS_Pos_Landmark)));
            MaxMarks = [MaxSmark MaxLmark];
            Scale = 1.2 * max(MaxMarks);
            
            % plot axes vectors (TCS)
            quiver3(0,0,0,Scale,0,0,0,'k');
            quiver3(0,0,0,0,Scale,0,0,'k');
            quiver3(0,0,0,0,0,Scale,0,'k');
            
            % plot axes names
            text(Scale,0,0,'X_T_C_S   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % X-axis
            text(0,Scale,0,'Y_T_C_S   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Y-axis
            text(0,0,Scale,'Z_T_C_S   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Z-axis
            
            % plot markers
            for i = 1 : NMarkers
                hold on;
                plot3(S.TCS_Pos_Markers(1,i), S.TCS_Pos_Markers(2,i), S.TCS_Pos_Markers(3,i), 'bo');
                text(S.TCS_Pos_Markers(1,i), S.TCS_Pos_Markers(2,i), S.TCS_Pos_Markers(3,i), ...
                    ['  ', S.LocalMarkers(i).Point.Name], 'FontSize', 10, 'Color', 'b');
                line([S.TCS_Pos_Markers(1,i),0],[S.TCS_Pos_Markers(2,i),0],[S.TCS_Pos_Markers(3,i),0],'Color','b');
            end
            % Plot landmarks referred to LCS
            for i = 1 : NLandMarks
                hold on;
                plot3(S.TCS_Pos_Landmark(1,i), S.TCS_Pos_Landmark(2,i), S.TCS_Pos_Landmark(3,i), 'mx');
                text(S.TCS_Pos_Landmark(1,i), S.TCS_Pos_Landmark(2,i), S.TCS_Pos_Landmark(3,i), ...
                    [' ', S.LocalALs(i).Point.Name, '  '], 'FontSize', 10, 'FontWeight', 'Bold', 'Color', 'm', ...
                    'HorizontalAlignment', 'Right');
                line([S.TCS_Pos_Landmark(1,i),0],[S.TCS_Pos_Landmark(2,i),0],[S.TCS_Pos_Landmark(3,i),0],'Color','m');
            end
            % Plot LCS axes
            quiver3(TCS_Pos_OrLCS(1),TCS_Pos_OrLCS(2),TCS_Pos_OrLCS(3),TCS_R_LCS(1,1),TCS_R_LCS(2,1),TCS_R_LCS(3,1),Scale,'r');
            quiver3(TCS_Pos_OrLCS(1),TCS_Pos_OrLCS(2),TCS_Pos_OrLCS(3),TCS_R_LCS(1,2),TCS_R_LCS(2,2),TCS_R_LCS(3,2),Scale,'g');
            quiver3(TCS_Pos_OrLCS(1),TCS_Pos_OrLCS(2),TCS_Pos_OrLCS(3),TCS_R_LCS(1,3),TCS_R_LCS(2,3),TCS_R_LCS(3,3),Scale,'b');
            % The title
            title(['The LCS of ',S.Name,' in TCS. Markers in blue; Landmarks in red'], 'FontWeight', 'Bold', 'FontSize', 10);
            
        end
        function drawPosture(S,Glob_R_Par,Glob_Pos_OrPar,Glob_R_Child,Glob_Pos_OrChild,TCSP_Pos_Markers,...
                              TCSC_PosMarkers,TCSP_Pos_AL,TCSC_Pos_AL,PName,Glob_Pos_PMarkers,Glob_Pos_CMarkers,Glo_Pos_AL,ParentMarkers)
            NPMarkers = size(ParentMarkers,1);
            NCMarkers = size(S.LocalMarkers,1);
            % figure settings
            figure('Name',['PosFrame ',PName,'-',S.Name],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
            hold on, axis equal, grid on
            % calculate figure scale
            MaxSmark = max(max(abs(Glob_Pos_CMarkers)));
            Scale = 0.3 * MaxSmark;

            % plot Global axes
            quiver3(0,0,0,Scale,0,0,0,'k');
            quiver3(0,0,0,0,Scale,0,0,'k');
            quiver3(0,0,0,0,0,Scale,0,'k');

            % plot axes names
            text(Scale,0,0,'X   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % X-axis
            text(0,Scale,0,'Y   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Y-axis
            text(0,0,Scale,'Z   ','FontSize',12,'FontWeight','Bold','Color','k','HorizontalAlignment','right'); % Z-axis

            % plot Global parent TCS axes
            quiver3(Glob_Pos_OrPar(1),Glob_Pos_OrPar(2),Glob_Pos_OrPar(3),...
                    Glob_R_Par(1,1),Glob_R_Par(2,1),Glob_R_Par(3,1),Scale,'r');
            quiver3(Glob_Pos_OrPar(1),Glob_Pos_OrPar(2),Glob_Pos_OrPar(3),...
                    Glob_R_Par(1,2),Glob_R_Par(2,2),Glob_R_Par(3,2),Scale,'g');
            quiver3(Glob_Pos_OrPar(1),Glob_Pos_OrPar(2),Glob_Pos_OrPar(3),...
                    Glob_R_Par(1,3),Glob_R_Par(2,3),Glob_R_Par(3,3),Scale,'b');

            % plot Global chil TCS axes
            quiver3(Glob_Pos_OrChild(1),Glob_Pos_OrChild(2),Glob_Pos_OrChild(3),...
                    Glob_R_Child(1,1),Glob_R_Child(2,1),Glob_R_Child(3,1),Scale,'m');
            quiver3(Glob_Pos_OrChild(1),Glob_Pos_OrChild(2),Glob_Pos_OrChild(3),...
                    Glob_R_Child(1,2),Glob_R_Child(2,2),Glob_R_Child(3,2),Scale,'y');
            quiver3(Glob_Pos_OrChild(1),Glob_Pos_OrChild(2),Glob_Pos_OrChild(3),...
                    Glob_R_Child(1,3),Glob_R_Child(2,3),Glob_R_Child(3,3),Scale,'c');

            % Plot parent markers in TCS
            for i=1:NPMarkers
                hold on;
                plot3(Glob_Pos_OrPar(1)+Glob_R_Par(1,:)*TCSP_Pos_Markers(:,i),...
                    Glob_Pos_OrPar(2)+Glob_R_Par(2,:)*TCSP_Pos_Markers(:,i),...
                    Glob_Pos_OrPar(3)+Glob_R_Par(3,:)*TCSP_Pos_Markers(:,i),'co')
                text( Glob_Pos_OrPar(1)+Glob_R_Par(1,:)*TCSP_Pos_Markers(:,i),...
                    Glob_Pos_OrPar(2)+Glob_R_Par(2,:)*TCSP_Pos_Markers(:,i),...
                    Glob_Pos_OrPar(3)+Glob_R_Par(3,:)*TCSP_Pos_Markers(:,i),...
                    ['  ', ParentMarkers(i).Point.Name], 'FontSize', 10, 'Color', 'c');
                line([Glob_Pos_OrPar(1)+Glob_R_Par(1,:)*TCSP_Pos_Markers(:,i),Glob_Pos_OrPar(1)],...
                    [Glob_Pos_OrPar(2)+Glob_R_Par(2,:)*TCSP_Pos_Markers(:,i),Glob_Pos_OrPar(2)],...
                    [Glob_Pos_OrPar(3)+Glob_R_Par(3,:)*TCSP_Pos_Markers(:,i),Glob_Pos_OrPar(3)],...
                    'Color','c');
            end
            % Plot parent markers in Globals
            for i=1:NPMarkers
                hold on;
                plot3(Glob_Pos_PMarkers(1,i), Glob_Pos_PMarkers(2,i), Glob_Pos_PMarkers(3,i), 'cx');
                text(Glob_Pos_PMarkers(1,i), Glob_Pos_PMarkers(2,i), Glob_Pos_PMarkers(3,i), ...
                    ['  ', ParentMarkers(i).Point.Name], 'FontSize', 10, 'Color', 'c');
                line([Glob_Pos_PMarkers(1,i),0],[Glob_Pos_PMarkers(2,i),0],[Glob_Pos_PMarkers(3,i),0],'Color','c');
            end
            % Plot child markers in TCS
            for i=1:NCMarkers
                hold on;
                plot3(Glob_Pos_OrChild(1)+Glob_R_Child(1,:)*TCSC_PosMarkers(:,i),...
                    Glob_Pos_OrChild(2)+Glob_R_Child(2,:)*TCSC_PosMarkers(:,i),...
                    Glob_Pos_OrChild(3)+Glob_R_Child(3,:)*TCSC_PosMarkers(:,i),'bo')
                text( Glob_Pos_OrChild(1)+Glob_R_Child(1,:)*TCSC_PosMarkers(:,i),...
                    Glob_Pos_OrChild(2)+Glob_R_Child(2,:)*TCSC_PosMarkers(:,i),...
                    Glob_Pos_OrChild(3)+Glob_R_Child(3,:)*TCSC_PosMarkers(:,i),...
                    ['  ', S.LocalMarkers(i).Point.Name], 'FontSize', 10, 'Color', 'b');
                line([Glob_Pos_OrChild(1)+Glob_R_Child(1,:)*TCSC_PosMarkers(:,i),Glob_Pos_OrChild(1)],...
                    [Glob_Pos_OrChild(2)+Glob_R_Child(2,:)*TCSC_PosMarkers(:,i),Glob_Pos_OrChild(2)],...
                    [Glob_Pos_OrChild(3)+Glob_R_Child(3,:)*TCSC_PosMarkers(:,i),Glob_Pos_OrChild(3)],...
                    'Color','b');
            end
            % Plot child markers in Globals
            for i=1:NCMarkers
                hold on;
                plot3(Glob_Pos_CMarkers(1,i), Glob_Pos_CMarkers(2,i), Glob_Pos_CMarkers(3,i), 'bx');
                text(Glob_Pos_CMarkers(1,i), Glob_Pos_CMarkers(2,i), Glob_Pos_CMarkers(3,i), ...
                    ['  ', S.LocalMarkers(i).Point.Name], 'FontSize', 10, 'Color', 'c');
                line([Glob_Pos_CMarkers(1,i),0],[Glob_Pos_CMarkers(2,i),0],[Glob_Pos_CMarkers(3,i),0],'Color','b');
            end
            % Plot AL in Global
            hold on;
            plot3(Glo_Pos_AL(1),Glo_Pos_AL(2),Glo_Pos_AL(3),'rx');
            line([Glo_Pos_AL(1),0],[Glo_Pos_AL(2),0],[Glo_Pos_AL(3),0],'Color','r');
            text(Glo_Pos_AL(1),Glo_Pos_AL(2),Glo_Pos_AL(3), ...
                ['  ', S.LocalPoints(1).Point.Name], 'FontSize', 10, 'Color', 'r');

            % Plot AL in parent TCS
            hold on;
            plot3(Glob_Pos_OrPar(1)+Glob_R_Par(1,:)*TCSP_Pos_AL,...
                Glob_Pos_OrPar(2)+Glob_R_Par(2,:)*TCSP_Pos_AL,...
                Glob_Pos_OrPar(3)+Glob_R_Par(3,:)*TCSP_Pos_AL,'go')
            text( Glob_Pos_OrPar(1)+Glob_R_Par(1,:)*TCSP_Pos_AL,...
                Glob_Pos_OrPar(2)+Glob_R_Par(2,:)*TCSP_Pos_AL,...
                Glob_Pos_OrPar(3)+Glob_R_Par(3,:)*TCSP_Pos_AL,...
                ['  ', S.LocalPoints(1).Point.Name], 'FontSize', 10, 'Color', 'g');
            line([Glob_Pos_OrPar(1)+Glob_R_Par(1,:)*TCSP_Pos_AL,Glob_Pos_OrPar(1)],...
                [Glob_Pos_OrPar(2)+Glob_R_Par(2,:)*TCSP_Pos_AL,Glob_Pos_OrPar(2)],...
                [Glob_Pos_OrPar(3)+Glob_R_Par(3,:)*TCSP_Pos_AL,Glob_Pos_OrPar(3)],...
                'Color','k');

            % Plot AL in child TCS
            hold on;
            plot3(Glob_Pos_OrChild(1)+Glob_R_Child(1,:)*TCSC_Pos_AL,...
                Glob_Pos_OrChild(2)+Glob_R_Child(2,:)*TCSC_Pos_AL,...
                Glob_Pos_OrChild(3)+Glob_R_Child(3,:)*TCSC_Pos_AL,'bo')
            text( Glob_Pos_OrChild(1)+Glob_R_Child(1,:)*TCSC_Pos_AL,...
                Glob_Pos_OrChild(2)+Glob_R_Child(2,:)*TCSC_Pos_AL,...
                Glob_Pos_OrChild(3)+Glob_R_Child(3,:)*TCSC_Pos_AL,...
                ['  ', S.LocalPoints(1).Point.Name], 'FontSize', 10, 'Color', 'b');
            line([Glob_Pos_OrChild(1)+Glob_R_Child(1,:)*TCSC_Pos_AL,Glob_Pos_OrChild(1)],...
                [Glob_Pos_OrChild(2)+Glob_R_Child(2,:)*TCSC_Pos_AL,Glob_Pos_OrChild(2)],...
                [Glob_Pos_OrChild(3)+Glob_R_Child(3,:)*TCSC_Pos_AL,Glob_Pos_OrChild(3)],...
                'Color','k');
            if strcmpi(S.Name,'LeftShank')
                view([78 8]);
            elseif strcmpi(S.Name,'LeftFoot')
                view([64 20]);
            end
        end
        function [Glb_Pos_Head, Loc_Pos_Head, Rad_Head, Rot_Lc_2_Gl_CS_Ori] = calcKJCinThighLCS_ULB(S,AL,B_S000_Left)
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%%  mcc -W lib:Ili_ACC_Regr_4_Als -T link:lib Ili_ACC_Regr_4_Als
             %%%  mcc -W cpplib:Ili_ACC_Regr_4_Als -T link:lib Ili_ACC_Regr_4_Als
             %%%  mcc -m Ili_ACC_Regr_4_Als
             %%%  mcc -W cpplib:LwLmb_Regr_DLL_Als -T link:lib Ili_ACC_Regr_4_Als Fem_HCND_Regr_3_Als
             %%%  mcc -W cpplib:Ili_ACC_Regr_4_Als -T link:lib Ili_ACC_Regr_4_Als Fem_HCND_Regr_3_Als
             
             if B_S000_Left == 1, AL(:,1) = -AL(:,1); end %%% Mirroring
             %%%  -------------- Head Position
             %%%  -------------- Cond Sulcus Position
             %%%  -------------- Lat Cond_Ell Position
             %%%  -------------- Med Cond_Ell Position
             Rgr_Pos_Head = [
                 1.01137      2.93269      1.96936      13.4606      14.6599     -29.7586
                 -0.0172652     0.201561     0.382503      11.4308      10.5479     -7.37479
                 -0.0067019    -0.31658     0.16603      7.9518      1.5578     0.73851
                 -0.014986     0.39484    0.017388      1.8079      3.1942     -5.2105
                 ];
             %%%  -------------- Head Radius
             %%%  -------------- Cond Sulcus Radius
             %%%  -------------- Lat Cond_Ell_X Semi_ax
             %%%  -------------- Lat Cond_Ell_Y Semi_ax
             %%%  -------------- Lat Cond_Ell_Z Semi_ax
             %%%  -------------- Med Cond_Ell_X Semi_ax
             %%%  -------------- Med Cond_Ell_Y Semi_ax
             %%%  -------------- Med Cond_Ell_Z Semi_ax
             Rgr_Rad_Head = [
                 0.011227    0.019676     0.33298      4.4669
                 0.00143702      1.16671    -0.716999      2.68848
                 0.046823      1.6553     -1.4928     -2.4381
                 0.010953     -0.8989      1.2417    -0.97728
                 0.036127    0.056774     0.37651      3.0915
                 0.0211177     0.191688   -0.0925776      11.4015
                 -0.0085972      1.2301    -0.79821      2.9883
                 0.035049     -1.0983      1.4625      3.7847
                 ];
             Ind_Regr_ALs = [1 2 3];
             
             Nmb_ALs = size(AL,1);
             P_1 = AL(1,:)'; %%% FTC (GT)
             P_2 = AL(2,:)'; %%% FME
             P_3 = AL(3,:)'; %%% FLE
             CS_Ori = (P_2+P_3)/2;
             %%%% IT WORKS !!!!!!!!!!!!!! Correction to ESCAPE Regr_M SINGULARITY
             CS_Ori = CS_Ori + [.1; .1; .1]; %%%% IT WORKS !!!!!!!!!!!!!!
             %%%% IT WORKS !!!!!!!!!!!!!! Correction to ESCAPE Regr_M SINGULARITY
             y_CS = (P_1-CS_Ori)/norm(P_1-CS_Ori);
             V_tmp = cross(P_2-CS_Ori,y_CS);
             x_CS = V_tmp/norm(V_tmp);
             z_CS = cross(x_CS,y_CS);
             Rot_Lc_2_Gl = [x_CS y_CS z_CS];
             Rot_Lc_2_Gl_CS_Ori = [Rot_Lc_2_Gl CS_Ori];
             %%% disp('num2str(Rot_Lc_2_Gl_CS_Ori)')
             %%% disp(num2str(Rot_Lc_2_Gl_CS_Ori))
             
             %%%%%%%%%%%%%%%%%%% Local CS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             CS_Ori_rot = CS_Ori'*Rot_Lc_2_Gl;
             Alum_2_CT1(:,1) =  AL*Rot_Lc_2_Gl(:,1)-CS_Ori_rot(1);
             Alum_2_CT1(:,2) =  AL*Rot_Lc_2_Gl(:,2)-CS_Ori_rot(2);
             Alum_2_CT1(:,3) =  AL*Rot_Lc_2_Gl(:,3)-CS_Ori_rot(3);
             %%% disp('Alum_2_CT1')
             %%% disp(num2str(Alum_2_CT1))
             Mult_0_95 = [1 1 1];
             if 0
                 %%% Mult_0_95 = [.75 .95 1];
                 %%% Mult_0_95 = [.6 1 1];
                 %%% ???
                 Mult_0_95 = [1 1 .8];
                 if 0
                     m_Alum_2_CT1 = mean(Alum_2_CT1);
                     C_Alum_2_CT1 = Alum_2_CT1 - repmat(m_Alum_2_CT1,Nmb_ALs,1);
                     C_Alum_2_CT1 = repmat(Mult_0_95,Nmb_ALs,1).*C_Alum_2_CT1;
                     Alum_2_CT1 = C_Alum_2_CT1 + repmat(m_Alum_2_CT1,Nmb_ALs,1);
                 else
                     Alum_2_CT1 = repmat(Mult_0_95,Nmb_ALs,1).*Alum_2_CT1;
                 end
                 %%% disp([num2str(Mult_0_95) '  Alum_2_CT1'])
                 %%% disp(num2str(Alum_2_CT1))
             end %%% if 1
             N_Pos_Coeffs = size(Rgr_Pos_Head,1);
             for k = 1:N_Pos_Coeffs
                 for i = 1:3 %%% XYZ by regr
                     Loc_Pos_Head(k,i) = Rgr_Pos_Head(k,1:Nmb_ALs)*Alum_2_CT1(Ind_Regr_ALs,i)+Rgr_Pos_Head(k,Nmb_ALs+i);
                     Loc_Pos_Head(k,i) = Loc_Pos_Head(k,i)/Mult_0_95(i);
                 end %%% for i = 1:3 %%% XYZ by regr
                 %%%%%%%%%%%%%%%%%%% Local CS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %%% disp(['Loc_Pos_Head XYZ, Mult_0_95 => ' num2str([Loc_Pos_Head(k,:) Mult_0_95])])
                 Glb_Pos_Head(k,:) = (Rot_Lc_2_Gl_CS_Ori*[Loc_Pos_Head(k,:)'; 1])';
                 if B_S000_Left == 1, Glb_Pos_Head(k,1) = -Glb_Pos_Head(k,1);   end %%% Back Mirroring
                 %%% disp(['Glb_Pos_Head XYZ => ' num2str(Glb_Pos_Head(k,:))])
             end %%% for k = 1:N_Pos_Coeffs
             N_Rad_Coeffs = size(Rgr_Rad_Head,1);
             for k = 1:N_Rad_Coeffs
                 Rad_Head(k) = Rgr_Rad_Head(k,:) * [sqrt(sum(Alum_2_CT1(Ind_Regr_ALs,:).^2,2)); 1];
                 %%% disp(['Rad_Head = ' num2str(Rad_Head(k))])
             end %%% for k = 1:N_Rad_Coeffs
             %%% function [Glb_Pos_Head, Loc_Pos_Head, Rad_Head, Rot_Lc_2_Gl_CS_Ori] = Fem_HCND_Regr_3_Als(AL, B_S000_Left)
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
        function getLocalPointsCoords(S)
            % extract transformation data from Body_T
            Glob_Pos_OrBody = S.Glob_Pos_OrSeg;
            Glob_R_Loc      = S.Glob_R_Seg;
            
            % Equations for changing coordinates
            % globP = Glob_Pos_OrBody + Glob_R_Loc * locP
            % locP  = Glob_R_Loc' * (globP - Glob_Pos_OrBody)
            
            % Get the local coordinates to all points, markers and CoM
            nPoints  = size(S.LocalPoints,1);
            nMarkers = size(S.LocalMarkers,1);
            for i=1:nPoints
                S.LocalPoints(i).LocCoord =  Glob_R_Loc' * (S.LocalPoints(i).Point.GlobalCoord - Glob_Pos_OrBody);
            end
            for i = 1:nMarkers
                S.LocalMarkers(i).LocCoord = Glob_R_Loc' * (S.LocalMarkers(i).Point.GlobalCoord - Glob_Pos_OrBody);
            end
            S.CoM.LocCoord = Glob_R_Loc' * (S.CoM.Point.GlobalCoord - Glob_Pos_OrBody);
        end
        function [R,Or] = getRd(S,q_i)
            if S.Fixed == 1 % If is Ground                                                    
            %if strcmp(S.Name,'Ground')
                Pini_loc  = [0;0;0];
                Pini_glob = [0;0;0];
                Rloc  = [1,0,0;0,1,0;0,0,1];
                Rglob = Rloc;
                Rinv_loc = inv(Rloc);
            
            else
                % local and global coord of PIni
                Pini_loc = S.LocalPoints(1).LocCoord;
                PosInqPIni = S.LocalPoints(1).Point.PosInq;
                if S.LocalPoints(1).Point.Fixed == 1
                    Pini_glob(1,1) = S.LocalPoints(1).Point.GlobalCoord(1);
                    Pini_glob(2,1) = S.LocalPoints(1).Point.GlobalCoord(2);
                    Pini_glob(3,1) = S.LocalPoints(1).Point.GlobalCoord(3);
                elseif isnan(q_i(PosInqPIni))
                    % This means the segment posture could not be reconstructed
                    % Then the segment is placed at the origin and scaled with a very small factor.
                    % Then the segment "dissappears" 
                    R = 0.0001*eye(3);
                    Or = [0; 0; 0];
                    return
                else
                    Pini_glob(1,1) = q_i(PosInqPIni);
                    Pini_glob(2,1) = q_i(PosInqPIni+1);
                    Pini_glob(3,1) = q_i(PosInqPIni+2);
                end
                % local and global coord of basis of Segment
                Rloc_u = S.LocalVectors(1).LocCoord;
                Rloc_v = S.LocalVectors(2).LocCoord;
                if S.LocalVectors(1).Vector.Fixed == 1
                    Rglob_u = Rloc_u;
                else
                   PosInqU = S.LocalVectors(1).Vector.PosInq;
                   Rglob_u(1,1) = q_i(PosInqU); Rglob_u(2,1) = q_i(PosInqU+1); Rglob_u(3,1) = q_i(PosInqU+2);
                   %Rglob_u(:,1) = q_i(PosInqU:PosInqU+2); 
                end
                if S.LocalVectors(2).Vector.Fixed == 1
                   Rglob_v = Rloc_v;
                else
                    PosInqV = S.LocalVectors(2).Vector.PosInq;
                    Rglob_v(1,1) = q_i(PosInqV); Rglob_v(2,1) = q_i(PosInqV+1); Rglob_v(3,1) = q_i(PosInqV+2);
                end
                if (S.BasisType == 1)
                    Rloc_w = S.LocalPoints(2).LocCoord - S.LocalPoints(1).LocCoord;
                    PosInqP2 = S.LocalPoints(2).Point.PosInq;
                    PosInqP1 = S.LocalPoints(1).Point.PosInq;
                    if S.LocalPoints(1).Point.Fixed == 1
                        Rglob_w(1,1) = q_i(PosInqP2)   - S.LocalPoints(1).Point.GlobalCoord(1);
                        Rglob_w(2,1) = q_i(PosInqP2+1) - S.LocalPoints(1).Point.GlobalCoord(2);
                        Rglob_w(3,1) = q_i(PosInqP2+2) - S.LocalPoints(1).Point.GlobalCoord(3);
                    else
                        Rglob_w(1,1) = q_i(PosInqP2)   - q_i(PosInqP1);
                        Rglob_w(2,1) = q_i(PosInqP2+1) - q_i(PosInqP1+1);
                        Rglob_w(3,1) = q_i(PosInqP2+2) - q_i(PosInqP1+2);
                    end
                elseif (S.BasisType == 0)
                    % Original Alex Valero
                    %Rloc_w = S.LocalVectors(3).LocCoord;
                    %PosInqW = S.LocalPoints(3).Point.PosInq;
                    %Rglob_w(1,1) = q_i(PosInqW); Rglob_v(2,1) = q_i(PosInqW+1); Rglob_v(3,1) = q_i(PosInqW+2);
                    % NEW VERSION. SAUSEJO
                    Rloc_w = S.LocalVectors(3).LocCoord;
                    if S.LocalVectors(3).Vector.Fixed == 1
                        Rglob_w = Rloc_w;
                    else
                        PosInqW = S.LocalVectors(3).Vector.PosInq;
                        Rglob_w(1,1) = q_i(PosInqW); Rglob_w(2,1) = q_i(PosInqW+1); Rglob_w(3,1) = q_i(PosInqW+2);
                    end
                    
                else
                    error('Incorrect BasisType in Segment %d',S.Name);
                end
                % normalized rotation matrix
                Rloc(:,1) = [Rloc_u(1)/norm(Rloc_u), Rloc_u(2)/norm(Rloc_u),Rloc_u(3)/norm(Rloc_u)];
                Rloc(:,2) = [Rloc_v(1)/norm(Rloc_v), Rloc_v(2)/norm(Rloc_v),Rloc_v(3)/norm(Rloc_v)];
                Rloc(:,3) = [Rloc_w(1)/norm(Rloc_w), Rloc_w(2)/norm(Rloc_w),Rloc_w(3)/norm(Rloc_w)];
                Rglob(:,1) = [Rglob_u(1)/norm(Rglob_u), Rglob_u(2)/norm(Rglob_u),Rglob_u(3)/norm(Rglob_u)];
                Rglob(:,2) = [Rglob_v(1)/norm(Rglob_v), Rglob_v(2)/norm(Rglob_v),Rglob_v(3)/norm(Rglob_v)];
                Rglob(:,3) = [Rglob_w(1)/norm(Rglob_w), Rglob_w(2)/norm(Rglob_w),Rglob_w(3)/norm(Rglob_w)];
                Rinv_loc = inv(Rloc);
            end
            % rotation matrix for segment
            R = Rglob * Rinv_loc;
            % origin in glob for segment
            Or = Pini_glob - R * Pini_loc;
        end
        function [Glb_Pos_Head, Loc_Pos_Head, Rad_Head, Rot_Lc_2_Gl_CS_Ori] = ...
                 calcHJCsInPelvisLCS_ULB(S,ALS_1245, Pelvis_2_Iliacs_Yes, B_S000_Left)%#ok<MANU>
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%  mcc -W lib:Ili_ACC_Regr_4_Als -T link:lib Ili_ACC_Regr_4_Als
            %%%  mcc -W cpplib:Ili_ACC_Regr_4_Als -T link:lib Ili_ACC_Regr_4_Als
            %%%  mcc -m Ili_ACC_Regr_4_Als
            %%%  mcc -W cpplib:LwLmb_Regr_DLL_Als -T link:lib Ili_ACC_Regr_4_Als Fem_HCND_Regr_3_Als
            %%%  mcc -W cpplib:Ili_ACC_Regr_4_Als -T link:lib Ili_ACC_Regr_4_Als Fem_HCND_Regr_3_Als
            
            if B_S000_Left == 1, ALS_1245(:,1) = -ALS_1245(:,1); end %%% Mirroring
            if  Pelvis_2_Iliacs_Yes == 1
                %%% 12 pelv (S026) 3ALs LIAS;LIPS;RIPS;
                Rgr_Pos_Head = [ 0.0142699565539886     0.104783654354478    0.0997309860279856     -12.8655957833115     -74.6798805256388      89.1529356127467];
                Rgr_Rad_Head = [ 0.0233550822068977    0.0599730210345513     0.107208643636248     -1.81778288982969];
            else %%% if  Pelvis_2_Iliacs_Yes == 1
                %%%%%%%%%%% S026 127 25oct09 Bone Corrected S026_L
                Rgr_Pos_Head = [  0.283452158527851    -0.548879148471912    0.0249463823865055     -2.06396900838662      15.7551115621845      19.8529047076718];
                Rgr_Rad_Head = [ 0.0818007669994074    0.0531096304550829    0.0501625165948081      2.04956753878832];
            end %%% if  Pelvis_2_Iliacs_Yes == 1
            Ind_Regr_ALs = [1 2 4];
            Nmb_ALs = size(ALS_1245,1);
            if  Pelvis_2_Iliacs_Yes == 1
                P_1 = ALS_1245(1,:)'; %%%  1_LIAS
                P_2 = ALS_1245(2,:)'; %%%  2_LIPS
                P_3 = ALS_1245(3,:)'; %%%  3_RIAS
                P_4 = ALS_1245(4,:)'; %%%  4_RIPS
                P_5 = (P_2+P_4)/2; %%%
                V_tmp = P_3-P_1;
                z_CS = V_tmp/norm(V_tmp);
                CS_Ori = (P_1+P_3)/2; %%% Local CS Ori
                V_tmp = cross(P_5-CS_Ori,z_CS);
                y_CS = V_tmp/norm(V_tmp);
                x_CS = cross(y_CS,z_CS);
            else %%% if  Pelvis_2_Iliacs_Yes == 1
                %%%%%%%%  A.	point 1: IAS 1
                %%%%%%%%  B.	point 2: IPS 2
                %%%%%%%%  C.	point 3: IPI NO
                %%%%%%%%  D.	point 4: IIT 3
                %%%%%%%%  E.	point 5: IPJ NO
                %%%%%%%%  F.	point 6: IPP 4 or NO
                %%%%%%%%  G.	point 7: ICT 4 or NO
                %%% point 1: IAS 1
                %%% point 2: IPS 2
                %%% point 4: IIT 3
                P_1 = ALS_1245(1,:)'; %%%  1_IAS_ili
                P_2 = ALS_1245(2,:)'; %%%  2_IPS_ili
                P_3 = ALS_1245(3,:)'; %%%  4_ITT_ili
                CS_Ori = (P_2+P_3)/2; %%% Local CS Ori
                V_tmp = cross(P_2-P_3,P_1-P_3);
                x_CS = V_tmp/norm(V_tmp);
                V_tmp = (P_1+P_2)/2-P_3;
                y_CS = V_tmp/norm(V_tmp);
                z_CS = cross(x_CS,y_CS);
            end %%% if  Pelvis_2_Iliacs_Yes == 1
            Rot_Lc_2_Gl = [x_CS y_CS z_CS]; %%% Local CS Rot_M
            %%% Gl = Rot_Lc_2_Gl*Loc+CS_Ori %%% Loc to Gl
            Rot_Lc_2_Gl_CS_Ori = [Rot_Lc_2_Gl CS_Ori];
            %%% disp('num2str(Rot_Lc_2_Gl_CS_Ori)')
            %%% disp(num2str(Rot_Lc_2_Gl_CS_Ori))
            
            %%%%%%%%%%%%%%%%%%% Local CS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            CS_Ori_rot = CS_Ori'*Rot_Lc_2_Gl;
            Alum_2_CT1(:,1) =  ALS_1245*Rot_Lc_2_Gl(:,1)-CS_Ori_rot(1);
            Alum_2_CT1(:,2) =  ALS_1245*Rot_Lc_2_Gl(:,2)-CS_Ori_rot(2);
            Alum_2_CT1(:,3) =  ALS_1245*Rot_Lc_2_Gl(:,3)-CS_Ori_rot(3);
            %%% disp('Alum_2_CT1')
            %%% disp(num2str(Alum_2_CT1))
            Mult_0_95 = [1 1 1];
            if 0
                %%% Mult_0_95 = [.75 .95 1];
                %%% Mult_0_95 = [.6 1 1];
                %%% ???
                Mult_0_95 = [1 1 .8];
                if 0
                    m_Alum_2_CT1 = mean(Alum_2_CT1);
                    C_Alum_2_CT1 = Alum_2_CT1 - repmat(m_Alum_2_CT1,Nmb_ALs,1);
                    C_Alum_2_CT1 = repmat(Mult_0_95,Nmb_ALs,1).*C_Alum_2_CT1;
                    Alum_2_CT1 = C_Alum_2_CT1 + repmat(m_Alum_2_CT1,Nmb_ALs,1);
                else
                    Alum_2_CT1 = repmat(Mult_0_95,Nmb_ALs,1).*Alum_2_CT1;
                end
                %%% disp([num2str(Mult_0_95) '  Alum_2_CT1'])
                %%% disp(num2str(Alum_2_CT1))
            end %%% if 1
            %%% point 1: IAS 1
            %%% point 2: IPS 2
            %%% point 7: ICT 4
            for i = 1:3 %%% XYZ by regr
                Loc_Pos_Head(i) = Rgr_Pos_Head(1:Nmb_ALs-1)*Alum_2_CT1(Ind_Regr_ALs,i)+Rgr_Pos_Head(Nmb_ALs-1+i);
                Loc_Pos_Head(i) = Loc_Pos_Head(i)/Mult_0_95(i);
            end %%% for i = 1:3 %%% XYZ by regr
            %%% disp(['Loc_Pos_Head XYZ, Mult_0_95 => ' num2str([Loc_Pos_Head Mult_0_95])])
            Glb_Pos_Head = Rot_Lc_2_Gl_CS_Ori*[Loc_Pos_Head'; 1];
            if B_S000_Left == 1, Glb_Pos_Head(1) = -Glb_Pos_Head(1);   end %%% Back Mirroring
            %%% disp(['Glb_Pos_Head XYZ => ' num2str(Glb_Pos_Head')])
            Rad_Head = Rgr_Rad_Head * [sqrt(sum(Alum_2_CT1(Ind_Regr_ALs,:).^2,2)); 1];
            %%% disp(['Rad_Head = ' num2str(Rad_Head)])
%             return %%% function [Glb_Pos_Head, Loc_Pos_Head, Rad_Head, Rot_Lc_2_Gl_CS_Ori] = Ili_ACC_Regr_4_Als(ALS_1245, Pelvis_2_Iliacs_Yes, B_S000_Left)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
        function Phi = mkAdditionalCtrs(S)
            NLocalVectorsAdded = size(S.LocalVectorsAdded,1);
            Phi = {};
            for i=1:NLocalVectorsAdded
                V = S.LocalVectorsAdded(i).Vector.CoordName;
                locV = S.LocalVectorsAdded(i).LocCoord;
                Phi = [Phi;S.mkCtrCharLinComb(V,locV)];
            end
        end
        function mkCharSegRef(S)
            %MKCHARSEGREF estimated char axes depending of type of basis
           
            Axis1 = S.LocalVectors(1).Vector.CoordName;
            Axis2 = S.LocalVectors(2).Vector.CoordName;

            % Check type of basis
            % Basis formed by 3 vectors
            if S.BasisType == 0
                Axis3 = S.LocalVectors(3).Vector.CoordName;
            else
                P1 = S.LocalPoints(1).Point.CoordName;
                P2 = S.LocalPoints(2).Point.CoordName;
                minus = ' - ';
                Axis3 = {[P2{1},minus,P1{1}];[P2{2},minus,P1{2}];[P2{3},minus,P1{3}]};
            end
            S.SegRef.Dir1.Global = Axis1;
            S.SegRef.Dir2.Global = Axis2;
            S.SegRef.Dir3.Global = Axis3;
            
        end
        function Phi = mkCtrCharLinComb(S,U,locU)
            %MKCTRCHARLINCOMB Summary of this function goes here
            %   Detailed explanation goes here
            if isnumeric(locU)
                locUtmp{1} = num2str(locU(1));
                locUtmp{2} = num2str(locU(2));
                locUtmp{3} = num2str(locU(3));
                locU = [];
                locU = locUtmp;
            end
            if isempty (S.InvR{1,1}) || strcmpi(locU{1},'0') A1 =[]; else A1 = ['(',S.InvR{1,1},')*(',locU{1},')']; end
            if isempty (S.InvR{1,2}) || strcmpi(locU{2},'0') B1 =[]; else B1 = ['+((',S.InvR{1,2},')*(',locU{2},'))']; end
            if isempty (S.InvR{1,3}) || strcmpi(locU{3},'0') C1 =[]; else C1 = ['+((',S.InvR{1,3},')*(',locU{3},'))']; end
            if isempty (S.InvR{2,1}) || strcmpi(locU{1},'0') A2 =[]; else A2 = ['(',S.InvR{2,1},')*(',locU{1},')']; end
            if isempty (S.InvR{2,2}) || strcmpi(locU{2},'0') B2 =[]; else B2 = ['+((',S.InvR{2,2},')*(',locU{2},'))']; end
            if isempty (S.InvR{2,3}) || strcmpi(locU{3},'0') C2 =[]; else C2 = ['+((',S.InvR{2,3},')*(',locU{3},'))']; end
            if isempty (S.InvR{3,1}) || strcmpi(locU{1},'0') A3 =[]; else A3 = ['(',S.InvR{3,1},')*(',locU{1},')']; end
            if isempty (S.InvR{3,2}) || strcmpi(locU{2},'0') B3 =[]; else B3 = ['+((',S.InvR{3,2},')*(',locU{2},'))']; end
            if isempty (S.InvR{3,3}) || strcmpi(locU{3},'0') C3 =[]; else C3 = ['+((',S.InvR{3,3},')*(',locU{3},'))']; end

            % Calc coefcients
            coef{1} = [A1,B1,C1];
            coef{2} = [A2,B2,C2];
            coef{3} = [A3,B3,C3];
            
            B {1} = ['- (',coef{1},')* (',S.SegRef.Dir1.Global{1},')'];
            B {2} = ['- (',coef{1},')* (',S.SegRef.Dir1.Global{2},')'];
            B {3} = ['- (',coef{1},')* (',S.SegRef.Dir1.Global{3},')'];
            C {1} = ['- (',coef{2},')* (',S.SegRef.Dir2.Global{1},')'];
            C {2} = ['- (',coef{2},')* (',S.SegRef.Dir2.Global{2},')'];
            C {3} = ['- (',coef{2},')* (',S.SegRef.Dir2.Global{3},')'];
            D {1} = ['- (',coef{3},')* (',S.SegRef.Dir3.Global{1},')'];
            D {2} = ['- (',coef{3},')* (',S.SegRef.Dir3.Global{2},')'];
            D {3} = ['- (',coef{3},')* (',S.SegRef.Dir3.Global{3},')'];
            
            if isempty(coef{1})  B{1} = []; B{2} = []; B{3} = []; end
            if isempty(coef{2})  C{1} = []; C{2} = []; C{3} = []; end
            if isempty(coef{3})  D{1} = []; D{2} = []; D{3} = []; end
            
            Phi = {[U{1},B{1},C{1},D{1}];[U{2},B{2},C{2},D{2}];[U{3},B{3},C{3},D{3}]};
                                 
        end
        function Phi = mkRelPointCtr(S,SegRef,InvR)
            PointRelGlob = mksymbolicXYZ(S.LocalPointsRel(1).Point.Name(1:end-2));
            PointRel = mkSymbolicXYZ(S.LocalPointsRel(1).Point.Name);
            PIni = mkSymbolicXYZ(S.LocalPoints(1).Point.Name);
            PLocIni  = mkSymbolicXYZ([S.Name,'_',S.LocalPoints(1).Point.Name]);
            LocU = PointRel - PLocIni;
            U    = PointRelGlob - PIni;
%             Phi = mkCtrLinComb(U, LocU, SegRef,InvR);
            Phi = PointRelGlob - PIni - [SegRef.Dir1.Global,SegRef.Dir2.Global,SegRef.Dir3.Global]*PointRel;
        end
        function Phi = mkCtrRigidSegment(S)
            SegName = S.Name;
            SegRef = mkSymbSegRef(S.BasisType,S.LocalPoints,S.LocalVectors,SegName);
            Phi = mkCtrAxisRef3D(SegRef);
            % projection of vector U on each axis of the reference.
            R = [SegRef.Dir1.Local, SegRef.Dir2.Local, SegRef.Dir3.Local];
            InvR = inv(R);
            if ~isnumeric(InvR)
                S.sym2CharInvR(InvR);
            else
                S.num2charInvR(InvR)
            end
            if ~isempty(S.LocalPointsRel)
                Phi =[Phi; S.mkRelPointCtr(SegRef,InvR)];
            end
            Phi = [Phi; mkCtrPointsNotInRef(S.BasisType,S.LocalPoints,SegRef,SegName,InvR)];
            if ~isempty(S.LocalMarkers)
                Phi = [Phi; mkCtrMarkersNotInRef(S.BasisType,S.LocalPoints,S.LocalMarkers,SegRef,SegName,InvR)];
            end
            Phi = [Phi; mkCtrVectorsNotInRef(S.BasisType,S.LocalVectors,SegRef,SegName,InvR)];
        end
        function num2charInvR(S,InvR)
            if InvR(1,1)==0 S.InvR {1,1}= []; else S.InvR {1,1} = num2str(InvR(1,1)); end
            if InvR(1,2)==0 S.InvR {1,2}= []; else S.InvR {1,2} = num2str(InvR(1,2)); end
            if InvR(1,3)==0 S.InvR {1,3}= []; else S.InvR {1,3} = num2str(InvR(1,3)); end
            if InvR(2,1)==0 S.InvR {2,1}= []; else S.InvR {2,1} = num2str(InvR(2,1)); end
            if InvR(2,2)==0 S.InvR {2,2}= []; else S.InvR {2,2} = num2str(InvR(2,2)); end
            if InvR(2,3)==0 S.InvR {2,3}= []; else S.InvR {2,3} = num2str(InvR(2,3)); end
            if InvR(3,1)==0 S.InvR {3,1}= []; else S.InvR {3,1} = num2str(InvR(3,1)); end
            if InvR(3,2)==0 S.InvR {3,2}= []; else S.InvR {3,2} = num2str(InvR(3,2)); end
            if InvR(3,3)==0 S.InvR {3,3}= []; else S.InvR {3,3} = num2str(InvR(3,3)); end
        end
        function sym2CharInvR(S,InvR)
            InvR_tmp = char(InvR);
            InvR_tmp2 = InvR_tmp(10:end-1);  % Remove 'matrix([[' in first and ')' at the end.
            InvRChar = strread(InvR_tmp2, '%s', 'delimiter', ',[]'); % Return vector with InvR values and two empty values inter rows
            if strcmpi(InvRChar{1},'0') S.InvR {1,1}= []; else S.InvR {1,1} = InvRChar{1}; end
            if strcmpi(InvRChar{2},'0') S.InvR {1,2}= []; else S.InvR {1,2} = InvRChar{2}; end
            if strcmpi(InvRChar{3},'0') S.InvR {1,3}= []; else S.InvR {1,3} = InvRChar{3}; end
            if strcmpi(InvRChar{6},'0') S.InvR {2,1}= []; else S.InvR {2,1} = InvRChar{6}; end
            if strcmpi(InvRChar{7},'0') S.InvR {2,2}= []; else S.InvR {2,2} = InvRChar{7}; end
            if strcmpi(InvRChar{8},'0') S.InvR {2,3}= []; else S.InvR {2,3} = InvRChar{8}; end
            if strcmpi(InvRChar{11},'0') S.InvR {3,1}= []; else S.InvR {3,1} = InvRChar{11}; end
            if strcmpi(InvRChar{12},'0') S.InvR {3,2}= []; else S.InvR {3,2} = InvRChar{12}; end
            if strcmpi(InvRChar{13},'0') S.InvR {3,3}= []; else S.InvR {3,3} = InvRChar{13}; end
            
%             
%             S.InvR {1,1} = if strcmp(InvRChar{1},'0')InvRChar{1};  S.InvR {1,2} = InvRChar{2};  S.InvR {1,3} = InvRChar{3};
%             S.InvR {2,1} = InvRChar{6};  S.InvR {2,2} = InvRChar{7};  S.InvR {2,3} = InvRChar{8};
%             S.InvR {3,1} = InvRChar{11}; S.InvR {3,2} = InvRChar{12}; S.InvR {3,3} = InvRChar{13};
        end
        function Phi = mkMarkerCtrs(S)
            S.mkCharSegRef();
            NMarkers = size(S.LocalMarkers,1);
            Phi = {};
            for i = 1:NMarkers
                Phi = [Phi;S.mkMarkerCtr(S.LocalMarkers(i))];
            end
        end
        function Phi = mkMarkerCtr(S,LocalMarker)%,MarkerLocCoord,OrLocCoord)
            %MKMARKERCONSTRAIN Summary of this function goes here
            %   Detailed explanation goes here
            Phi ={};
            PIni    = S.LocalPoints(1).Point.CoordName;
            PFin    = LocalMarker.Point.CoordName;
            PLocIni = S.LocalPoints(1).CoordName;
            PLocFin = LocalMarker.CoordName;
%             PLocIni = S.LocalPoints(1).LocCoord;
%             PLocFin = LocalMarker.LocCoord;
            minus = '-';
            LocU = {[PLocFin{1},minus,PLocIni{1}];[PLocFin{2},minus,PLocIni{2}];[PLocFin{3},minus,PLocIni{3}]};
%             LocU = PLocFin - PLocIni;
            U    = {[PFin{1},minus,PIni{1}];[PFin{2},minus,PIni{2}];[PFin{3},minus,PIni{3}]};
            Phi  = [Phi; S.mkCtrCharLinComb(U, LocU)];
        end
        function calcMassMatrix(S,J)
            X = [S.LocalPoints(2).LocCoord - S.LocalPoints(1).LocCoord, S.LocalVectors(1).LocCoord,S.LocalVectors(2).LocCoord];
            a = inv(X)*(S.CoM.LocCoord-S.LocalPoints(1).LocCoord);
            Z = inv(X)*J*inv(X');
            M1  = (S.Mass-2*S.Mass*a(1)+Z(1,1))*eye(3);
            M2  = (S.Mass*a(1)-Z(1,1))*eye(3);
            M3  = (S.Mass*a(2)-Z(1,2))*eye(3);
            M4  = (S.Mass*a(3)-Z(1,3))*eye(3);
            M5  = M2;            
            M6  = Z(1,1)*eye(3);
            M7  = Z(1,2)*eye(3);
            M8  = Z(1,3)*eye(3);
            M9  = (S.Mass*a(2)-Z(2,1))*eye(3);
            M10 = Z(2,1)*eye(3);
            M11 = Z(2,2)*eye(3);
            M12 = Z(2,3)*eye(3);
            M13 = (S.Mass*a(3)-Z(3,1))*eye(3);
            M14 = Z(3,1)*eye(3);
            M15 = Z(3,2)*eye(3);
            M16 = Z(3,3)*eye(3);
            S.M =[ M1,  M2,  M3,  M4;
                   M5,  M6,  M7,  M8;
                   M9, M10, M11, M12;
                  M13, M14, M15, M16];                    
        end
        function setGraphicWireFrame(S,DrawSeq,Radius,Color)
            PosSep = findstr(',',DrawSeq);
            NSep = length(PosSep);
            S.WireFrame.DrawSeq = {DrawSeq(1:(PosSep(1)-1))};
            i = 0;
            if NSep>1
                for i=1:NSep-1
                    Name = DrawSeq((PosSep(i)+1):(PosSep(i+1)-1));
                    S.WireFrame.DrawSeq = [S.WireFrame.DrawSeq;Name];
                end
            end
            i=i+1;
            Name = DrawSeq((PosSep(i)+1):end);
            S.WireFrame.DrawSeq = [S.WireFrame.DrawSeq;Name];
            S.WireFrame.Radius = Radius;
            S.WireFrame.Color = Color;
        end
        function setGraphicFiles(S,GraphicName,GraphicTras,GraphicRot,GraphicScal)
            Graphic{1,1} = GraphicName;
%             Tras = calcCoord(GraphicTras);
%             Rot  = calcCoord(GraphicRot);
%             Scal = calcCoord(GraphicScal);
            Graphic{1,2} = GraphicTras;
            Graphic{1,3} = GraphicRot;
            Graphic{1,4} = GraphicScal;
            S.Graphic = [S.Graphic;Graphic];
        end
        function setMarkersAndALsInLCS(S,TCS_Pos_Or,LCS_R_TCS)
            NMarkers = size(S.LocalMarkers,1);
            NALandmarks = size (S.LocalALs,1);
            LCS_Pos_Markers = changeCoordSys(S.TCS_Pos_Markers,TCS_Pos_Or,LCS_R_TCS);
            LCS_Pos_AL      = changeCoordSys(S.TCS_Pos_Landmark,TCS_Pos_Or,LCS_R_TCS);
            for i=1:NMarkers
                S.LocalMarkers(i).LocCoord = LCS_Pos_Markers(:,i);
            end
            for i=1:NALandmarks
                S.LocalALs(i).LocCoord = LCS_Pos_AL(:,i);
            end
        end
        function setTypeOfBasis(S)
            % Calc number of points and vectors
            NVectors = size(S.LocalVectors,1);
            NPoints  = size(S.LocalPoints,1);
            if NPoints >= 2 && NVectors>=2
                S.BasisType = 1; % Basis formed by 2 vectors and 2 points
            elseif NVectors >= 3  && NPoints >=1
                S.BasisType = 0; % Basis formed by 3 vectors and 1 point
            else
                error(['Segment "',S.Name,'" must have at least 2 points & 2 vectors OR 1 point & 3 vectors']);
            end
            
        end
        function setGraphicPars(S)
            if ~isempty(S.Graphic)
                S.Graphic{1,2} = S.setGraphicCell(S.Graphic{1,2}); % Tras
                S.Graphic{1,3} = S.setGraphicCell(S.Graphic{1,3}); % Rot
                S.Graphic{1,4} = S.setGraphicCell(S.Graphic{1,4}); % Scal
                
                % store data in numeric format.
                S.Graphic{1,2} = str2num(S.Graphic{1,2});
                S.Graphic{1,3} = str2num(S.Graphic{1,3});
                S.Graphic{1,4} = str2num(S.Graphic{1,4});
            end
            
        end
        function GraphicCell = setGraphicCell(S,GraphicCell)
            if ~isnumeric(GraphicCell(2:end-1))
                str = GraphicCell;
                NLocalPoints = size(S.LocalPoints,1);
                for i=1:NLocalPoints
                    for j=1:3
                        StrCoord = num2str(S.LocalPoints(i).LocCoord(j));
                        str = strrep(str,S.LocalPoints(i).CoordName{j},StrCoord);
                    end
                end
                GraphicCell = str;
            end
        end
        function SegmentElement = writeSubjectParsIntermed(S,docNode)
            NAl= size(S.LocalALs,1);
            NPoints = size(S.LocalPoints,1);
            NMarkers = size(S.LocalMarkers,1);
            % Fill the segmetn attributes
            SegmentElement = docNode.createElement('Segment');
            MAttribute = docNode.createAttribute('Mass');
            NAttribute = docNode.createAttribute('Name');
            MAttribute.setValue(num2str(S.Mass));
            NAttribute.setValue(S.Name);
            SegmentElement.setAttributeNode(MAttribute);
            SegmentElement.setAttributeNode(NAttribute);
            % Fill all AL elements in segment and its attributes
            for i=1:NAl
                ALElement = SegmentElement.appendChild(docNode.createElement('AL'));
                ALNAttribute = docNode.createAttribute('Name');
                ALCoAttribute = docNode.createAttribute('XYZ_in_AL-LCS');
                ALNAttribute.setValue(S.LocalALs(i).Point.Name);
%                 B = [0;0.5;-12];
%                 SCoord= ['[',num2str(B(1)),';',num2str(B(2)),';',num2str(B(3)),']'];
                SCoord= ['[',num2str(S.LocalALs(i).LocCoord(1),'%-6.2f'),';',...
                    num2str(S.LocalALs(i).LocCoord(2),'%-6.2f'),';',...
                    num2str(S.LocalALs(i).LocCoord(3),'%-6.2f'),']'];
                ALCoAttribute.setValue(SCoord);
                ALElement.setAttributeNode(ALNAttribute);
                ALElement.setAttributeNode(ALCoAttribute);
            end
            % Fill all Points in segment and its attributes
            for i=1:NPoints
                PointElement = SegmentElement.appendChild(docNode.createElement('AL'));
                PointNAttribute = docNode.createAttribute('Name');
                PointCoAttribute = docNode.createAttribute('XYZ_in_AL-LCS');
                PointNAttribute.setValue(S.LocalPoints(i).Point.Name);
%                 B = [0;0.5;-12];
%                 SCoord= ['[',num2str(B(1)),';',num2str(B(2)),';',num2str(B(3)),']'];
                SCoord= ['[',num2str(S.LocalPoints(i).LocCoord(1),'%-6.2f'),';',...
                    num2str(S.LocalPoints(i).LocCoord(2),'%-6.2f'),';',...
                    num2str(S.LocalPoints(i).LocCoord(3),'%-6.2f'),']'];
                PointCoAttribute.setValue(SCoord);
                PointElement.setAttributeNode(PointNAttribute);
                PointElement.setAttributeNode(PointCoAttribute);
            end
            % Fill all Marker elements in segment and its attributes
            for i=1:NMarkers
                MarElement = SegmentElement.appendChild(docNode.createElement('Marker'));
                MarNAttribute = docNode.createAttribute('Name');
                MarCoAttribute = docNode.createAttribute('XYZ_in_AL-LCS');
                MarNAttribute.setValue(S.LocalMarkers(i).Point.Name);
%                 B = [0;0.5;-12];
%                 SCoord= ['[',num2str(B(1)),';',num2str(B(2)),';',num2str(B(3)),']'];
                SCoord= ['[',num2str(S.LocalMarkers(i).LocCoord(1),'%-6.2f'),';',...
                     num2str(S.LocalMarkers(i).LocCoord(2),'%-6.2f'),';',...
                     num2str(S.LocalMarkers(i).LocCoord(3),'%-6.2f'),']'];
                MarCoAttribute.setValue(SCoord);
                MarElement.setAttributeNode(MarNAttribute);
                MarElement.setAttributeNode(MarCoAttribute);
            end
            % Fill the CoM element
            CoMElement = SegmentElement.appendChild(docNode.createElement('CoM'));
            CoMCoAttribute = docNode.createAttribute('XYZ_in_AL-LCS');
%             B = [0;0.5;-12];
%             SCoord= ['[',num2str(B(1)),';',num2str(B(2)),';',num2str(B(3)),']'];
            SCoord= ['[',num2str(S.CoM.LocCoord(1),'%-6.2f'),';',...
                    num2str(S.CoM.LocCoord(2),'%-6.2f'),';',...
                    num2str(S.CoM.LocCoord(3),'%-6.2f'),']'];
            CoMCoAttribute.setValue(SCoord);
            CoMElement.setAttributeNode(CoMCoAttribute);
            % Fill the MoI element
%             S.I=ones(3,3);
            MoIElement = SegmentElement.appendChild(docNode.createElement('MoI'));
            MoIxxAttribute = docNode.createAttribute('Ixx'); MoIxxAttribute.setValue(num2str(S.I(1,1),'%-6.2f'));MoIElement.setAttributeNode(MoIxxAttribute);
            MoIxyAttribute = docNode.createAttribute('Ixy'); MoIxyAttribute.setValue(num2str(S.I(1,2),'%-6.2f'));MoIElement.setAttributeNode(MoIxyAttribute);
            MoIxzAttribute = docNode.createAttribute('Ixz'); MoIxzAttribute.setValue(num2str(S.I(1,3),'%-6.2f'));MoIElement.setAttributeNode(MoIxzAttribute);
            MoIyxAttribute = docNode.createAttribute('Iyx'); MoIyxAttribute.setValue(num2str(S.I(2,1),'%-6.2f'));MoIElement.setAttributeNode(MoIyxAttribute);
            MoIyyAttribute = docNode.createAttribute('Iyy'); MoIyyAttribute.setValue(num2str(S.I(2,2),'%-6.2f'));MoIElement.setAttributeNode(MoIyyAttribute);
            MoIyzAttribute = docNode.createAttribute('Iyz'); MoIyzAttribute.setValue(num2str(S.I(2,3),'%-6.2f'));MoIElement.setAttributeNode(MoIyzAttribute);
            MoIzxAttribute = docNode.createAttribute('Izx'); MoIzxAttribute.setValue(num2str(S.I(3,1),'%-6.2f'));MoIElement.setAttributeNode(MoIzxAttribute);
            MoIzyAttribute = docNode.createAttribute('Izy'); MoIzyAttribute.setValue(num2str(S.I(3,2),'%-6.2f'));MoIElement.setAttributeNode(MoIzyAttribute);
            MoIzzAttribute = docNode.createAttribute('Izz'); MoIzzAttribute.setValue(num2str(S.I(3,3),'%-6.2f'));MoIElement.setAttributeNode(MoIzzAttribute);                      
        end
    end
    
end

