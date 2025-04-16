classdef HUMAN_PARSER_MAT < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        StructMat         % Parameters of model in Local coord.  STRUCT
        mm2m = 0.001;     % Constant
        deg2rad = pi/180; % Constant
    end
    
    methods
        function HPM = HUMAN_PARSER_MAT()
            % Constructor of class
        end
        function parsePoints(HPM,Segments)
            NSegments = size(Segments,1);
            for i=1:NSegments
                if Segments(i).Fixed ~= 1 % If is not Ground
                %if ~strcmpi(Segments(i).Name,'Ground')
                    NPoints = size(Segments(i).LocalPoints,1);
                    Segments(i).LocalPoints(1).LocCoord = [0;0;0];
                    for j=2:NPoints
                        PName = Segments(i).LocalPoints(j).Point.Name;
                        Segments(i).LocalPoints(j).LocCoord = HPM.mm2m * (HPM.StructMat.(PName))';
                    end
                end
            end
        end
        function parseMarkers(HPM,Human)
            % Table of RAMSIS-RPx names
            BodyTranslator( 1,:) = {'pelvis','BEC'};
            BodyTranslator( 2,:) = {'lower lumbar spine','ULW'};
            BodyTranslator( 3,:) = {'upper lumbar spine','OLW'};
            BodyTranslator( 4,:) = {'lower thoracic spine','UBW'};
            BodyTranslator( 5,:) = {'upper thoracic spine','OBW'};
            BodyTranslator( 6,:) = {'lower cervical spine','UHW'};
            BodyTranslator( 7,:) = {'upper cervical spine','OHW'};
            BodyTranslator( 8,:) = {'head','KO'};
            BodyTranslator( 9,:) = {'chest','BRK'};
            BodyTranslator(10,:) = {'right clavicle','SBR'};
            BodyTranslator(11,:) = {'left clavicle','SBL'};
            BodyTranslator(12,:) = {'right upper arm','OAR'};
            BodyTranslator(13,:) = {'left upper arm','OAL'};
            BodyTranslator(14,:) = {'right lower arm','UAR'};
            BodyTranslator(15,:) = {'left lower arm','UAL'};
            BodyTranslator(16,:) = {'right hand','HAR'};
            BodyTranslator(17,:) = {'left hand','HAL'};
            BodyTranslator(18,:) = {'right thigh','OSR'};
            BodyTranslator(19,:) = {'left thigh','OSL'};
            BodyTranslator(20,:) = {'right lower leg','USR'};
            BodyTranslator(21,:) = {'left lower leg','USL'};
            BodyTranslator(22,:) = {'right foot','FUR'};
            BodyTranslator(23,:) = {'left foot','FUL'};
            BodyTranslator(24,:) = {'right ball of foot','FBR'};
            BodyTranslator(25,:) = {'left ball of foot','FBL'};
            % PUFO. Fingers are not currently in the model.
            % Markers on fingers are added to body Hand
            BodyTranslator(26,:) = {'right thumb finger 3','HAR'};
            BodyTranslator(27,:) = {'right index finger 3','HAR'};
            BodyTranslator(28,:) = {'right middle finger 3','HAR'};
            BodyTranslator(29,:) = {'left index finger 3','HAL'};
            
            NMarkers = length(HPM.StructMat.Marker);
            if (NMarkers == 0)
                error('There are not markers in Mat file');
            end
            for i=1:NMarkers
                MarkerName   = upper(HPM.StructMat.Marker(i).Name);
                SegmentIndex = getCellIndex(BodyTranslator,HPM.StructMat.Marker(i).BodyName);
                SegmentName  = BodyTranslator{SegmentIndex,2};
                SegmentIndex = getVecIndex(SegmentName,Human.Segments);
                if SegmentIndex>0
                    % add marker in list of markers in human 
                    AddedMarker = Human.addMarker(MarkerName);
                                        
                    %  Add Local Point to Segment and its local coordinates                                                 
                    Human.Segments(SegmentIndex).addMarker(AddedMarker);
                    NLocalMarker = size(Human.Segments(SegmentIndex).LocalMarkers,1);
                    Human.Segments(SegmentIndex).LocalMarkers(NLocalMarker).LocCoord = HPM.mm2m * (HPM.StructMat.Marker(i).LocalPosition)';
                else
%                     warning('Body %s is not in Human Model',SegmentName);
                end
            end
            % add marker in q
            Human.addMarkerInq();
        end
        function parseXMat(HPM,Human)
            
            Angle.GHZ_a1  = 0*HPM.deg2rad; Angle.GHZ_a2  = 0*HPM.deg2rad; Angle.GHZ_a3  = 0*HPM.deg2rad;
            Angle.GHUR_a3 = 0*HPM.deg2rad; Angle.GHUR_a2 = 0*HPM.deg2rad; Angle.GHUR_a1 = 0*HPM.deg2rad;
            Angle.GKNR_a2 = 0*HPM.deg2rad; Angle.GKNR_a1 = 0*HPM.deg2rad;
            Angle.GSPR_a3 = 0*HPM.deg2rad; Angle.GSPR_a1 = 0*HPM.deg2rad;
            Angle.GFBR_a3 = 0*HPM.deg2rad;
            Angle.GHUL_a3 = 0*HPM.deg2rad; Angle.GHUL_a2 = 0*HPM.deg2rad; Angle.GHUL_a1 = 0*HPM.deg2rad;
            Angle.GKNL_a2 = 0*HPM.deg2rad; Angle.GKNL_a1 = 0*HPM.deg2rad;
            Angle.GSPL_a3 = 0*HPM.deg2rad; Angle.GSPL_a1  = 0*HPM.deg2rad;
            Angle.GFBL_a3 = 0*HPM.deg2rad;
            Angle.GLK_a3 = 0*HPM.deg2rad; Angle.GLK_a2= 0*HPM.deg2rad;  Angle.GLK_a1 = 0*HPM.deg2rad;
            Angle.GLL_a3 = 0*HPM.deg2rad; Angle.GLL_a2 = 0*HPM.deg2rad; Angle.GLL_a1 = 0*HPM.deg2rad;
            Angle.GBL_a3 = 0*HPM.deg2rad; Angle.GBL_a2 = 0*HPM.deg2rad; Angle.GBL_a1 = 0*HPM.deg2rad;
            Angle.GBB_a3 = 0*HPM.deg2rad; Angle.GBB_a2 = 0*HPM.deg2rad; Angle.GBB_a1 = 0*HPM.deg2rad; 
            Angle.GHB_a3 = 0*HPM.deg2rad; Angle.GHB_a2 = 0*HPM.deg2rad; Angle.GHB_a1= 0*HPM.deg2rad;
            Angle.GHH_a3 = 0*HPM.deg2rad; Angle.GHH_a2 = 0*HPM.deg2rad; Angle.GHH_a1= 0*HPM.deg2rad;
            Angle.GKH_a3 = 0*HPM.deg2rad; Angle.GKH_a2 = 0*HPM.deg2rad; Angle.GKH_a1= 0*HPM.deg2rad;
            Angle.GBRK_a3 = 0*HPM.deg2rad;
            Angle.GSBR_a3 = 0*HPM.deg2rad; Angle.GSBR_a2 = 0*HPM.deg2rad; Angle.GSBR_a1 = 0*HPM.deg2rad;
            Angle.GSBL_a3 = 0*HPM.deg2rad; Angle.GSBL_a2 = 0*HPM.deg2rad; Angle.GSBL_a1 = 0*HPM.deg2rad;
            Angle.GSR_a3  = 0*HPM.deg2rad;  Angle.GSR_a2 = 0*HPM.deg2rad; Angle.GSR_a1  = 0*HPM.deg2rad;
            Angle.GSL_a3  = 0*HPM.deg2rad;  Angle.GSL_a2 = 0*HPM.deg2rad; Angle.GSL_a1  = 0*HPM.deg2rad;
            Angle.GELR_a2 = 0*HPM.deg2rad; Angle.GELR_a1 = 0*HPM.deg2rad;
            Angle.GELL_a2 = 0*HPM.deg2rad; Angle.GELL_a1 = 0*HPM.deg2rad;
            Angle.GHAR_a3 = 0*HPM.deg2rad; Angle.GHAR_a2 = 0*HPM.deg2rad;
            Angle.GHAL_a3 = 0*HPM.deg2rad; Angle.GHAL_a2 = 0*HPM.deg2rad;
            
            %             Human.Angles(1).a1 = GHZ_a1;    Human.Angles(1).a2 = GHZ_a2;    Human.Angles(1).a3 = GHZ_a3;
            %             Human.Angles(2).a1 = GHUR_a1;    Human.Angles(2).a2 = GHUR_a2;    Human.Angles(2).a3 = GHUR_a3;
            %             Human.Angles(3).a1 = GKNR_a1;    Human.Angles(3).a2 = GKNR_a2;
            %             Human.Angles(4).a1 = GSPR_a1;    Human.Angles(4).a3 = GSPR_a3;
            %             Human.Angles(5).a1 = GHUL_a1;    Human.Angles(5).a2 = GHUL_a2;    Human.Angles(5).a3 = GHUL_a3;
            %             Human.Angles(6).a1 = GKNL_a1;    Human.Angles(6).a2 = GKNL_a2;
            %             Human.Angles(7).a1 = GSPL_a1;    Human.Angles(7).a3 = GSPL_a3;
            
            NAngles = size(Human.Angles,1);
            for i =1:NAngles
                if(strcmpi(Human.Angles(i).Joint.Type,'SPH'))
                   Human.Angles(i).a1 = Angle.(Human.Angles(i).Name1);
                   Human.Angles(i).a2 = Angle.(Human.Angles(i).Name2);
                   Human.Angles(i).a3 = Angle.(Human.Angles(i).Name3);
                elseif (strcmpi(Human.Angles(i).Joint.Type,'UNI'))
                   TypeAngle1  = Human.Angles(i).Name1(end-1:end);
                   TypeAngle2  = Human.Angles(i).Name2(end-1:end);
                   if (strcmpi(TypeAngle1,'a1'))
                       Human.Angles(i).a1  = Angle.(Human.Angles(i).Name1);
                   elseif (strcmpi(TypeAngle1,'a2'))
                       Human.Angles(i).a2  = Angle.(Human.Angles(i).Name1);
                   elseif (strcmpi(TypeAngle1,'a3')) 
                       Human.Angles(i).a3  = Angle.(Human.Angles(i).Name1);
                   end
                   if (strcmpi(TypeAngle2,'a1'))
                       Human.Angles(i).a1  = Angle.(Human.Angles(i).Name2);
                   elseif (strcmpi(TypeAngle2,'a2'))
                       Human.Angles(i).a2  = Angle.(Human.Angles(i).Name2);
                   elseif (strcmpi(TypeAngle2,'a3')) 
                       Human.Angles(i).a3  = Angle.(Human.Angles(i).Name2);
                   end
                elseif (strcmpi(Human.Angles(i).Joint.Type,'REV'))
                   TypeAngle  = Human.Angles(i).Name1(end-1:end);
                   if (strcmpi(TypeAngle,'a1'))
                       Human.Angles(i).a1  = Angle.(Human.Angles(i).Name1);
                   elseif (strcmpi(TypeAngle,'a2'))
                       Human.Angles(i).a2  = Angle.(Human.Angles(i).Name1);
                   elseif (strcmpi(TypeAngle,'a3')) 
                       Human.Angles(i).a3  = Angle.(Human.Angles(i).Name1);
                   end
                end
            end

            
            Glob_R_BEC = rot123s(Angle.GHZ_a1,Angle.GHZ_a2,Angle.GHZ_a3);
            BEC_R_OSR = rot321s(Angle.GHUR_a3,Angle.GHUR_a2,Angle.GHUR_a1);
            OSR_R_USR = rot21(Angle.GKNR_a2,Angle.GKNR_a1);
            USR_R_FUR = rot31(Angle.GSPR_a3,Angle.GSPR_a1);
            FUR_R_FBR = rot_z(Angle.GFBR_a3);
            BEC_R_OSL = rot321s(Angle.GHUL_a3,Angle.GHUL_a2,Angle.GHUL_a1);
            OSL_R_USL = rot21(Angle.GKNL_a2,Angle.GKNL_a1);
            USL_R_FUL = rot31(Angle.GSPL_a3,Angle.GSPL_a1);
            FUL_R_FBL = rot_z(Angle.GFBL_a3);
%             BEC_R_ULW = rot32(Angle.GLK_a3,Angle.GLK_a2);
            BEC_R_ULW = rot321s(Angle.GLK_a3,Angle.GLK_a2,Angle.GLK_a1);
            OLW_R_ULW = rot321s(Angle.GLL_a3,Angle.GLL_a2,Angle.GLL_a1);
            ULW_R_UBW = rot321s(Angle.GBL_a3,Angle.GBL_a2,Angle.GBL_a1);
            UBW_R_OBW = rot321s(Angle.GBB_a3,Angle.GBB_a2,Angle.GBB_a1);
            OBW_R_UHW = rot321s(Angle.GHB_a3,Angle.GHB_a2,Angle.GHB_a1);
            UHW_R_OHW = rot321s(Angle.GHH_a3,Angle.GHH_a2,Angle.GHH_a1);
            OHW_R_KO  = rot321s(Angle.GKH_a3,Angle.GKH_a2,Angle.GKH_a1);
            UHW_R_BRK = rot_z(Angle.GBRK_a3);
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
            Glob_R_FBR = Glob_R_FUR * FUR_R_FBR;                                RigidTrans = [RigidTrans; {'FBR',Glob_R_FBR}];
            Glob_R_FBL = Glob_R_FUL * FUL_R_FBL;                                RigidTrans = [RigidTrans; {'FBL',Glob_R_FBL}];
            Glob_R_ULW = Glob_R_BEC * BEC_R_ULW;                                RigidTrans = [RigidTrans; {'ULW',Glob_R_ULW}];
            Glob_R_OLW = Glob_R_ULW * OLW_R_ULW;                                RigidTrans = [RigidTrans; {'OLW',Glob_R_OLW}];
            Glob_R_UBW = Glob_R_OLW * ULW_R_UBW;                                RigidTrans = [RigidTrans; {'UBW',Glob_R_UBW}];
            Glob_R_OBW = Glob_R_UBW * UBW_R_OBW;                                RigidTrans = [RigidTrans; {'OBW',Glob_R_OBW}];
            Glob_R_UHW = Glob_R_OBW * OBW_R_UHW;                                RigidTrans = [RigidTrans; {'UHW',Glob_R_UHW}];
            Glob_R_OHW = Glob_R_UHW * UHW_R_OHW;                                RigidTrans = [RigidTrans; {'OHW',Glob_R_OHW}];
            Glob_R_KO  = Glob_R_OHW * OHW_R_KO;                                 RigidTrans = [RigidTrans; {'KO',Glob_R_KO}];
            Glob_R_BRK = Glob_R_UHW *  rot_z(pi/2) * rot_x(pi) * UHW_R_BRK;     RigidTrans = [RigidTrans; {'BRK',Glob_R_BRK}];
            Glob_R_SBR = Glob_R_BRK * rot_y(-pi/2) * rot_x(-pi/2) * BRK_R_SBR;  RigidTrans = [RigidTrans; {'SBR',Glob_R_SBR}];
            Glob_R_SBL = Glob_R_BRK * rot_y(pi/2) * rot_x(pi/2) *   BRK_R_SBL;  RigidTrans = [RigidTrans; {'SBL',Glob_R_SBL}];
            Glob_R_OAR = Glob_R_SBR * SBR_R_OAR;                                RigidTrans = [RigidTrans; {'OAR',Glob_R_OAR}];
            Glob_R_OAL = Glob_R_SBL * SBL_R_OAL;                                RigidTrans = [RigidTrans; {'OAL',Glob_R_OAL}];
            Glob_R_UAR = Glob_R_OAR * rot_x(-pi/2) * OAR_R_UAR * rot_x(pi/2);   RigidTrans = [RigidTrans; {'UAR',Glob_R_UAR}];
            Glob_R_UAL = Glob_R_OAL * rot_x(-pi/2) * OAL_R_UAL * rot_x(pi/2);   RigidTrans = [RigidTrans; {'UAL',Glob_R_UAL}];
            Glob_R_HAR = Glob_R_UAR * UAR_R_HAR;                                RigidTrans = [RigidTrans; {'HAR',Glob_R_HAR}];
            Glob_R_HAL = Glob_R_UAL * UAL_R_HAL;                                RigidTrans = [RigidTrans; {'HAL',Glob_R_HAL}];
            
            NTrans = size(RigidTrans);
            for i=1:NTrans
                SegmentIndex = getVecIndex(RigidTrans{i,1},Human.Segments);
                if SegmentIndex ~= 0
                    Human.Segments(SegmentIndex).Glob_R_Seg =  RigidTrans{i,2};
                end
            end
        end
        function readMat(HPM,CalibPath,CalibFile)
            Tmp = load([CalibPath, CalibFile]);
            HPM.StructMat = Tmp.dim;
        end

        
    end
    
end

