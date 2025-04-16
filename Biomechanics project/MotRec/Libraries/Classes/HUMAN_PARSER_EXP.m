classdef HUMAN_PARSER_EXP < handle
    %HUMAN_PARSER_EXP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        StructExp         % Parameters of model in global coord     STRUCT
        mm2m = 0.001;     % Constant
        deg2rad = pi/180; % Constant
    end
    
    methods
        function HPexp = HUMAN_PARSER_EXP()
        end
        function checkPointsData(HPexp,Points,File)
            % CHECKPOINTSDATA check if all points of the model have value in exp file
            NModelPoints = size(Points,1);
            for i=1:NModelPoints
                if (Points(i).Fixed == 0)
                    PointName = Points(i).Name;
                    if ~isfield (HPexp.StructExp.Sklt,PointName)
                        error(['Point "',Points(i).Name,'" is missing in file: ',File])
                    end
                end
            end
        end
        function parseAngles(HPexp,Angles)
            % go over about angle list in human and asign its value
            NModelAngles = size(Angles,1);
            for i=1:NModelAngles
                AngleName = Angles(i).Name1;
                PosOf_    = findstr('_',AngleName);
                JointName = AngleName(1:(PosOf_-1));
                Angles(i).a1 = HPexp.StructExp.Jnt.(JointName).rx * HPexp.deg2rad;
                Angles(i).a2 = HPexp.StructExp.Jnt.(JointName).ry * HPexp.deg2rad;
                Angles(i).a3 = HPexp.StructExp.Jnt.(JointName).rz * HPexp.deg2rad;
            end
        end
        function parsePoints(HPexp,Points)
            % We find point from poin list in human in EXP and asign its global coordinate value to point 
            NModelPoints = size(Points,1);
            for i=1:NModelPoints
                if (Points(i).Fixed == 0)
                    PointName = Points(i).Name;
                    Points(i).GlobalCoord = HPexp.mm2m * [HPexp.StructExp.Sklt.(PointName).x; ...
                                                          HPexp.StructExp.Sklt.(PointName).y; ...
                                                          HPexp.StructExp.Sklt.(PointName).z];
                end
            end
        end
        function parseVectors(HPexp,Vectors)
            % go through vector list in HUMAN and asign its global coordinate value from file EXP.
            % We also need global position of x and y axes global coordinates.
            % Hence we add vector to vector list in HUMAN from file EXP.
            NModelVectors = size(Vectors,1);
            for i=1:NModelVectors
                if (Vectors(i).Fixed == 0)
                    VectorName = Vectors(i).Name;
                    [JointName, TypeOfVector] = strtok(VectorName,'_');
                    TypeOfVector = TypeOfVector(2:end);
                    
                    if strcmpi(TypeOfVector,'v1') || strcmpi(TypeOfVector,'v2') || strcmpi(TypeOfVector,'v3') || ...
                       strcmpi(TypeOfVector,'x')  || strcmpi(TypeOfVector,'y')  || strcmpi(TypeOfVector,'z')
                        
                        Vx=HPexp.StructExp.Jnt.(JointName).([TypeOfVector,'x']);
                        Vy=HPexp.StructExp.Jnt.(JointName).([TypeOfVector,'y']);
                        Vz=HPexp.StructExp.Jnt.(JointName).([TypeOfVector,'z']);
                        Norm = norm([Vx;Vy;Vz]);
                        Vx=Vx/Norm;
                        Vy=Vy/Norm;
                        Vz=Vz/Norm;
                        Vectors(i).GlobalCoord = [Vx;Vy;Vz];
                        
                    else
                        warning('Vector "%s" is not in EXP file, it will be initialised to [0;0;0]\n',VectorName);
                        Vectors(i).GlobalCoord = [0;0;0];
                    end
                end
            end
           
        
        end
        function parseMarkers(HPexp,HM)
            NExpMarkers = size(HPexp.StructExp.listMrk,2);
            if (NExpMarkers == 0)
                error('There are not markers in EXP file');
            end
            for i=1:NExpMarkers
                % find if segmet is in modelposition 
                MarkerName = HPexp.StructExp.listMrk{i};
                SegmentName = HPexp.StructExp.Mrk.(MarkerName).body;
                SegmentIndex = getVecIndex(SegmentName,HM.Segments);
                if SegmentIndex>0
                    % add marker in list of markers in human and asign his global coordinate value
                    AddedMarker = HM.addMarker(MarkerName);
                    NModelMarkers = size(HM.Markers,1);
                    HM.Markers(NModelMarkers).GlobalCoord = HPexp.mm2m * [HPexp.StructExp.Mrk.(MarkerName).x;...
                                                                          HPexp.StructExp.Mrk.(MarkerName).y;...
                                                                          HPexp.StructExp.Mrk.(MarkerName).z];
                    %  Add Local Point to Segment                                                  
                    HM.Segments(SegmentIndex).addMarker(AddedMarker);
                else
%                     warning('Body %s is not in Human Model',SegmentName);
                end
            end
            % add marker in q
            HM.addMarkerInq();
        end
        function parseXMat(HPexp,HM)
            NSegments = size(HM.Segments,1);
            for i=1:NSegments
                % Not include Ground, there are not necesary
                if HM.Segments(i).Fixed ~= 1 % If is not Ground
                %if ~strcmpi(HM.Segments(i).Name,'Ground')
                    % All the time the first LocalPoints is the origin
                    OrName = HM.Segments(i).LocalPoints(1).Point.Name;
                    Glob_Pos_OrSegment = HM.Segments(i).LocalPoints(1).Point.GlobalCoord;
                    % get global coordinates of vectors
                    Glob_Vec_x = [HPexp.StructExp.Jnt.(OrName).xx;HPexp.StructExp.Jnt.(OrName).xy;HPexp.StructExp.Jnt.(OrName).xz];
                    Glob_Vec_y = [HPexp.StructExp.Jnt.(OrName).yx;HPexp.StructExp.Jnt.(OrName).yy;HPexp.StructExp.Jnt.(OrName).yz];
                    HM.Segments(i).calcRd(Glob_Pos_OrSegment, Glob_Vec_x, Glob_Vec_y);
                end
            end
        end
        function parseSegmentsParameters(HPexp,HM)
            % We find Segment from Segments list in Human and assign mass and Global coordinates of CoM of Segment
            NModelSeg = size(HM.Segments,1);
            for i=1:NModelSeg
                if HM.Segments(i).Fixed ~= 1 % If is not Ground
                %if ~strcmpi(HM.Segments(i).Name,'Ground')
                    SegmentName = HM.Segments(i).Name;
                    EXPSegmentName = [SegmentName,'cog']; % all body names end with "cog" e.g. BECcog
                    HM.Segments(i).Mass = HPexp.StructExp.Body.(EXPSegmentName).mass;
                    MoIxx = HPexp.StructExp.Body.(EXPSegmentName).MoIxx;
                    MoIxy = HPexp.StructExp.Body.(EXPSegmentName).MoIxy;
                    MoIxz = HPexp.StructExp.Body.(EXPSegmentName).MoIxz;
                    MoIyx = HPexp.StructExp.Body.(EXPSegmentName).MoIyx;
                    MoIyy = HPexp.StructExp.Body.(EXPSegmentName).MoIyy;
                    MoIyz = HPexp.StructExp.Body.(EXPSegmentName).MoIyz;
                    MoIzx = HPexp.StructExp.Body.(EXPSegmentName).MoIzx;
                    MoIzy = HPexp.StructExp.Body.(EXPSegmentName).MoIzy;
                    MoIzz = HPexp.StructExp.Body.(EXPSegmentName).MoIzz;
                    HM.Segments(i).I = [MoIxx MoIxy MoIxz; MoIyx MoIyy MoIyz; MoIzx MoIzy MoIzz]*0.0001;% Pass from kg.cm2 to kg.m2
                    % Add point in list of points in Human
                    AddedPoint = HM.addCoM(SegmentName);
                    NModelCoM = size(HM.CoMs,1);
                    HM.CoMs(NModelCoM).GlobalCoord = HPexp.mm2m * [HPexp.StructExp.Body.(EXPSegmentName).x;...
                                                                  HPexp.StructExp.Body.(EXPSegmentName).y;...
                                                                  HPexp.StructExp.Body.(EXPSegmentName).z];
                    % Add point in Segment
                    HM.Segments(i).CoM = LOCAL_POINT(AddedPoint,HM.Segments(i).Name);
                end
            end
        end
        function readExp(HPexp,Path,Filename)
            fichExp = [Path,Filename];
            
            exp=struct([]);
            if(nargin==0),  help loadExp;       return, end,
            %lecture fichier
            fid = fopen(fichExp, 'r'); str = char(fread(fid)'); fclose(fid);
            if ~isequal(str(1:6),'RAMSIS'), exp = loadXml( fichExp ); return, end,
            
            try,    tgtSubj =   textread(fichExp,'%s','delimiter','\n','whitespace','');
            catch,  disp(['pb de lecture de : ',fichExp]);   return,
            end,
            
            nLignes=length(tgtSubj);
            
            sec.SktPnt.Deb=0;	sec.SknPnt.Deb=0;	sec.Jnt.Deb=0;	sec.Mrk.Deb=0;
            sec.Cogs.Deb=0;	    sec.Mass.Deb=0;     sec.JCS.Deb=0;
            sec.SktPnt.Fin=0;	sec.SknPnt.Fin=0;	sec.Jnt.Fin=0;	sec.Mrk.Fin=0;
            sec.Cogs.Fin=0;	    sec.Mass.Fin=0;     sec.JCS.Fin=0;
            sec.MoI.Deb=0;      sec.MoI.Fin=0;
            
            for iLigne=1:nLignes,
                if( isequal( sscanf(tgtSubj{iLigne},'%s',1) , 'SECTION' ) ),
                    switch(sscanf(tgtSubj{iLigne},'%*s %s',1)),
                        case 'SKELETON_POINTS',
                            sec.SktPnt.Deb=iLigne;      sec.SktPnt.Fin= HPexp.scanEndSections(iLigne,tgtSubj);
                            iLigne=sec.SktPnt.Fin;
                        case 'JOINTS',
                            sec.Jnt.Deb=iLigne;      sec.Jnt.Fin= HPexp.scanEndSections(iLigne,tgtSubj);
                            iLigne=sec.Jnt.Fin;
                        case 'JOINTS_COORDINATE_SYSTEMS',
                            sec.JCS.Deb=iLigne;      sec.JCS.Fin= HPexp.scanEndSections(iLigne,tgtSubj);
                            iLigne=sec.JCS.Fin;
                        case 'BODY_ELEMENT_COORDINATE_SYSTEMS',
                            sec.BCS.Deb=iLigne;      sec.BCS.Fin= HPexp.scanEndSections(iLigne,tgtSubj);
                            iLigne=sec.BCS.Fin;
                        case 'SKIN_POINTS',
                            sec.SknPnt.Deb=iLigne;      sec.SknPnt.Fin= HPexp.scanEndSections(iLigne,tgtSubj);
                            iLigne=sec.SknPnt.Fin;
                        case 'MARKER',
                            sec.Mrk.Deb=iLigne;      sec.Mrk.Fin= HPexp.scanEndSections(iLigne,tgtSubj);
                            iLigne=sec.Mrk.Fin;
                        case 'CENTRES_OF_GRAVITY',
                            sec.Cogs.Deb=iLigne;      sec.Cogs.Fin= HPexp.scanEndSections(iLigne,tgtSubj);
                            iLigne=sec.Cogs.Fin;
                        case 'MASS',
                            sec.Mass.Deb=iLigne;      sec.Mass.Fin= HPexp.scanEndSections(iLigne,tgtSubj);
                            iLigne=sec.Mass.Fin;
                        case 'MOMENT_OF_INERTIA'
                            sec.MoI.Deb = iLigne;     sec.MoI.Fin= HPexp.scanEndSections(iLigne,tgtSubj);
                            iLigne=sec.MoI.Fin;
                    end,
                end,
                if(iLigne==-1), error(['loadExp :: wrong SECTION format for : ', ...
                        fichExp]);   return, end,
            end,
            clear exp;
            exp.LoadFileName=fichExp;
            exp.versionInfo=tgtSubj{1};
            
            valGrp='Sklt';
            valFmts={' %s', ' %f %f %f'};
            valNams={{'fullName'}, {'x', 'y', 'z'}};
            exp = HPexp.fillStructFmt(exp, tgtSubj, sec.SktPnt, valGrp, valFmts, valNams);
            %
            valGrp='Jnt';
            valFmts={' %f %f %f', ' %f %f %f'};
            valNams={{'x', 'y', 'z'},   {'rx', 'ry', 'rz'}};
            exp = HPexp.fillStructFmt(exp, tgtSubj, sec.Jnt, valGrp, valFmts, valNams);
            valFmts={' %f %f %f', ' %f %f %f', ' %f %f %f', ' %f %f %f'};
            valNams={{}, {'v1x', 'v1y', 'v1z'},{'v2x', 'v2y', 'v2z'},{'v3x', 'v3y', 'v3z'}};
            exp = HPexp.fillStructFmt(exp, tgtSubj, sec.JCS, valGrp, valFmts, valNams);
            valFmts={' %f %f %f', ' %f %f %f', ' %f %f %f', ' %f %f %f'};
            valNams={{}, {'xx', 'xy', 'xz'},{'yx', 'yy', 'yz'},{'zx', 'zy', 'zz'}};
            exp = HPexp.fillStructFmt(exp, tgtSubj, sec.BCS, valGrp, valFmts, valNams);
            %  Lecture de la section skin implementee mais desactive pour raison de
            %  performance de loadExp (section très longue à lire)
            valGrp='Skin';
            valFmts={' %f %f %f'};
            valNams={{'x', 'y', 'z'}};
            exp = HPexp.fillStructFmt(exp, tgtSubj, sec.SknPnt, valGrp, valFmts, valNams);
            %
            valGrp='Mrk';
%             valFmts={' %s', ' %s', ' %f %f %f'};
%             valNams={{'body'}, {'surfSec'}, {'x', 'y', 'z'}};
            valFmts={' %s',  ' %f %f %f'};
            valNams={{'body'},  {'x', 'y', 'z'}};
            exp = HPexp.fillStructFmt(exp, tgtSubj, sec.Mrk, valGrp, valFmts, valNams);
            %
            valGrp='Body';
            valFmts={' %f %f %f'};
            valNams={{'x' 'y' 'z'}};
            exp = HPexp.fillStructFmt(exp, tgtSubj, sec.Cogs, valGrp, valFmts, valNams);
            valFmts={' %f'};
            valNams={{'mass'}};
            exp = HPexp.fillStructFmt(exp, tgtSubj, sec.Mass, valGrp, valFmts, valNams);
            %
            valFmts={' %f %f %f %f %f %f %f %f %f'};
            valNams={{'MoIxx', 'MoIxy', 'MoIxz', 'MoIyx', 'MoIyy', 'MoIyz', 'MoIzx', 'MoIzy', 'MoIzz'}};
            exp = HPexp.fillStructFmt(exp, tgtSubj, sec.MoI, valGrp, valFmts, valNams);
            HPexp.StructExp = exp;
            
            return,
        end
        function [iLigneFin]= scanEndSections(HPexp,iLigneDeb,fileText) %#ok<MANU>
            iLigne=iLigneDeb;
            %
            flagFound=0;
            while(~flagFound & iLigne<length(fileText)),
                iLigne=iLigne+1;
                if( isequal( sscanf(fileText{iLigne},'%s',1) , 'ENDSECTION' ) ),
                    flagFound=1;
                end,
            end,
            %
            if (~flagFound),
                iLigneFin=-1;
            else,   iLigneFin=iLigne;   end,
            
            return,
        end
        function [skipFmts] = makeSkipFmt(HPexp,prevFmts) %#ok<MANU>
            indParts=strfind(prevFmts,'%');
            skipFmts='';    lastInd=0;
            for iPart=1:length(indParts),
                skipFmts=[skipFmts prevFmts(lastInd+1:indParts(iPart)) '*'];
                lastInd=indParts(iPart);
            end,
            skipFmts=[skipFmts prevFmts(lastInd+1:end)];
            
            return,
        end            
        function [exp]= fillStructFmt(HPexp,exp, tgtSubj, curSec, valGrp, valFmts, valNams)
            % Lecture des blocks 'SECTION' et enregistrement dans la structure exp
            valLst=['list' valGrp]; exp.(valLst)={};
            %
            readyStruct=struct('temp','tt');    readyStruct=rmfield(readyStruct,'temp');
            %
            if(~isfield(exp,valGrp)),    exp.(valGrp)=readyStruct;  end,
            for iLigne=curSec.Deb+1:curSec.Fin-1,
                curLigne=tgtSubj{iLigne};
                valName     =sscanf(curLigne,'%s',1);
                if isequal(valGrp,'Body') | isequal(valGrp,'mass'), valName=[valName 'cog']; end,
                evL=exp.(valLst);
                evL{end+1}=valName;
                exp.(valLst) =evL;
                %     exp.(valLst) ={exp.(valLst){:} valName};
                evG=exp.(valGrp);
                %     try,
                %         vals=exp.(valGrp).(valName);
                %         if(isempty(vals)),  vals=readyStruct;   end,
                %     catch,
                %         vals=readyStruct;
                %     end,
                
                if(isfield(evG,valName)),
                    vals=evG.(valName);
                    if(isempty(vals)),     vals=readyStruct;   end,
                else,   vals=readyStruct;   end,
                
                for iFmt=1:length(valFmts),
                    skipFmts    = HPexp.makeSkipFmt(['%s' valFmts{1:iFmt-1}]);
                    fullFmts    = [skipFmts valFmts{iFmt}];
                    nVals   = length(valNams{iFmt});
                    readVals    =sscanf(curLigne, fullFmts, nVals);
                    switch(nVals),
                        case 1, % une seule valeur numérique ou une chaine de caractères
                            vals.(valNams{iFmt}{1})=readVals;
                        otherwise   % plusieurs valeurs numériques
                            for iVal=1:nVals,
                                vals.(valNams{iFmt}{iVal})=readVals(iVal);
                            end,
                    end,
                end,
                exp.(valGrp).(valName)=vals;
            end,
            
            return,
        end

                        
    end
    
end
