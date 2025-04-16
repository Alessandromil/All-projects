classdef RECONSTRUCTION <handle
    %RECONSTRUCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Subject         % Model with parameters                                         SUBJECT
        Motion          % Motion dats                                                   Motion.Path/ Motion.File (CSV)
        SolverType      % Type of Solver                                                char
        Markers         % Markers for recontruc.                                        POINT[]
        DataCSV         % Contain all the data of the CSV
                        % DataCSV.MarkerTraj contain MarkerNames,MarkerCoords,Frecuency
                        % DataCSV.Analog contain Frecuency, and structure with Forces, Torques and EMG
                        %                Each structure have units and values. Example: DataCSV.Analog.Fx1
                        % DataCSV.PMap contain Units, Seat, and Back, Seat and Back contain values
        PoszInq         % Position of guided variables z in q                           double[]
        ExperimentName  % The name of experiment                                        char
        ExperimentPath  % The path of experiment                                        char
        ModelEqPath     % Model Equations path                                          char
        AddCtrPath      % Path of the additional constrains                             char
        PressurePath    % Path of the files of seat pressure maps                       char
        Deltat          % the sample time of the motion capture                         double [seconds]             
        g_t             % Vector with driven variables values in t                      double[nFrames x nGuidedVariables]
        gs_t            % Vector with value of driven coord in in equality constraints
        FilePhiName     % Name of the file Phi.m name                                   char
        FilePhiqName    % Name of the file Phiq.m name                                  char
        Wm              % Weighting factor of guided variables                          cell{}
        Ws              % Weighting factor to each driven coord include in equality constraints
        z               % Cell with driven variables                                    char{}
        Par             % Vector with point and marker coord value of the subject       double[]
        q0              % Value of  q in initial instant                                double[]
        q0User          % Value of q in initial instant provided by user                double[]
        q_t             % Value of eachgeneralized coordinate at each sample time.      double [NFrames x NVars]
        qdot_t          % The first time derivative of the generalized coordinates.     double [NFrames x NVars]
        qdot2_t         % The 2nd time derivative of the generalized coordinates.       double [NFrames x NVars]
        ExperLogFileId  % The files .log are text files that containt information about the experiment loop.
                        % It is a history of the reconstruction and contains also warning and errors ocurred during.
        FileLogName     % The Name of log file        
        ErrorID = 0;    % If ErrorID = 1, imposible to do ID.
        ResultType      % Stores whick calculations have been done. Options: 'IK' or 'ID'
        RawMarkerCoords % The raw coordinates of all the markers in the motion file
        RawMarkerNames  % The names of all the markers in motion file
        MarkerSensor    %###
        Settings        % Struct that contains settings for motion reconstruction
                        % Settings.Type. Inverse Kinematics or Inverse Dynamics. Possible options: IK, ID for IK+ID
                        % Settings.Display. Define ammount of feedback to user: 0(minimum) 1(standard) 2(complete information)
                        %   ALL RESULTS VARIABLES CAN HAVE 2 POSSIBLE INTEGER VALUES: 1(Yes)/0(No).
                        % Settings.Results.Ramsis. Turns on/off Ramsis generation results. 
                        % Settings.Results.PAM. Turns on/off Ramsis generation results.
                        % Settings.Results.InitPosture. Posture used for initialization in Compamm format (*_t0.pb & *_t0.sim)
                        % Settings.Results.CompPlayback. Reconstructed motion in Compamm format (*.pb & *.sim)
                        % Settings.Results.Position. Position of all model elements (markers, vectors, points) in Compamm format (*.pos)
                        % Settings.Results.Velocity. Velocity of all model elements (markers, vectors, points) in Compamm format (*.vel)
                        % Settings.Results.Acceleration. Acceleration of all model elements (markers, vectors, points) in Compamm format (*.acc)
                        % Settings.Results.MarkerError. Distance between experimental-markers and model-markers in Compamm format (*.dis)
                        % Settings.Results.RawMarkerTrajectory. Experimental marker trajectories in Compamm format (extracted from oringal data) (*.traj)
                        % Settings.Results.AcondMarkerTrajectory = 0; % Aconditioned (smoothed & gaps filled) marker trajectories in Compamm format (*.amt)
                        % Settings.Results.Sensor. Variables measured by sensors in Compamm format (*.sen)
                        % Settings.Results.JointEffort. Forces & torques at all the joints in the model. Includes reactions & motor efforts
                        % Settings.Interpolation.Method. Method for interpolating gaps in marker trajectories. Options:`'linear'
                        % Settings.Interpolation.NInterFrames. Threshold (in frames) for interpolating missing markers. Options: any positive integer
                        % Settings.Smoothing.Method. Smoothing method for marker trajectories. Options: 'butter' (butterworth filter)
                        % Settings.Smoothing.CutFreq. Cutt-off frequency in Hz for butterworth filter
        Model           % Model information                 Model.Path
                        %                                   Model.File / Name of model.m

    end
    
    methods
        function R = RECONSTRUCTION(Subject,Motion,PressurePath,SolverType,ModelEqPath,AddCtrPath,ExperimentPath,ExperimentName,ExperVar,ExperLogFileId,Settings,Model)
            R.Subject = Subject;
            R.Motion  = Motion;
            R.PressurePath = PressurePath;
            R.SolverType = SolverType;
            R.ModelEqPath = ModelEqPath;
            R.AddCtrPath = AddCtrPath;
            R.ExperimentPath = ExperimentPath;
            R.ExperimentName = ExperimentName;
            R.z =  ExperVar.z;
            R.Markers = ExperVar.Markers;
            R.Wm = ExperVar.Wm;
            R.Ws = ExperVar.Ws;
            R.PoszInq = ExperVar.PoszInq;
            R.gs_t = ExperVar.gs_t;
            R.FilePhiName = ExperVar.FilePhiName;
            R.FilePhiqName = ExperVar.FilePhiqName;
            R.Par = ExperVar.Par;
            R.ExperLogFileId = ExperLogFileId;
            R.FileLogName = fopen(ExperLogFileId);
            R.Settings = Settings;
            R.Settings.ExType = [];
            R.Model = Model;
                       
        end
        function addRecontructionMarkers(R,GuidedMarkers)
            % Add Markers for reconstruction if they are defined
            if ~isempty(GuidedMarkers)
                if strcmpi(GuidedMarkers{1},'AllMarkers')
                    R.Markers = R.Subject.Markers;
                    NMarkers = size(R.Markers,1);
                    for i=1:NMarkers
                        R.Markers(i).Wm = GuidedMarkers{1,2};
                    end
                else
                    NGuidedMarkers = size(GuidedMarkers,1); % markers defined in _guidedVars file
                    NMarkers = size(R.Subject.Markers,1);   % markers in  the model


                    if NMarkers ~= NGuidedMarkers
                        ModelFilename = getFilenameAndExt(R.Model.File);
                        ModelFile = [ModelFilename,'.xml'];
                        GuidedVarsFile = [R.Model.Guided,'.m'];
                        
                        error(['These files have a different number of markers:\n'...
                               '     ',ModelFile,' defines ',num2str(NMarkers),' markers\n'...
                               '     ',GuidedVarsFile,' defines ',num2str(NGuidedMarkers),' markers\n'...
                               '   The number of markers should be the same']);
                    end
                    for i=1:NGuidedMarkers
                        SegmentIndex = getVecIndex(GuidedMarkers{i,1},R.Subject.Segments);
                        if SegmentIndex > 0
                            NSegMarker = size(R.Subject.Segments(SegmentIndex).LocalMarkers);
                            for j=1:NSegMarker
                                R.Subject.Segments(SegmentIndex).LocalMarkers(j).Point.Wm = GuidedMarkers{i,2};
                                R.Markers = [R.Markers;R.Subject.Segments(SegmentIndex).LocalMarkers(j).Point];
                            end
                        else
                            MarkerIndex = getVecIndex(GuidedMarkers{i,1},R.Subject.Markers);
                            if MarkerIndex > 0
                                R.Subject.Markers(MarkerIndex).Wm = GuidedMarkers{i,2};
                                if size(GuidedMarkers,2)== 3
                                    R.Subject.Markers(MarkerIndex).MarkerTrajectory = GuidedMarkers{i,3};
                                end
                                R.Markers = [R.Markers;R.Subject.Markers(MarkerIndex)];
                            else 
                                error('The %s Segment or Marker does not belong to the model',GuidedMarkers{i,1});
                            end
                        end
                    end
                end
            end
        end
        function calcSegmentForces(R,DataCSV,MarkerNames,MarkerCoords,Sensor,SensorName,SegmentName,Frequency)
            MarkersIndex = [];
            MissingMarkers = {};
            Glob_Pos_Markers = [];
            Sen1Fix_Pos_MarkerSensors = [];
            Glob_Pos_MarkerSensor = [];
            if strcmp(SensorName,'Sensor1')
                % Guardo los indices y las coordenadas tanto en locales en el sistema de los markes fijos entre experimentos como las globales
                % para dichos marcadores
                SenMarNames = fieldnames(Sensor);
                for i=1:8
                    Index = getCellIndex(MarkerNames,Sensor.(SenMarNames{i}).Label);
                    SumNaN = sum(isnan(MarkerCoords(:,3*(Index-1)+1)));
                    if isempty(Index) || SumNaN == size(MarkerCoords,1)
                        str =(['  WARNING: The Marker ',Sensor.(SenMarNames{i}).Label,' is missing all the time in motion file.']);
                        fprintf(R.ExperLogFileId, '%s\n', str);
                        MissingMarkers = [MissingMarkers;Sensor.(SenMarNames{i}).Label];
                    else
                        MarkersIndex = [MarkersIndex; Index];
                        if isfield(Sensor.(SenMarNames{i}),'LocCoord')
                            Sen1Fix_Pos_MarkerSensors = [Sen1Fix_Pos_MarkerSensors,Sensor.(SenMarNames{i}).LocCoord];
                            MarkerSensor_i = MarkerCoords(:,3*(Index-1)+1:3*(Index-1)+3);
                            Glob_Pos_MarkerSensor = [Glob_Pos_MarkerSensor,MarkerSensor_i];
                        end
                    end
                end
                % Utilidad para dibujar los marcadores del sensor
                NMarkers = size(MarkersIndex,1);
                for i=1:NMarkers%###
                    MarkerSensor_i = MarkerCoords(:,3*(MarkersIndex(i)-1)+1:3*(MarkersIndex(i)-1)+3);
                    R.MarkerSensor = [R.MarkerSensor,MarkerSensor_i];
                   
                end
                for i=1:8-NMarkers
                    MarkerSensor_i = NaN(size(MarkerCoords,1),3);
                    R.MarkerSensor = [R.MarkerSensor,MarkerSensor_i];
                end
                % Calculo la R y la d en el sistema con los marcadores fijos entre experimentos por Söderkvist para todo el
                % movimiento
                SegIndex = getVecIndex(SegmentName,R.Subject.Segments);
%                 R.Subject.Segments(SegIndex).F_Ext.Sys = [];
%                 R.Subject.Segments(SegIndex).F_Ext.Value = [];
%                 R.Subject.Segments(SegIndex).F_Ext.Pos   = [];
                Glob_R_Sen1Fix = [];
                if NMarkers<4
                        str=('  WARNING: Imposible definition of Sensor 1 LCS. 5 to 8 markers missing all the time in motion file');
                        disp(str)
                        fprintf(R.ExperLogFileId, '%s\n', str);
                        R.ErrorID = 1;
                        R.MarkerSensor = [];
                else
                    for j=1:size(Glob_Pos_MarkerSensor,1)
                        Glob_Pos_MarkerSensors_j = [];
                        Sen1Fix_Pos_MarkerSensors_j = [];
                        for k=1:(size(Glob_Pos_MarkerSensor,2))/3
                            if ~isnan(Glob_Pos_MarkerSensor(j,3*k-2:3*k))
                                Glob_Pos_MarkerSensors_j = [Glob_Pos_MarkerSensors_j,Glob_Pos_MarkerSensor(j,3*k-2:3*k)'];
                                Sen1Fix_Pos_MarkerSensors_j = [Sen1Fix_Pos_MarkerSensors_j,Sen1Fix_Pos_MarkerSensors(:,k)];
                            end
                        end
                        [Glob_R_Sen1Fix_j, Glob_Pos_OrSen1Fix(:,j)] = calcOptPose(Glob_Pos_MarkerSensors_j, Sen1Fix_Pos_MarkerSensors_j);
                        Glob_R_Sen1Fix = [Glob_R_Sen1Fix;Glob_R_Sen1Fix_j];
                    end
                end
                % Miro si hay un marcador perdido durante todo el experimento
                if(NMarkers<8)&& R.ErrorID == 0
                    MissXM = 0;
                    MissYM = 0;
                    for i=1:size(MissingMarkers)
                        if strcmpi('CPD_XM',MissingMarkers{i})
                            MissXM = 1;
                        elseif strcmpi('CPD_YM',MissingMarkers{i})
                            MissYM = 1;
                        end
                        % Calculo un valor en globales para el primer frame de los marcadores que están fijos.
                        if strcmpi('CPD_XP',MissingMarkers{i})
                            Glob_Pos_XP = Glob_Pos_OrSen1Fix(:,1)+ Glob_R_Sen1Fix(1:3,:)*Sensor.Marker_X.LocCoord;
                        elseif strcmpi('CPD_YF',MissingMarkers{i})
                            Glob_Pos_YF = Glob_Pos_OrSen1Fix(:,1)+ Glob_R_Sen1Fix(1:3,:)*Sensor.Marker_Y_Neg.LocCoord;
                        elseif strcmpi('CPD_YB',MissingMarkers{i})
                            Glob_Pos_YB = Glob_Pos_OrSen1Fix(:,1)+ Glob_R_Sen1Fix(1:3,:)*Sensor.Marker_Y_Pos.LocCoord;
                        elseif strcmpi('CPD_ZP',MissingMarkers{i})
                            Glob_Pos_ZP = Glob_Pos_OrSen1Fix(:,1)+ Glob_R_Sen1Fix(1:3,:)*Sensor.Marker_Z.LocCoord;
                        elseif strcmpi('CPD_FU',MissingMarkers{i})
                            Glob_Pos_FU = Glob_Pos_OrSen1Fix(:,1)+ Glob_R_Sen1Fix(1:3,:)*Sensor.Marker_Rot1.LocCoord;
                        elseif strcmpi('CPD_BU',MissingMarkers{i})
                            Glob_Pos_BU = Glob_Pos_OrSen1Fix(:,1)+ Glob_R_Sen1Fix(1:3,:)*Sensor.Marker_Rot2.LocCoord;
                        end
                    end
                    % Si XM e YM están perdidos no se puede resolver
                    if MissXM == 1 && MissYM == 1
                        str=('  WARNING: Imposible definition of Sensor 1 LCS. CPD_XM and CPD_YM missing all the time in motion file');
                        disp(str)
                        fprintf(R.ExperLogFileId, '%s\n', str);
                        R.ErrorID = 1;
                        R.MarkerSensor = [];
                    % Si sólo uno de los dos está perdido, hallo un valor en globales para el primer frame.
                    % Calculo la local en el SC fijo del otro marcador, con este y la distancia entre ambos cte para todos los movimientos calculo su local y luego la global. 
                    elseif MissXM == 1 && MissYM == 0
                        IndexYM =  getCellIndex(MarkerNames,Sensor.Marker_Y_Axis_Neg.Label);
                        for j=1:size(MarkerCoords,1)
                            Glob_Pos_YM = MarkerCoords(j,3*(IndexYM-1)+1:3*(IndexYM-1)+3);
                            if ~isnan(Glob_Pos_YM)
                                Sen1Fix_Pos_YM = Glob_R_Sen1Fix(3*j-2:3*j,:)' * (Glob_Pos_YM' - Glob_Pos_OrSen1Fix(:,j));
                                break;
                            end
                        end
                        Sen1Fix_Pos_XM = Sen1Fix_Pos_YM + Sensor.DifMarkers.LocCoord;
                        Glob_Pos_XM = (Glob_Pos_OrSen1Fix (:,1)+ Glob_R_Sen1Fix(1:3,:)*Sen1Fix_Pos_XM);
                        
                    elseif MissYM == 1 && MissXM == 0
                        IndexXM =  getCellIndex(MarkerNames,Sensor.Marker_X_Axis_Neg.Label);
                        for j=1:size(MarkerCoords,1)
                            Glob_Pos_XM = MarkerCoords(j,3*(IndexXM-1)+1:3*(IndexXM-1)+3);
                            if ~isnan(Glob_Pos_XM)
                                Sen1Fix_Pos_XM = Glob_R_Sen1Fix(3*j-2:3*j,:)' * (Glob_Pos_XM' - Glob_Pos_OrSen1Fix(:,j));
                                break;
                            end
                        end
                        Sen1Fix_Pos_YM = Sen1Fix_Pos_XM - Sensor.DifMarkers.LocCoord;
                        Glob_Pos_YM = Glob_Pos_OrSen1Fix (:,1)+ Glob_R_Sen1Fix(1:3,:)*Sen1Fix_Pos_YM;
                        
                    end
                    
                end
%                 % Utilizamos los 20 primeros frames media
%                 Glob_Pos_Markers = [];
%                 for i=1:NMarkers
%                     Glob_Pos_Marker_i = MarkerCoords(1:20,3*(MarkersIndex(i)-1)+1:3*(MarkersIndex(i)-1)+3);
%                     Glob_Pos_Marker_i = sum(Glob_Pos_Marker_i)/size(Glob_Pos_Marker_i,1);
%                     Glob_Pos_Markers =  [Glob_Pos_Markers,(Glob_Pos_Marker_i)'];
%                 end
%                 % Utilizamos los 20 primeros frames mediana
%                 Glob_Pos_Markers = [];
%                 for i=1:NMarkers
%                     Glob_Pos_Marker_i = MarkerCoords(1:20,3*(MarkersIndex(i)-1)+1:3*(MarkersIndex(i)-1)+3);
%                     Glob_Pos_Marker_i = median(Glob_Pos_Marker_i);
%                     Glob_Pos_Markers =  [Glob_Pos_Markers,(Glob_Pos_Marker_i)'];
%                 end
                
%                 if R.ErrorID == 0
                    % Encontrar un frame donde se ven todos los markers
%                     while (size(Glob_Pos_Markers,2) ~= 6 && j~=size(MarkerCoords,1))
%                         Glob_Pos_Markers = [];
%                         j = j+1;
%                         for i=1:NMarkers
%                             Glob_Pos_Marker_i = MarkerCoords(j,3*(MarkersIndex(i)-1)+1:3*(MarkersIndex(i)-1)+3);
%                             if isnan(Glob_Pos_Marker_i)
%                                 break;
%                             else
%                                 Glob_Pos_Markers =  [Glob_Pos_Markers,(Glob_Pos_Marker_i)'];
%                             end
%                         end
%                         if j==size(MarkerCoords,1)
%                             str=('  WARNING: Imposible definition of Sensor 1 LCS. Some markers Nan all the time');
%                             disp(str)
%                             fprintf(R.ExperLogFileId, '%s\n', str);
%                             R.ErrorID = 1;
%                         end
%                     end
                    % Obtengo un valor para el primer frame de cada marcador, siguiendo el siguiente orden
                    % 1º global que marca el csv.
                    % 2º si el primer frame está perdido, mediante su coordenada local calculo la global
                    % 3º si está todo el rato perdido, el calculado anteriormente
                    if R.ErrorID == 0
                        IndexXP =  getCellIndex(MarkerNames,Sensor.Marker_X.Label);
                        if isempty(IndexXP)
                            Glob_Pos_XP = Glob_Pos_XP;
                        else
                            if isnan(MarkerCoords(1,3*(IndexXP-1)+1:3*(IndexXP-1)+3))
                                Glob_Pos_XP = Glob_Pos_OrSen1Fix (:,1)+ Glob_R_Sen1Fix(1:3,:)*Sensor.Marker_X.LocCoord;
                            else
                                Glob_Pos_XP = MarkerCoords(1,3*(IndexXP-1)+1:3*(IndexXP-1)+3)';
                            end
                        end
                        IndexYF =getCellIndex(MarkerNames,Sensor.Marker_Y_Neg.Label);
                        if isempty(IndexYF)
                            Glob_Pos_YF = Glob_Pos_YF;
                        else
                            if isnan(MarkerCoords(1,3*(IndexYF-1)+1:3*(IndexYF-1)+3))
                                Glob_Pos_YF = Glob_Pos_OrSen1Fix (:,1)+ Glob_R_Sen1Fix(1:3,:)*Sensor.Marker_Y_Neg.LocCoord;
                            else
                                Glob_Pos_YF = MarkerCoords(1,3*(IndexYF-1)+1:3*(IndexYF-1)+3)';
                            end
                        end
                        IndexYB =getCellIndex(MarkerNames,Sensor.Marker_Y_Pos.Label);
                        if isempty(IndexYB)
                            Glob_Pos_YB = Glob_Pos_YB;
                        else
                            if isnan(MarkerCoords(1,3*(IndexYB-1)+1:3*(IndexYB-1)+3))
                                Glob_Pos_YB = Glob_Pos_OrSen1Fix (:,1)+ Glob_R_Sen1Fix(1:3,:)*Sensor.Marker_Y_Pos.LocCoord;
                            else
                                Glob_Pos_YB = MarkerCoords(1,3*(IndexYB-1)+1:3*(IndexYB-1)+3)';
                            end
                        end
                        IndexZP =getCellIndex(MarkerNames,Sensor.Marker_Z.Label);
                        if isempty(IndexZP)
                            Glob_Pos_ZP = Glob_Pos_ZP;
                        else
                            if isnan(MarkerCoords(1,3*(IndexZP-1)+1:3*(IndexZP-1)+3))
                                Glob_Pos_ZP = Glob_Pos_OrSen1Fix (:,1)+ Glob_R_Sen1Fix(1:3,:)*Sensor.Marker_Z.LocCoord;
                            else
                                Glob_Pos_ZP = MarkerCoords(1,3*(IndexZP-1)+1:3*(IndexZP-1)+3)';
                            end
                        end
                        IndexFU =getCellIndex(MarkerNames,Sensor.Marker_Rot1.Label);
                        if isempty(IndexFU)
                            Glob_Pos_FU = Glob_Pos_FU;
                        else
                            if isnan(MarkerCoords(1,3*(IndexFU-1)+1:3*(IndexFU-1)+3))
                                Glob_Pos_FU = Glob_Pos_OrSen1Fix (:,1)+ Glob_R_Sen1Fix(1:3,:)*Sensor.Marker_Z.LocCoord;
                            else
                                Glob_Pos_FU = MarkerCoords(1,3*(IndexFU-1)+1:3*(IndexFU-1)+3)';
                            end
                        end
                        IndexBU =getCellIndex(MarkerNames,Sensor.Marker_Rot2.Label);
                        if isempty(IndexBU)
                            Glob_Pos_BU = Glob_Pos_BU;
                        else
                            if isnan(MarkerCoords(1,3*(IndexBU-1)+1:3*(IndexBU-1)+3))
                                Glob_Pos_BU = Glob_Pos_OrSen1Fix (:,1)+ Glob_R_Sen1Fix(1:3,:)*Sensor.Marker_Z.LocCoord;
                            else
                                Glob_Pos_BU = MarkerCoords(1,3*(IndexBU-1)+1:3*(IndexBU-1)+3)';
                            end
                        end
                        IndexXM =getCellIndex(MarkerNames,Sensor.Marker_X_Axis_Neg.Label);
                        if isempty(IndexXM)
                            Glob_Pos_XM = Glob_Pos_XM;
                        else
                            if isnan(MarkerCoords(1,3*(IndexXM-1)+1:3*(IndexXM-1)+3))
                                 
                                 for j=1:size(MarkerCoords,1)
                                     Glob_Pos_XM = MarkerCoords(j,3*(IndexXM-1)+1:3*(IndexXM-1)+3);
                                     if ~isnan(Glob_Pos_XM)
                                         Sen1Fix_Pos_XM = Glob_R_Sen1Fix(3*j-2:3*j,:)' * (Glob_Pos_XM' - Glob_Pos_OrSen1Fix(:,j));
                                         break;
                                     end
                                 end
                                Glob_Pos_XM = Glob_Pos_OrSen1Fix (:,1)+ Glob_R_Sen1Fix(1:3,:)*Sen1Fix_Pos_XM;
                            else
                                Glob_Pos_XM = MarkerCoords(1,3*(IndexXM-1)+1:3*(IndexXM-1)+3)';
                            end
                        end
                        IndexYM =getCellIndex(MarkerNames,Sensor.Marker_Y_Axis_Neg.Label);
                        if isempty(IndexYM)
                            Glob_Pos_YM = Glob_Pos_YM;
                        else
                            if isnan(MarkerCoords(1,3*(IndexYM-1)+1:3*(IndexYM-1)+3))
                                 for j=1:size(MarkerCoords,1)
                                     Glob_Pos_YM = MarkerCoords(j,3*(IndexYM-1)+1:3*(IndexYM-1)+3);
                                     if ~isnan(Glob_Pos_YM)
                                         Sen1Fix_Pos_YM = Glob_R_Sen1Fix(3*j-2:3*j,:)' * (Glob_Pos_YM' - Glob_Pos_OrSen1Fix(:,j));
                                         break;
                                     end
                                 end
                                Glob_Pos_YM = Glob_Pos_OrSen1Fix (:,1)+ Glob_R_Sen1Fix(1:3,:)*Sen1Fix_Pos_YM;
                            else
                                Glob_Pos_YM = MarkerCoords(1,3*(IndexYM-1)+1:3*(IndexYM-1)+3)';
                            end
                        end
                        % Calucular el sistema local del Sensor
                        Glob_Pos_Markers = [Glob_Pos_XP,Glob_Pos_YF,Glob_Pos_YB,Glob_Pos_ZP,Glob_Pos_FU,Glob_Pos_BU,Glob_Pos_XM,Glob_Pos_YM];
                        [Glob_R_Sensor1,Glob_Pos_OrSensor1]= R.calcSensor1LCS(Glob_Pos_Markers);
                        % Obtener las coordenadas locales de los marcadores en el sensor
                        for i=1:size(Glob_Pos_Markers,2)
                            Sensor1_Pos_Markeri(:,i) = Glob_R_Sensor1' * (Glob_Pos_Markers(:,i) - Glob_Pos_OrSensor1);
                        end
                        % Obtener la posición de aplicación de la fuerza y el valor de la fuerza del sistema en globales
                        Glob_Pos_Markers = [];
                        Sen1_Pos_Markers = [];
                        % check Number of Fext in segment
                        if isfield(R.Subject.Segments(SegIndex).F_Ext,'Sys')
                            NFExt = size(R.Subject.Segments(SegIndex).F_Ext,1) + 1;
                        else
                            NFExt = 1;
                        end
                        for i=1:size(MarkerCoords,1)
                            for j=1:NMarkers
                                Glob_Pos_Marker_i = MarkerCoords(i,3*(MarkersIndex(j)-1)+1:3*(MarkersIndex(j)-1)+3);
                                if ~isnan(Glob_Pos_Marker_i)
                                    Glob_Pos_Markers = [Glob_Pos_Markers,(Glob_Pos_Marker_i)'];
                                    Sen1_Pos_Markers = [Sen1_Pos_Markers,Sensor1_Pos_Markeri(:,j)];
                                end
                            end
                            [Glob_R_Sen1, Glob_Pos_OrSen1] = calcOptPose(Glob_Pos_Markers, Sen1_Pos_Markers);
                            Sen1_For_Fx = DataCSV.Analog.(Sensor.Label_Fx).Values((DataCSV.Analog.Frequency/Frequency)*i,1); 
                            Sen1_For_Fy = DataCSV.Analog.(Sensor.Label_Fy).Values((DataCSV.Analog.Frequency/Frequency)*i,1);
                            Sen1_For_Fz = DataCSV.Analog.(Sensor.Label_Fz).Values((DataCSV.Analog.Frequency/Frequency)*i,1);
                            Glob_Forces = Glob_R_Sen1 * [Sen1_For_Fx;Sen1_For_Fy;Sen1_For_Fz];
                            R.Subject.Segments(SegIndex).F_Ext(NFExt).Sys = 1;
                            R.Subject.Segments(SegIndex).F_Ext(NFExt).Value(i,:) = Glob_Forces;
                            R.Subject.Segments(SegIndex).F_Ext(NFExt).Pos(i,:)   = Glob_Pos_OrSen1;
                            Glob_Pos_Markers = [];
                            Sen1_Pos_Markers = [];
                        end
                        %                 tiempoTotal = 765;%size(MarkerCoords,1);
                        %                 for k =1:tiempoTotal
                        %                     tiempo(k)=k/100;
                        %                 end
                        %                 figure('Name','Reducido F2z','NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
                        %                 % figure settings
                        %                 grid on
                        %                 plot(tiempo,R.Subject.Segments(SegIndex).F_Ext.Value(:,3));
                    end
%                 end
            elseif strcmp(SensorName,'SeatBack')
                % set data
                NFrames = size(DataCSV.Back,2);
                SegIndex = getVecIndex(SegmentName,R.Subject.Segments);
                % set data from the sensor
                IndexDirXP1 = getCellIndex(MarkerNames,Sensor.DirX.P1Name);
                IndexDirXP2 = getCellIndex(MarkerNames,Sensor.DirX.P2Name);
                IndexDirXP3 = getCellIndex(MarkerNames,Sensor.DirX.P3Name);
                IndexDirXP4 = getCellIndex(MarkerNames,Sensor.DirX.P4Name);
                IndexDirYP1 = getCellIndex(MarkerNames,Sensor.DirY.P1Name);
                IndexDirYP2 = getCellIndex(MarkerNames,Sensor.DirY.P2Name);
                CoordDirXP1 = MarkerCoords(1,3*(IndexDirXP1-1)+1:3*(IndexDirXP1-1)+3)';
                CoordDirXP2 = MarkerCoords(1,3*(IndexDirXP2-1)+1:3*(IndexDirXP2-1)+3)';
                CoordDirXP3 = MarkerCoords(1,3*(IndexDirXP3-1)+1:3*(IndexDirXP3-1)+3)';
                CoordDirXP4 = MarkerCoords(1,3*(IndexDirXP4-1)+1:3*(IndexDirXP4-1)+3)';
                CoordDirYP1 = MarkerCoords(1,3*(IndexDirYP1-1)+1:3*(IndexDirYP1-1)+3)';
                CoordDirYP2 = MarkerCoords(1,3*(IndexDirYP2-1)+1:3*(IndexDirYP2-1)+3)';
                % Norm in the BackRest
                UpPoint  = (CoordDirXP1+CoordDirXP2)/2;
                DownPoint = (CoordDirXP3+CoordDirXP4)/2;
                AxeX = (UpPoint - DownPoint)/norm(UpPoint - DownPoint);
                AxeX = rot_y(Sensor.DirX.AngleY)*AxeX;
                AxeY = CoordDirYP1 - CoordDirYP2/norm(CoordDirYP1 - CoordDirYP2);
                Norm = cross(AxeX,AxeY)/norm(cross(AxeX,AxeY));
                % Sensor Initial position
                IndexAplP = getCellIndex(MarkerNames,Sensor.AplPoint.Name);
                CoordAplP = MarkerCoords(1,3*(IndexAplP-1)+1:3*(IndexAplP-1)+3)';
                SenIni = CoordAplP + Sensor.AplPoint.Coord;
                % check Number of Fext in segment
                if isfield(R.Subject.Segments(SegIndex).F_Ext,'Sys')
                    NFExt = size(R.Subject.Segments(SegIndex).F_Ext,1) + 1;
                else
                    NFExt = 1;
                end
                for i=1:NFrames
                    if strcmpi(DataCSV.Units,'PSI')
                        F_Value = sum(sum(DataCSV.Back(i).Values))*0.5*0.5*4.44822162; % squares are 05x05 inch/ 4.44822162 from pound to N.
                        PosX = DataCSV.Back(i).CoPRow*0.5*0.0254; % 1 inch 0.0254m
                        PosY = (DataCSV.Back(i).CoPCol - 24)*0.5*0.0254;
                    elseif strcmpi(DataCSV.Units,'N/cm2')
                        SumP = sum(sum(DataCSV.Back(i).Values));
                        SumPRow = sum(DataCSV.Back(i).Values,2);
                        SumPCol = sum(DataCSV.Back(i).Values,1);
                        F_Value = sum(sum(DataCSV.Back(i).Values))*0.5*0.5*2.54*2.54; % squares are 05x05 inch/ 1 inch 2.54cm.
                        PRow_x = 0;
                        PCol_y = 0;
                        for j=1:size(SumPRow,1)
                            PRow_x = PRow_x + SumPRow(j)*j;
                        end
                        for j=1:size(SumPCol,2)
                            PCol_y = PCol_y + SumPCol(j)*j;
                        end
                        CoPRow = PRow_x/SumP;
                        CoPCol = PCol_y/SumP;
                        PosX = CoPRow*0.5*0.0254; % 1 inch 0.0254m
                        PosY = (CoPCol - 24)*0.5*0.0254;
                    else
                        error(['The units ',DataCSV.Units,' of the pressure maps are not correct.'])
                    end
                    Glob_Pos = SenIni + AxeX*PosX - AxeY*PosY;
                    Glob_F = F_Value*Norm;
                    R.Subject.Segments(SegIndex).F_Ext(NFExt).Sys = 1;
                    R.Subject.Segments(SegIndex).F_Ext(NFExt).Value(i,:) = Glob_F;
                    R.Subject.Segments(SegIndex).F_Ext(NFExt).Pos(i,:)   = Glob_Pos;
                end
            elseif strcmp(SensorName,'Seat')
                % set data
                NFrames = size(DataCSV.Back,2);
                SegIndex = getVecIndex(SegmentName,R.Subject.Segments);
                % set data from the sensor
                IndexDirXP1 = getCellIndex(MarkerNames,Sensor.DirX.P1Name);
                IndexDirXP2 = getCellIndex(MarkerNames,Sensor.DirX.P2Name);
                IndexDirYP1 = getCellIndex(MarkerNames,Sensor.DirY.P1Name);
                IndexDirYP2 = getCellIndex(MarkerNames,Sensor.DirY.P2Name);
                CoordDirXP1 = MarkerCoords(1,3*(IndexDirXP1-1)+1:3*(IndexDirXP1-1)+3)';
                CoordDirXP2 = MarkerCoords(1,3*(IndexDirXP2-1)+1:3*(IndexDirXP2-1)+3)';
                CoordDirYP1 = MarkerCoords(1,3*(IndexDirYP1-1)+1:3*(IndexDirYP1-1)+3)';
                CoordDirYP2 = MarkerCoords(1,3*(IndexDirYP2-1)+1:3*(IndexDirYP2-1)+3)';
                % Norm in the seat
                AxeX = (CoordDirXP1 - CoordDirXP2)/norm(CoordDirXP1 - CoordDirXP2);
                AxeX = rot_y(Sensor.DirX.AngleY)*AxeX;
                AxeY = (CoordDirYP1 - CoordDirYP2)/norm(CoordDirYP1 - CoordDirYP2);
                Norm = cross(AxeX,AxeY)/norm(cross(AxeX,AxeY));
                % Sensor Initial position
                IndexAplP = getCellIndex(MarkerNames,Sensor.AplPoint.Name);
                CoordAplP = MarkerCoords(1,3*(IndexAplP-1)+1:3*(IndexAplP-1)+3)';
                SenIni = CoordAplP + Sensor.AplPoint.Coord;
                % check Number of Fext in segment
                if isfield(R.Subject.Segments(SegIndex).F_Ext,'Sys')
                    NFExt = size(R.Subject.Segments(SegIndex).F_Ext,1) + 1;
                else
                    NFExt = 1;
                end
                for i=1:NFrames
                    ZeroRow = DataCSV.Seat(i).ZeroRow;
                    ZeroCol = DataCSV.Seat(i).ZeroCol;
                    if strcmpi(SegmentName,'BEC')||strcmpi(SegmentName,'HM50KR_posture')
                        SumP = sum(sum(DataCSV.Seat(i).Values(1:ZeroRow,:)));
                        SumPRow = sum(DataCSV.Seat(i).Values(1:ZeroRow,:),2); % with 2 sum function sum rows
                        SumPCol = sum(DataCSV.Seat(i).Values(1:ZeroRow,:),1); % with 1 sum function sum columns
                    elseif strcmpi(SegmentName,'OSL')||strcmpi(SegmentName,'Left_leg')
                        SumP = sum(sum(DataCSV.Seat(i).Values(ZeroRow+1:end,ZeroCol:end)));
                        SumPRow = sum(DataCSV.Seat(i).Values(ZeroRow+1:end,ZeroCol:end),2); % with 2 sum function sum rows
                        SumPCol = sum(DataCSV.Seat(i).Values(ZeroRow+1:end,ZeroCol:end),1); % with 1 sum function sum columns
                    end
                    if strcmpi(DataCSV.Units,'PSI') 
                        F_Value = SumP*0.5*0.5*4.44822162; % squares are 05x05 inch/ 4.44822162 from pound to N.
                    elseif strcmpi(DataCSV.Units,'N/cm2')
                        F_Value = SumP*0.5*0.5*2.54*2.54; % squares are 05x05 inch/ 1 inch 2.54cm.
                    else
                        error(['The units ',DataCSV.Units,' of the pressure maps are not correct.'])
                    end    
                    PRow_x = 0;
                    PCol_y = 0;
                    for j=1:size(SumPRow,1)
                        PRow_x = PRow_x + SumPRow(j)*j;
                    end
                    for j=1:size(SumPCol,2)
                        PCol_y = PCol_y + SumPCol(j)*j;
                    end
                    CoPRow = PRow_x/SumP;
                    CoPCol = PCol_y/SumP;
                    if strcmpi(SegmentName,'OSL')||strcmpi(SegmentName,'Left_leg')
                        CoPRow = CoPRow + ZeroRow;
                        CoPCol = CoPCol + ZeroCol;
                    end
                    PosX = CoPRow*0.5*0.0254; % 1 inch 0.0254m
                    PosY = (CoPCol - 24)*0.5*0.0254;
                    Glob_Pos = SenIni + AxeX*PosX + AxeY*PosY;
                    Glob_F = F_Value*Norm;
                    
                    R.Subject.Segments(SegIndex).F_Ext(NFExt).Sys = 1;
                    R.Subject.Segments(SegIndex).F_Ext(NFExt).Value(i,:) = Glob_F;
                    R.Subject.Segments(SegIndex).F_Ext(NFExt).Pos(i,:)   = Glob_Pos;
                end
            end
        end
        function [R, OX, OY, cosang] = calcSensor1LCS (R,Glob_Pos_Markers)
            
            % Definition of the local coordinate system for the pedal force sensor
            % Input:
            % YF, YB         : global coordinates of the two markers used to define
            %                  the Y direction
            % XP             : idem for the X direction
            % ZP             : idem for the Z direction, left handed system
            % XM, ZM :         global coordinates of the markers used to define
            %                  the origin of the LCS (centre of the force sensor)
            % Output:
            % R: rotation matrix of the sensor LCS with respect the GCS. Its columns
            % are the vectors of the LCS. LEFT HANDED
            % OX: global coordinates of the centre of the sensor, calculated using XM
            % OY: idem, calculated using YM
            % cosang: cosinus of the angle between the planes defined by XP and ZP
            
            YF = Glob_Pos_Markers(:,2);
            YB = Glob_Pos_Markers(:,3);
            XP = Glob_Pos_Markers(:,1);
            ZP = Glob_Pos_Markers(:,4);
            XM = Glob_Pos_Markers(:,7);
            YM = Glob_Pos_Markers(:,8);
            % Unit vector in the Y direction
            V=(YB-YF)/norm(YB-YF);

            
            % Unit vector in the Z direction 
%             W=-cross(XP-YF,V); %(left handed)
            W=cross(XP-YF,V); %(right handed)
            W=W/norm(W);
            
            % Unit vector in the X direction
%             U=cross(W,V); % Levogiro
            U=-cross(W,V); % Destrogiro
            
            % Alternate system
            U2=-cross(V,ZP-YF);% Levogiro
%             U2=cross(V,ZP-YF);% Destrogiro
            U2=U2/norm(U2);
            W2=-cross(U2,V); % Levogiro
%             W2=cross(U2,V); % Destrogiro 
            
            % Angle between both systems
            cosang=W2'*U;
            
            % Orientation matrix
            R=[U V W];
            
            % Origin of the sensor LCS
            % Intersection of lines (XM, U) and (YM, V)
            % Equation to solve: XM + U * l = YM +  V * m
            % In form:   U * l - V * m = YM - XM
            A=[U -V];
            B=YM-XM;
            lm=A\B;
            
            % Intersection point, calculated from the two lines
            OX=XM+lm(1)*U;
            OY=YM+lm(2)*V;
%             OX = YF;
%             OY = 0;
            
            % Distance between the two lines
            dist=det([A B]);
            
            % Dibujar para comprobar
%             figure('Name','Sistema coord','NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%             U=R(:,1);
%             V=R(:,2);
%             W=R(:,3);
% %             OR=0.5*(OX+OY);
%             OR = OX;
%             
%             A=[U -V];
%             B=YM-XM;
%             AB=[A B];
%             
%             % Dibujo
%             
%             % Sistema de ejes LCS
%             UF=OR+U;
%             VF=OR+V;
%             WF=OR+W;
%             line([OR(1) UF(1)],[OR(2) UF(2)],[OR(3) UF(3)],'Color','r');
%             line([OR(1) VF(1)],[OR(2) VF(2)],[OR(3) VF(3)],'Color','g');
%             line([OR(1) WF(1)],[OR(2) WF(2)],[OR(3) WF(3)],'Color','b');
%             
%             % Markers
%             line([YF(1) YF(1)],[YF(2) YF(2)],[YF(3) YF(3)],'Marker','O','Color','g')
%             line([YB(1) YB(1)],[YB(2) YB(2)],[YB(3) YB(3)],'Marker','O','Color','g')
%             line([XP(1) XP(1)],[XP(2) XP(2)],[XP(3) XP(3)],'Marker','O','Color','r')
%             line([ZP(1) ZP(1)],[ZP(2) ZP(2)],[ZP(3) ZP(3)],'Marker','O','Color','b')
%             
%             line([XM(1) OX(1)],[XM(2) OX(2)],[XM(3) OX(3)],'Marker','x','Color','k')
%             line([YM(1) OY(1)],[YM(2) OY(2)],[YM(3) OY(3)],'Marker','O','Color','k')
%             
%             % Ejes generales
%             line([0 10],[0 0],[0 0],'Color','r');
%             line([0 0],[0 10],[0 0],'Color','g');
%             line([0 0],[0 0],[0 10],'Color','b');
%             view(120,15);
%             axis equal;
            
        end
        function calcUpperBodyInertialParameters(R)
            global PathBar
            % Calc of total mass, heigh and gender from .alm
            IntSubFileDOM = xmlread([R.ExperimentPath,'Subjects',PathBar,R.Subject.SubjectName,PathBar,R.Subject.SubjectName,'.xml']);
            SubjectMass = str2num(IntSubFileDOM.getElementsByTagName('AM25').item(0).getAttribute('Value'));% kg
            SubjectHeigh = str2num(IntSubFileDOM.getElementsByTagName('AM01').item(0).getAttribute('Value'))/10;% cm
            SubjectGender = IntSubFileDOM.getElementsByTagName('AM33').item(0).getAttribute('Value');
            % Calc the CoM and I of each segment
            if strcmpi(SubjectGender,'MALE')
                HandMass = -0.1165 + 0.0036*SubjectMass + 0.00175*SubjectHeigh;
                ForearmMass = 0.3185 + 0.01445*SubjectMass - 0.00114*SubjectHeigh;
                UpperarmMass = 0.250 + 0.03012*SubjectMass - 0.0027*SubjectHeigh;
                HeadNeckMass = 1.296 + 0.0171*SubjectMass + 0.0143*SubjectHeigh;
                UpTrunkMass = 8.2144 + 0.1862*SubjectMass - 0.0584*SubjectHeigh;
                MidTrunkMass = 7.181 + 0.2234*SubjectMass - 0.0663*SubjectHeigh;
                HandCoM = (4.11 + 0.026*SubjectMass + 0.033*SubjectHeigh)*0.01; % cm to m
                ForearmCoM = (0.192 - 0.028*SubjectMass + 0.093*SubjectHeigh)*0.01;
                UpperarmCoM = (1.67 + 0.03*SubjectMass + 0.054*SubjectHeigh)*0.01;
                HeadNeckCoM = (8.357 - 0.0025*SubjectMass + 0.023*SubjectHeigh)*0.01;
                UpTrunkCoM = (3.32 + 0.0076*SubjectMass + 0.047*SubjectHeigh)*0.01;
                MidTrunkCoM = (1.398 + 0.0058*SubjectMass + 0.045*SubjectHeigh)*0.01;
                HandIAntPos = -19.5 + 0.17*SubjectMass + 0.116*SubjectHeigh;
                ForearmIAntPos = -64 + 0.95*SubjectMass + 0.34*SubjectHeigh;
                UpperarmIAntPos = -250.7 + 1.56*SubjectMass + 1.512*SubjectHeigh;
                HeadNeckIAntPos = -78 + 1.171*SubjectMass + 1.519*SubjectHeigh;
                UpTrunkIAntPos = 81.2 + 36.73*SubjectMass - 5.97*SubjectHeigh;
                MidTrunkIAntPos = 618.5 + 39.8*SubjectMass - 12.87*SubjectHeigh;
                HandIMedLat = -13.68 + 0.088*SubjectMass + 0.092*SubjectHeigh;
                ForearmIMedLat = -67.9 + 0.855*SubjectMass + 0.376*SubjectHeigh;
                UpperarmIMedLat = -232 + 1.525*SubjectMass + 1.343*SubjectHeigh;
                HeadNeckIMedLat = -112 + 1.43*SubjectMass + 1.73*SubjectHeigh;
                UpTrunkIMedLat = 367 + 18.3*SubjectMass - 5.73*SubjectHeigh;
                MidTrunkIMedLat = 263 + 26.7*SubjectMass - 8.0*SubjectHeigh;
                HandILong = -6.26 + 0.0762*SubjectMass + 0.0347*SubjectHeigh;
                ForearmILong = 5.66 + 0.306*SubjectMass - 0.088*SubjectHeigh;
                UpperarmILong = -16.9 + 0.662*SubjectMass + 0.0435*SubjectHeigh;
                HeadNeckILong = 61.6 + 1.72*SubjectMass + 0.0814*SubjectHeigh;
                UpTrunkILong = 561 + 36.03*SubjectMass - 9.98*SubjectHeigh;
                MidTrunkILong = 1501 + 43.14*SubjectMass - 19.8*SubjectHeigh;
            elseif strcmpi(SubjectGender,'FEMALE')
                HandMass = -0.116 + 0.0017*SubjectMass + 0.002*SubjectHeigh;
                ForearmMass = 0.295 + 0.009*SubjectMass + 0.0003*SubjectHeigh;
                UpperarmMass = 0.206 + 0.0053*SubjectMass + 0.0066*SubjectHeigh;
                HeadNeckMass = 2.388 - 0.001*SubjectMass + 0.015*SubjectHeigh;
                UpTrunkMass = -16.593 + 0.14*SubjectMass + 0.0995*SubjectHeigh;
                MidTrunkMass = -2.741 + 0.031*SubjectMass + 0.056*SubjectHeigh;
                HandCoM = (41.74 -0.120*SubjectMass + 0.172*SubjectHeigh)*0.01; % cm to m
                ForearmCoM = (61.40 + 0.096*SubjectMass - 0.062*SubjectHeigh)*0.01; 
                UpperarmCoM = (44.96 + 0.034*SubjectMass + 0.051*SubjectHeigh)*0.01;
                HeadNeckCoM = (21.50 + 0.181*SubjectMass - 0.085*SubjectHeigh)*0.01;
                UpTrunkCoM = (34.5 + 0.012*SubjectMass + 0.084*SubjectHeigh)*0.01;
                MidTrunkCoM = (36.68 + 0.025*SubjectMass + 0.037*SubjectHeigh)*0.01;
                HandIAntPos = -5.71 + 0.122*SubjectMass + 0.035*SubjectHeigh;
                ForearmIAntPos = -132.1 + 0.62*SubjectMass + 0.825*SubjectHeigh;
                UpperarmIAntPos = -151.4 + 0.107*SubjectMass + 1.554*SubjectHeigh;
                HeadNeckIAntPos = 217.8 -0.032*SubjectMass + 0.059*SubjectHeigh;
                UpTrunkIAntPos = -4038.5 + 28.6*SubjectMass + 20*SubjectHeigh;
                MidTrunkIAntPos = -368.7 -6.22*SubjectMass + 8.86*SubjectHeigh;
                HandIMedLat = -5.79 + 0.087*SubjectMass + 0.034*SubjectHeigh;
                ForearmIMedLat = -138.5 + 0.533*SubjectMass + 0.887*SubjectHeigh;
                UpperarmIMedLat = -330.4 - 0.461*SubjectMass + 2.67*SubjectHeigh;
                HeadNeckIMedLat = 66.4 - 0.447*SubjectMass + 1.29*SubjectHeigh;
                UpTrunkIMedLat = -2075 + 15.6*SubjectMass + 9.4*SubjectHeigh;
                MidTrunkIMedLat = -546 + 2.87*SubjectMass + 5.1*SubjectHeigh;
                HandILong = -2.138 + 0.053*SubjectMass + 0.0073*SubjectHeigh;
                ForearmILong = 7.4 + 0.21*SubjectMass - 0.08*SubjectHeigh;
                UpperarmILong = -118.6 + 1.19*SubjectMass + 0.44*SubjectHeigh;
                HeadNeckILong = -35.48 + 2.43*SubjectMass + 0.237*SubjectHeigh;
                UpTrunkILong = -2823.2 + 25.8*SubjectMass + 12.8*SubjectHeigh;
                MidTrunkILong = -672.9 + 1.47*SubjectMass + 7.53*SubjectHeigh;
            end
            % Calc the position of the joints
            % get the Glob_Pos_GLK & Glob_R_BEC
            SegIndex = getVecIndex('BEC',R.Subject.Segments);
            GLKIndex = getVecIndex('GLK',R.Subject.Segments(SegIndex).LocalPoints);
            GLKPosInq = R.Subject.Segments(SegIndex).LocalPoints(GLKIndex).Point.PosInq;
            Glob_Pos_GLK = R.q_t(1,GLKPosInq:2+GLKPosInq)';
            Glob_R_BEC = R.Subject.Segments(SegIndex).getRd(R.q_t(1,:));
            % local positions of Joints
            % Angles from EXP
            deg2rad = pi/180;
            ParserExp = HUMAN_PARSER_EXP();
            ParserExp.readExp([R.ExperimentPath,'Subjects',PathBar,R.Subject.SubjectName,PathBar],[R.Subject.SubjectName,'.exp']);
            Joint.GHZ_a3 = ParserExp.StructExp.Jnt.GHZ.rz*deg2rad;
            Joint.GHZ_a2 = ParserExp.StructExp.Jnt.GHZ.ry*deg2rad;
            Joint.GHZ_a1 = ParserExp.StructExp.Jnt.GHZ.rx*deg2rad;
            Joint.GLK_a3 = ParserExp.StructExp.Jnt.GLK.rz*deg2rad;
            Joint.GLK_a2 = ParserExp.StructExp.Jnt.GLK.ry*deg2rad;
            Joint.GLK_a1 = ParserExp.StructExp.Jnt.GLK.rx*deg2rad;
            Joint.GLL_a3 = ParserExp.StructExp.Jnt.GLL.rz*deg2rad;
            Joint.GLL_a2 = ParserExp.StructExp.Jnt.GLL.ry*deg2rad;
            Joint.GLL_a1 = ParserExp.StructExp.Jnt.GLL.rx*deg2rad;
            Joint.GBL_a3 = ParserExp.StructExp.Jnt.GBL.rz*deg2rad;
            Joint.GBL_a2 = ParserExp.StructExp.Jnt.GBL.ry*deg2rad;
            Joint.GBL_a1 = ParserExp.StructExp.Jnt.GBL.rx*deg2rad;
            Joint.GBB_a3 = ParserExp.StructExp.Jnt.GBB.rz*deg2rad;
            Joint.GBB_a2 = ParserExp.StructExp.Jnt.GBB.ry*deg2rad;
            Joint.GBB_a1 = ParserExp.StructExp.Jnt.GBB.rx*deg2rad;
            Joint.GHB_a3 = ParserExp.StructExp.Jnt.GHB.rz*deg2rad;
            Joint.GHB_a2 = ParserExp.StructExp.Jnt.GHB.ry*deg2rad;
            Joint.GHB_a1 = ParserExp.StructExp.Jnt.GHB.rx*deg2rad;
            Joint.GHH_a3 = ParserExp.StructExp.Jnt.GHH.rz*deg2rad;
            Joint.GHH_a2 = ParserExp.StructExp.Jnt.GHH.ry*deg2rad;
            Joint.GHH_a1 = ParserExp.StructExp.Jnt.GHH.rx*deg2rad;
            Joint.GKH_a3 = ParserExp.StructExp.Jnt.GKH.rz*deg2rad;
            Joint.GKH_a2 = ParserExp.StructExp.Jnt.GKH.ry*deg2rad;
            Joint.GKH_a1 = ParserExp.StructExp.Jnt.GKH.rx*deg2rad;
            Joint.GBRK_a3 = ParserExp.StructExp.Jnt.GBRK.rx*deg2rad;
            Joint.GSBR_a3 = ParserExp.StructExp.Jnt.GSBR.rz*deg2rad;
            Joint.GSBL_a2 = ParserExp.StructExp.Jnt.GSBR.ry*deg2rad;
            Joint.GSBR_a1 = ParserExp.StructExp.Jnt.GSBR.rx*deg2rad;
            Joint.GSBL_a3 = ParserExp.StructExp.Jnt.GSBL.rz*deg2rad;
            Joint.GSBR_a2 = ParserExp.StructExp.Jnt.GSBL.ry*deg2rad;
            Joint.GSBL_a1 = ParserExp.StructExp.Jnt.GSBL.rx*deg2rad;
            Joint.GSR_a3 = ParserExp.StructExp.Jnt.GSR.rz*deg2rad;
            Joint.GSR_a2 = ParserExp.StructExp.Jnt.GSR.ry*deg2rad;
            Joint.GSR_a1 = ParserExp.StructExp.Jnt.GSR.rx*deg2rad;
            Joint.GSL_a3 = ParserExp.StructExp.Jnt.GSL.rz*deg2rad;
            Joint.GSL_a2 = ParserExp.StructExp.Jnt.GSL.ry*deg2rad;
            Joint.GSL_a1 = ParserExp.StructExp.Jnt.GSL.rx*deg2rad;
            Joint.GELR_a2 = ParserExp.StructExp.Jnt.GELR.ry*deg2rad;
            Joint.GELR_a1 = ParserExp.StructExp.Jnt.GELR.rx*deg2rad;
            Joint.GELL_a2 = ParserExp.StructExp.Jnt.GLK.ry*deg2rad;
            Joint.GELL_a1 = ParserExp.StructExp.Jnt.GLK.rx*deg2rad;
            Joint.GHAR_a3 = ParserExp.StructExp.Jnt.GHAR.rz*deg2rad;
            Joint.GHAR_a2 = ParserExp.StructExp.Jnt.GHAR.ry*deg2rad;
            Joint.GHAL_a3 = ParserExp.StructExp.Jnt.GHAL.rz*deg2rad;
            Joint.GHAL_a2 = ParserExp.StructExp.Jnt.GHAL.ry*deg2rad;

            Glob_REXP_BEC = rot123s(Joint.GHZ_a1,Joint.GHZ_a2,Joint.GHZ_a3);
            BEC_REXP_ULW = rot321s(Joint.GLK_a3,Joint.GLK_a2,Joint.GLK_a1);
            ULW_REXP_OLW = rot321s(Joint.GLL_a3,Joint.GLL_a2,Joint.GLL_a1);
            OLW_REXP_UBW = rot321s(Joint.GBL_a3,Joint.GBL_a2,Joint.GBL_a1);
            UBW_REXP_OBW = rot321s(Joint.GBB_a3,Joint.GBB_a2,Joint.GBB_a1);
            OBW_REXP_UHW = rot321s(Joint.GHB_a3,Joint.GHB_a2,Joint.GHB_a1);
            UHW_REXP_OHW = rot321s(Joint.GHH_a3,Joint.GHH_a2,Joint.GHH_a1);
            OHW_REXP_KO  = rot321s(Joint.GKH_a3,Joint.GKH_a2,Joint.GKH_a1);
            UHW_REXP_BRK = rot_z(Joint.GBRK_a3);
            BRK_REXP_SBR = rot321s(Joint.GSBR_a3,Joint.GSBR_a2,Joint.GSBR_a1);
            BRK_REXP_SBL = rot321s(Joint.GSBL_a3,Joint.GSBL_a2,Joint.GSBL_a1);
            SBR_REXP_OAR = rot321s(Joint.GSR_a3,Joint.GSR_a2,Joint.GSR_a1);
            SBL_REXP_OAL = rot321s(Joint.GSL_a3,Joint.GSL_a2,Joint.GSL_a1);
            OAR_REXP_UAR = rot21(Joint.GELR_a2,Joint.GELR_a1);
            OAL_REXP_UAL = rot21(Joint.GELL_a2,Joint.GELL_a1);
            UAR_REXP_HAR = rot32(Joint.GHAR_a3,Joint.GHAR_a2);
            UAL_REXP_HAL = rot32(Joint.GHAL_a3,Joint.GHAL_a2);
            Glob_REXP_BEC = rot_y(-pi/2)*rot_x(pi/2)*Glob_REXP_BEC;
            Glob_REXP_ULW = Glob_REXP_BEC * BEC_REXP_ULW;                                
            Glob_REXP_OLW = Glob_REXP_ULW * ULW_REXP_OLW;                                
            Glob_REXP_UBW = Glob_REXP_OLW * OLW_REXP_UBW;                                
            Glob_REXP_OBW = Glob_REXP_UBW * UBW_REXP_OBW;                                
            Glob_REXP_UHW = Glob_REXP_OBW * OBW_REXP_UHW;                                
            Glob_REXP_OHW = Glob_REXP_UHW * UHW_REXP_OHW;                                
            Glob_REXP_KO  = Glob_REXP_OHW * OHW_REXP_KO;                                 
            Glob_REXP_BRK = Glob_REXP_UHW *  rot_z(pi/2) * rot_x(pi) * UHW_REXP_BRK;     
            Glob_REXP_SBR = Glob_REXP_BRK * rot_y(-pi/2) * rot_x(-pi/2) * BRK_REXP_SBR;  
            Glob_REXP_SBL = Glob_REXP_BRK * rot_y(pi/2) * rot_x(pi/2) *   BRK_REXP_SBL;  
            Glob_REXP_OAR = Glob_REXP_SBR * SBR_REXP_OAR;                                
            Glob_REXP_OAL = Glob_REXP_SBL * SBL_REXP_OAL;                                
            Glob_REXP_UAR = Glob_REXP_OAR * rot_x(-pi/2) * OAR_REXP_UAR * rot_x(pi/2);   
            Glob_REXP_UAL = Glob_REXP_OAL * rot_x(-pi/2) * OAL_REXP_UAL * rot_x(pi/2);   
            Glob_REXP_HAR = Glob_REXP_UAR * UAR_REXP_HAR;                                
            Glob_REXP_HAL = Glob_REXP_UAL * UAL_REXP_HAL;
            
            % Global position of Joints from EXP
            GlobEXP_Pos_GLK = [ParserExp.StructExp.Jnt.GLK.x;ParserExp.StructExp.Jnt.GLK.y;ParserExp.StructExp.Jnt.GLK.z;]*0.001;% mm to m
            GlobEXP_Pos_GLL = [ParserExp.StructExp.Jnt.GLL.x;ParserExp.StructExp.Jnt.GLL.y;ParserExp.StructExp.Jnt.GLL.z;]*0.001;
            GlobEXP_Pos_GBL = [ParserExp.StructExp.Jnt.GBL.x;ParserExp.StructExp.Jnt.GBL.y;ParserExp.StructExp.Jnt.GBL.z;]*0.001;
            GlobEXP_Pos_GBB = [ParserExp.StructExp.Jnt.GBB.x;ParserExp.StructExp.Jnt.GBB.y;ParserExp.StructExp.Jnt.GBB.z;]*0.001;
            GlobEXP_Pos_GHB = [ParserExp.StructExp.Jnt.GHB.x;ParserExp.StructExp.Jnt.GHB.y;ParserExp.StructExp.Jnt.GHB.z;]*0.001;
            GlobEXP_Pos_GHH = [ParserExp.StructExp.Jnt.GHH.x;ParserExp.StructExp.Jnt.GHH.y;ParserExp.StructExp.Jnt.GHH.z;]*0.001;
            GlobEXP_Pos_GKH = [ParserExp.StructExp.Jnt.GKH.x;ParserExp.StructExp.Jnt.GKH.y;ParserExp.StructExp.Jnt.GKH.z;]*0.001;
            GlobEXP_Pos_PKSP = [ParserExp.StructExp.Sklt.PKSP.x;ParserExp.StructExp.Sklt.PKSP.y;ParserExp.StructExp.Sklt.PKSP.z;]*0.001;
            GlobEXP_Pos_GBRK = [ParserExp.StructExp.Jnt.GBRK.x;ParserExp.StructExp.Jnt.GBRK.y;ParserExp.StructExp.Jnt.GBRK.z;]*0.001;
            GlobEXP_Pos_GSBR = [ParserExp.StructExp.Jnt.GSBR.x;ParserExp.StructExp.Jnt.GSBR.y;ParserExp.StructExp.Jnt.GSBR.z;]*0.001;
            GlobEXP_Pos_GSBL = [ParserExp.StructExp.Jnt.GSBL.x;ParserExp.StructExp.Jnt.GSBL.y;ParserExp.StructExp.Jnt.GSBL.z;]*0.001;
            GlobEXP_Pos_GSR = [ParserExp.StructExp.Jnt.GSR.x;ParserExp.StructExp.Jnt.GSR.y;ParserExp.StructExp.Jnt.GSR.z;]*0.001;
            GlobEXP_Pos_GELR = [ParserExp.StructExp.Jnt.GELR.x;ParserExp.StructExp.Jnt.GELR.y;ParserExp.StructExp.Jnt.GELR.z;]*0.001;
            GlobEXP_Pos_GHAR = [ParserExp.StructExp.Jnt.GHAR.x;ParserExp.StructExp.Jnt.GHAR.y;ParserExp.StructExp.Jnt.GHAR.z;]*0.001;
            GlobEXP_Pos_PZSR = [ParserExp.StructExp.Sklt.PZSR.x;ParserExp.StructExp.Sklt.PZSR.y;ParserExp.StructExp.Sklt.PZSR.z;]*0.001;
            GlobEXP_Pos_GSL = [ParserExp.StructExp.Jnt.GSL.x;ParserExp.StructExp.Jnt.GSL.y;ParserExp.StructExp.Jnt.GSL.z;]*0.001;
            GlobEXP_Pos_GELL = [ParserExp.StructExp.Jnt.GELL.x;ParserExp.StructExp.Jnt.GELL.y;ParserExp.StructExp.Jnt.GELL.z;]*0.001;
            GlobEXP_Pos_GHAL = [ParserExp.StructExp.Jnt.GHAL.x;ParserExp.StructExp.Jnt.GHAL.y;ParserExp.StructExp.Jnt.GHAL.z;]*0.001;
            GlobEXP_Pos_PZSL = [ParserExp.StructExp.Sklt.PZSL.x;ParserExp.StructExp.Sklt.PZSL.y;ParserExp.StructExp.Sklt.PZSL.z;]*0.001;
            
            % Local position of Joints
            ULW_Pos_GLL = Glob_REXP_ULW * (GlobEXP_Pos_GLL - GlobEXP_Pos_GLK);
            OLW_Pos_GBL = Glob_REXP_OLW * (GlobEXP_Pos_GBL - GlobEXP_Pos_GLL);
            UBW_Pos_GBB = Glob_REXP_UBW * (GlobEXP_Pos_GBB - GlobEXP_Pos_GBL);
            OBW_Pos_GHB = Glob_REXP_OBW * (GlobEXP_Pos_GHB - GlobEXP_Pos_GBB);
            UHW_Pos_GHH = Glob_REXP_UHW * (GlobEXP_Pos_GHH - GlobEXP_Pos_GHB);
            OHW_Pos_GKH = Glob_REXP_OHW * (GlobEXP_Pos_GKH - GlobEXP_Pos_GHH);
            KO_Pos_PKSP = Glob_REXP_KO * (GlobEXP_Pos_PKSP - GlobEXP_Pos_GKH);
            UHW_Pos_GBRK = Glob_REXP_UHW * (GlobEXP_Pos_GBRK - GlobEXP_Pos_GHB);
            BRK_Pos_GSBR = Glob_REXP_BRK * (GlobEXP_Pos_GSBR - GlobEXP_Pos_GBRK);
            BRK_Pos_GSBL = Glob_REXP_BRK * (GlobEXP_Pos_GSBL - GlobEXP_Pos_GBRK);
            SBR_Pos_GSR = Glob_REXP_SBR * (GlobEXP_Pos_GSR - GlobEXP_Pos_GSBR);
            OAR_Pos_GELR = Glob_REXP_OAR * (GlobEXP_Pos_GELR - GlobEXP_Pos_GSR);
            UAR_Pos_GHAR = Glob_REXP_UAR * (GlobEXP_Pos_GHAR - GlobEXP_Pos_GELR);
            HAR_Pos_PZSR = Glob_REXP_HAR * (GlobEXP_Pos_PZSR - GlobEXP_Pos_GHAR);
            SBL_Pos_GSL = Glob_REXP_SBL * (GlobEXP_Pos_GSL - GlobEXP_Pos_GSBL);
            OAL_Pos_GELL = Glob_REXP_OAL * (GlobEXP_Pos_GELL - GlobEXP_Pos_GSL);
            UAL_Pos_GHAL = Glob_REXP_UAL * (GlobEXP_Pos_GHAL - GlobEXP_Pos_GELL);
            HAL_Pos_PZSL = Glob_REXP_HAL * (GlobEXP_Pos_PZSL - GlobEXP_Pos_GHAL);

            % set the angle posture for the upper body in CP operation
            Angle.GLK_a3 = -14.69*deg2rad; Angle.GLK_a2= 0*deg2rad;  Angle.GLK_a1 = 0*deg2rad;
            Angle.GLL_a3 = 16.60*deg2rad; Angle.GLL_a2 = 0*deg2rad; Angle.GLL_a1 = 0*deg2rad;
            Angle.GBL_a3 = 7.90*deg2rad; Angle.GBL_a2 = 0*deg2rad; Angle.GBL_a1 = 0*deg2rad;
            Angle.GBB_a3 = 9.60*deg2rad; Angle.GBB_a2 = 0*deg2rad; Angle.GBB_a1 = 0*deg2rad; 
            Angle.GHB_a3 = 12.50*deg2rad; Angle.GHB_a2 = 0*deg2rad; Angle.GHB_a1= 0*deg2rad;
            Angle.GHH_a3 = -5.70*deg2rad; Angle.GHH_a2 = 0*deg2rad; Angle.GHH_a1= 0*deg2rad;
            Angle.GKH_a3 = -1.40*deg2rad; Angle.GKH_a2 = 0*deg2rad; Angle.GKH_a1= 0*deg2rad;
            Angle.GBRK_a3 = 0*deg2rad;
            Angle.GSBR_a3 = 5.70*deg2rad; Angle.GSBR_a2 = 7.80*deg2rad; Angle.GSBR_a1 = 0*deg2rad;
            Angle.GSBL_a3 = 5.70*deg2rad; Angle.GSBL_a2 = 7.80*deg2rad; Angle.GSBL_a1 = 0*deg2rad;
            Angle.GSR_a3  = 75.60*deg2rad;  Angle.GSR_a2 = 33.0*deg2rad; Angle.GSR_a1  = 66.90*deg2rad;
            Angle.GSL_a3  = 75.60*deg2rad;  Angle.GSL_a2 = -33.0*deg2rad; Angle.GSL_a1  = -66.90*deg2rad;
            Angle.GELR_a2 = -54.0*deg2rad; Angle.GELR_a1 = -4.60*deg2rad;
            Angle.GELL_a2 = -54.0*deg2rad; Angle.GELL_a1 = 4.60*deg2rad;
            Angle.GHAR_a3 = 8.50*deg2rad; Angle.GHAR_a2 = -6.60*deg2rad;
            Angle.GHAL_a3 = 8.50*deg2rad; Angle.GHAL_a2 = 6.60*deg2rad;
            % rotation matrix between segments
            BEC_R_ULW = rot321s(Angle.GLK_a3,Angle.GLK_a2,Angle.GLK_a1);
            ULW_R_OLW = rot321s(Angle.GLL_a3,Angle.GLL_a2,Angle.GLL_a1);
            OLW_R_UBW = rot321s(Angle.GBL_a3,Angle.GBL_a2,Angle.GBL_a1);
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
            % Global rotation matrix
            Glob_R_BEC = rot_y(-pi/2)*rot_x(pi/2)*Glob_R_BEC;
            Glob_R_ULW = Glob_R_BEC * BEC_R_ULW;                                
            Glob_R_OLW = Glob_R_ULW * ULW_R_OLW;                                
            Glob_R_UBW = Glob_R_OLW * OLW_R_UBW;                                
            Glob_R_OBW = Glob_R_UBW * UBW_R_OBW;                                
            Glob_R_UHW = Glob_R_OBW * OBW_R_UHW;                                
            Glob_R_OHW = Glob_R_UHW * UHW_R_OHW;                                
            Glob_R_KO  = Glob_R_OHW * OHW_R_KO;                                 
            Glob_R_BRK = Glob_R_UHW *  rot_z(pi/2) * rot_x(pi) * UHW_R_BRK;     
            Glob_R_SBR = Glob_R_BRK * rot_y(-pi/2) * rot_x(-pi/2) * BRK_R_SBR;  
            Glob_R_SBL = Glob_R_BRK * rot_y(pi/2) * rot_x(pi/2) *   BRK_R_SBL;  
            Glob_R_OAR = Glob_R_SBR * SBR_R_OAR;                                
            Glob_R_OAL = Glob_R_SBL * SBL_R_OAL;                                
            Glob_R_UAR = Glob_R_OAR * rot_x(-pi/2) * OAR_R_UAR * rot_x(pi/2);   
            Glob_R_UAL = Glob_R_OAL * rot_x(-pi/2) * OAL_R_UAL * rot_x(pi/2);   
            Glob_R_HAR = Glob_R_UAR * UAR_R_HAR;                                
            Glob_R_HAL = Glob_R_UAL * UAL_R_HAL;
            % get Glob_Pos of joints
            Glob_Pos_GLL = Glob_Pos_GLK + Glob_R_ULW * ULW_Pos_GLL;
            Glob_Pos_GBL = Glob_Pos_GLL + Glob_R_OLW * OLW_Pos_GBL;
            Glob_Pos_GBB = Glob_Pos_GBL + Glob_R_UBW * UBW_Pos_GBB;
            Glob_Pos_GHB = Glob_Pos_GBB + Glob_R_OBW * OBW_Pos_GHB;
            Glob_Pos_GHH = Glob_Pos_GHB + Glob_R_UHW * UHW_Pos_GHH;
            Glob_Pos_GKH = Glob_Pos_GHH + Glob_R_OHW * OHW_Pos_GKH;
            Glob_Pos_PKSP = Glob_Pos_GKH + Glob_R_KO *  KO_Pos_PKSP;
            Glob_Pos_GBRK = Glob_Pos_GHB + Glob_R_UHW * UHW_Pos_GBRK;
            Glob_Pos_GSBR = Glob_Pos_GBRK + Glob_R_BRK * BRK_Pos_GSBR;
            Glob_Pos_GSBL = Glob_Pos_GBRK + Glob_R_BRK * BRK_Pos_GSBL;
            Glob_Pos_GSR = Glob_Pos_GSBR + Glob_R_SBR * SBR_Pos_GSR;
            Glob_Pos_GELR = Glob_Pos_GSR + Glob_R_OAR * OAR_Pos_GELR;
            Glob_Pos_GHAR = Glob_Pos_GELR + Glob_R_UAR * UAR_Pos_GHAR;
            Glob_Pos_PZSR = Glob_Pos_GHAR + Glob_R_HAR * HAR_Pos_PZSR;
            Glob_Pos_GSL = Glob_Pos_GSBL + Glob_R_SBL * SBL_Pos_GSL;
            Glob_Pos_GELL = Glob_Pos_GSL + Glob_R_OAL * OAL_Pos_GELL;
            Glob_Pos_GHAL = Glob_Pos_GELL + Glob_R_UAL * UAL_Pos_GHAL;
            Glob_Pos_PZSL = Glob_Pos_GHAL + Glob_R_HAL * HAL_Pos_PZSL;
            
            % CoM of the segments in CP posture
            Glob_Pos_VecCoM{1} = Glob_R_BEC*R.Subject.Segments(SegIndex).CoM.LocCoord;
            % Hacerlos con la superficie
            Glob_Pos_VecCoM{2} = Glob_Pos_GBB + ((Glob_Pos_GLL - Glob_Pos_GBB)/norm(Glob_Pos_GLL - Glob_Pos_GBB))* MidTrunkCoM;
            Glob_Pos_VecCoM{3} = Glob_Pos_GHH + ((Glob_Pos_GBB - Glob_Pos_GHH)/norm(Glob_Pos_GBB - Glob_Pos_GHH))* UpTrunkCoM;
            Glob_Pos_VecCoM{4} = Glob_Pos_PKSP + ((Glob_Pos_GHH - Glob_Pos_PKSP)/norm(Glob_Pos_GHH - Glob_Pos_PKSP))* HeadNeckCoM;
            %
            Glob_Pos_VecCoM{5} = Glob_Pos_PZSR + ((Glob_Pos_GHAR - Glob_Pos_PZSR)/norm(Glob_Pos_GHAR - Glob_Pos_PZSR))* HandCoM;
            Glob_Pos_VecCoM{6} = Glob_Pos_PZSL + ((Glob_Pos_GHAL - Glob_Pos_PZSL)/norm(Glob_Pos_GHAL - Glob_Pos_PZSL))* HandCoM;
            Glob_Pos_VecCoM{7} = Glob_Pos_GHAR + ((Glob_Pos_GELR - Glob_Pos_GHAR)/norm(Glob_Pos_GELR - Glob_Pos_GHAR))* ForearmCoM;
            Glob_Pos_VecCoM{8} = Glob_Pos_GHAL + ((Glob_Pos_GELL - Glob_Pos_GHAL)/norm(Glob_Pos_GELL - Glob_Pos_GHAL))* ForearmCoM;
            Glob_Pos_VecCoM{9} = Glob_Pos_GELR + ((Glob_Pos_GSR - Glob_Pos_GELR)/norm(Glob_Pos_GSR - Glob_Pos_GELR))* UpperarmCoM;
            Glob_Pos_VecCoM{10} = Glob_Pos_GELL + ((Glob_Pos_GSL - Glob_Pos_GELL)/norm(Glob_Pos_GSL - Glob_Pos_GELL))* UpperarmCoM;
            VecMass{1} = R.Subject.Segments(SegIndex).Mass;
            VecMass{2} = MidTrunkMass;
            VecMass{3} = UpTrunkMass;
            VecMass{4} = HeadNeckMass;
            VecMass{5} = HandMass;
            VecMass{6} = HandMass;
            VecMass{7} = ForearmMass;
            VecMass{8} = ForearmMass;
            VecMass{9} = UpperarmMass;
            VecMass{10} = UpperarmMass;
            % Calc of the CoM of the Upperbody
            SumVecCoMMass = 0;
            SumMass = 0;
            for i=1:size(Glob_Pos_VecCoM,1)
                SumVecCoMMass = SumVecCoMMass + Glob_Pos_VecCoM{i} * VecMass{i};
                SumMass = SumMass + VecMass{i};
            end
            Glob_Pos_CoM = SumVecCoMMass/SumMass;
            GHZIndex = getVecIndex('GHZ',R.Subject.Segments(SegIndex).LocalPoints);
            GHZPosInq = R.Subject.Segments(SegIndex).LocalPoints(GHZIndex).Point.PosInq;
            Glob_Pos_GHZ = R.q_t(1,GHZPosInq:2+GHZPosInq)';
            R.Subject.Segments(SegIndex).CoM.LocCoord = Glob_R_BEC * (Glob_Pos_CoM - Glob_Pos_GHZ);
            R.Subject.Segments(SegIndex).Mass = SumMass;
            % Calc of the I of the Upperbody
            % For each axis Iaxis = Icm + mass*r2
            Anat_R_Glob = [-1,0,0;0,-1,0;0,0,1];
            Glob_VecI{1} =[];
            Glob_VecI{2} = (Anat_R_Glob'*[MidTrunkIAntPos,0,0;0,MidTrunkIMedLat,0;0,0,MidTrunkILong]*Anat_R_Glob)*0.0001;  % to kg*m2
            Glob_VecI{3} = (Anat_R_Glob'*[UpTrunkIAntPos,0,0;0,UpTrunkIMedLat,0;0,0,UpTrunkILong]*Anat_R_Glob)*0.0001;
            Glob_VecI{4} = (Anat_R_Glob'*[HeadNeckIAntPos,0,0;0,HeadNeckIMedLat,0;0,0,HeadNeckILong]*Anat_R_Glob)*0.0001;
            Glob_VecI{5} = (Anat_R_Glob'*[HandIAntPos,0,0;0,HandIMedLat,0;0,0,HandILong]*Anat_R_Glob)*0.0001;
            Glob_VecI{6} = (Anat_R_Glob'*[HandIAntPos,0,0;0,HandIMedLat,0;0,0,HandILong]*Anat_R_Glob)*0.0001;
            Glob_VecI{7} = (Anat_R_Glob'*[ForearmIAntPos,0,0;0,ForearmIMedLat,0;0,0,ForearmILong]*Anat_R_Glob)*0.0001;
            Glob_VecI{8} = (Anat_R_Glob'*[ForearmIAntPos,0,0;0,ForearmIMedLat,0;0,0,ForearmILong]*Anat_R_Glob)*0.0001;
            Glob_VecI{9} = (Anat_R_Glob'*[UpperarmIAntPos,0,0;0,UpperarmIMedLat,0;0,0,UpperarmILong]*Anat_R_Glob)*0.0001;
            Glob_VecI{10} = (Anat_R_Glob'*[UpperarmIAntPos,0,0;0,UpperarmIMedLat,0;0,0,UpperarmILong]*Anat_R_Glob)*0.0001;
            Sum_ICoM = zeros(3,3);
            for i=2:length(Glob_VecI)
                SteinerMat(1,1)= (Glob_Pos_VecCoM{i}(2) - Glob_Pos_CoM(2))^2 + (Glob_Pos_VecCoM{i}(3) - Glob_Pos_CoM(3))^2;
                SteinerMat(2,2)= (Glob_Pos_VecCoM{i}(1) - Glob_Pos_CoM(1))^2 + (Glob_Pos_VecCoM{i}(3) - Glob_Pos_CoM(3))^2;
                SteinerMat(3,3)= (Glob_Pos_VecCoM{i}(1) - Glob_Pos_CoM(1))^2 + (Glob_Pos_VecCoM{i}(2) - Glob_Pos_CoM(2))^2;
                SteinerMat(1,2)= -(Glob_Pos_VecCoM{i}(1) - Glob_Pos_CoM(1))*(Glob_Pos_VecCoM{i}(2) - Glob_Pos_CoM(2));
                SteinerMat(1,3)= -(Glob_Pos_VecCoM{i}(1) - Glob_Pos_CoM(1))*(Glob_Pos_VecCoM{i}(3) - Glob_Pos_CoM(3));
                SteinerMat(2,3)= -(Glob_Pos_VecCoM{i}(2) - Glob_Pos_CoM(2))*(Glob_Pos_VecCoM{i}(3) - Glob_Pos_CoM(3));
                SteinerMat(2,1)= SteinerMat(1,2);
                SteinerMat(3,1)= SteinerMat(1,3);
                SteinerMat(3,2)= SteinerMat(2,3);
                Sum_ICoM = Sum_ICoM + Glob_VecI{i}+ VecMass{i}*SteinerMat;
            end
            BEC_ICoM = Glob_R_BEC'*Sum_ICoM*Glob_R_BEC;
            BEC_IGHZ = R.Subject.Segments(SegIndex).I;
            SteinerMat(1,1)= (Glob_Pos_GHZ(2) - Glob_Pos_CoM(2))^2 + (Glob_Pos_GHZ(3) - Glob_Pos_CoM(3))^2;
            SteinerMat(2,2)= (Glob_Pos_GHZ(1) - Glob_Pos_CoM(1))^2 + (Glob_Pos_GHZ(3) - Glob_Pos_CoM(3))^2;
            SteinerMat(3,3)= (Glob_Pos_GHZ(1) - Glob_Pos_CoM(1))^2 + (Glob_Pos_GHZ(2) - Glob_Pos_CoM(2))^2;
            SteinerMat(1,2)= -(Glob_Pos_GHZ(1) - Glob_Pos_CoM(1))*(Glob_Pos_GHZ(2) - Glob_Pos_CoM(2));
            SteinerMat(1,3)= -(Glob_Pos_GHZ(1) - Glob_Pos_CoM(1))*(Glob_Pos_GHZ(3) - Glob_Pos_CoM(3));
            SteinerMat(2,3)= -(Glob_Pos_GHZ(2) - Glob_Pos_CoM(2))*(Glob_Pos_GHZ(3) - Glob_Pos_CoM(3));
            SteinerMat(2,1)= SteinerMat(1,2);
            SteinerMat(3,1)= SteinerMat(1,3);
            SteinerMat(3,2)= SteinerMat(2,3);
            R.Subject.Segments(SegIndex).I = (BEC_ICoM + BEC_IGHZ + VecMass{1}*SteinerMat);
            
        end
        function MarkerCoordsInterp = checkMarkers(R,MarkerNames, MarkerCoords) 
            % CHECKMARKER check if there are missing markers , writes a report and interpolate with few frames
            %
            %   checkMarker(MarkerNames, MarkerCoords)
            %   Inputs:
            %     + MarkerNames is a cell (nMarker x 1) with the names of the skin-markers. Each
            %       element of the cell is a string.
            %     + MarkerCoords is a double array(nFrames x 3*nMarker) with the 3D coordinates of
            %       skin-markers. The coordinates x,y,z of skin-marker i are in position
            %       MarkerCoords(3*i-2,3*i-1,3*i)
            %   Outputs:
            %     +  MarkerCoordsInterp is a double array(nFrames x 3*nSmark) with the interpolated 3D 
            %       coordinates of skin-markers. The coordinates x,y,z of skin-marker i are in position
          
            % ----------------------------------------------------------------------------
            % Check the missing skin-markers and interpolate
            % ----------------------------------------------------------------------------

            
            %  Display the Interpolation information
            if strcmpi(R.Settings.Interpolation.Method,'linear')
                str = ['  Interpolating gaps in marker trajectories... Method: linear, GapInterpolated: ',num2str(R.Settings.Interpolation.NInterFrames),' frames'];
                if(R.Settings.Display == 1 || R.Settings.Display == 2), disp(str); end
                fprintf(R.ExperLogFileId, '%s\n', str);                
            elseif isempty(R.Settings.Interpolation.Method)
                str = '  No interpolation applied to gaps in marker trajectories';
                if(R.Settings.Display == 1 || R.Settings.Display == 2), disp(str); end
                fprintf(R.ExperLogFileId, '%s\n', str);                
            else
                str = ['  Interpolation method "',R.Settings.Interpolation.Method,'" is unknown'];
                fprintf(R.ExperLogFileId, '%s\n', str);
                error(str);                
            end
            
            
            % sizes
            [nFrames,nCoords] = size(MarkerCoords);
            MarkerCoordsInterp = MarkerCoords;
            
            % initialize variables
            occlfr = [];    % The missing frames numbers will be stored here temporarily
            History = {};   % The msssing frames numbers summary will be stored in this cell
            % here for each frame so that then they are monitored in order
            % This cell size is (number of missing markers x number of missing
            % intervals). The missing intervals could be a 1x2 vector that
            % contains the first and the last missing frames of the interval
            % or a 1x1 scalar; the missing frame
            k = 1;          % This counter is for the missing frames storage
            q = 1;          % This counter is for the rows of the History cell
            
            for j = 1:3:nCoords   % The matrix will be scanned by columns of three (x y z of each marker)
                p = 1;      % This counter is for determine whether it is the first time an NaN is
                % encountered for each marker or not
                r = 3;      % This counter is for the variable columns of the History cell
                
                % Get the number of the marker corresponding to the scanned column
                MarkerName = MarkerNames{(j-1)/3+1};
                
                for i = 1:nFrames % i is the counter for the rows of the column being scanned
                    if isnan( MarkerCoords(i,j) ) && i ~= nFrames   % If a NaN is found it is missing skin-marker frame
                        occlfr(k) = i;                             % The missing frame number is stored
                        k = k + 1;
                    elseif (~isnan( MarkerCoords(i,j) ) && norm(occlfr) ~= 0) || (isnan( MarkerCoords(i,j) ) && i == nFrames)
                        % If an non-NaN value is encountered and there have been missing frames
                        % then the History cell must be loaded
                        if isnan( MarkerCoords(i,j) ) && i == nFrames
                            occlfr(k) = i;
                            k = k + 1;
                        end
                        
                        if k == 2           % If there is only one missing frame then the History cell will be stored here
                            if p == 1       % If it is the first encountered NaN for this marker save the number of the
                                % marker and the missing frame
                                History(q,1:2) = {MarkerName, occlfr(1)};
                                q = q+1;
                            else            % If it isn't the first encountered missing frame
                                if q == 1   % If it is the first row of History
                                    History{q,r} = occlfr(1);
                                    r = r+1;
                                else        % If it isn't the first row of History
                                    History{q-1,r} = occlfr(1); % Add the missing frame to the marker line of History
                                    % (notice that when the line was created the q value was
                                    % increased so now to refer to the marker
                                    % line it must be referred to q-1 not q)
                                    r = r+1;
                                end
                            end
                            p = p+1;
                            if R.Settings.Interpolation.NInterFrames > 0
                                if occlfr(1) == 1
                                    MarkerCoordsInterp(i-1,j:j+2) = MarkerCoordsInterp(i,j:j+2);
                                elseif occlfr(end) == nFrames
                                    MarkerCoordsInterp(i,j:j+2) = MarkerCoordsInterp(i-1,j:j+2);
                                else
                                    MarkerCoordsInterp(i-1,j:j+2) = (MarkerCoordsInterp(i-2,j:j+2) + MarkerCoordsInterp(i,j:j+2))/2;
                                    %                                 SPLINE interpolation
                                    %                                 time = 1:nFrames;
                                    %                                 Frames = [1:occlfr(1)-1,occlfr(end)+1:nFrames];
                                    %                                 MarkerCoordsValue = [MarkerCoordsInterp(1:occlfr(1)-1,j:j+2);MarkerCoordsInterp(occlfr(end)+1:nFrames,j:j+2)];
                                    %                                 pp_pos = spline(Frames,MarkerCoordsValue');
                                    %                                 OclusedMarker = ppval(pp_pos,time)';
                                    % %                                 plot(Frames,MarkerCoordsValue(:,1)); hold on,
                                    % %                                 plot(time,MarkerCoordsInterpEste(:,1),'r');
                                    %                                 MarkerCoordsInterp(:,j:j+2)= OclusedMarker;
                                end
                            end
%                             str = (['The marker ',MarkerName, ' is missing and is interpolate in frame',num2str(i-1)]);
%                             fprintf(R.ExperLogFileId,'%s\n', str);
%                             disp(str);
                            occlfr = [];
                            k = 1;
                            
                        % The algorithm is the same as above but in this case the first and the last missing
                        % frames are stored in History
                        else
                            if p == 1
                                History(q,1:2) = {MarkerName, [occlfr(1) occlfr(end)]};
                                q = q+1;
                            else
                                if q == 1
                                    History{q,r} = [occlfr(1) occlfr(end)];
                                    r = r+1;
                                else
                                    History{q-1,r} = [occlfr(1) occlfr(end)];
                                    r = r+1;
                                end
                            end
                            p = p+1;
                            % ---------------------------------------------------
                            % Linear Interpolation for Nan < NInterPolated Frames
                            % ---------------------------------------------------                            
                            nOcclFr = length(occlfr);
                            
                            if nOcclFr < R.Settings.Interpolation.NInterFrames
                                if occlfr(1) == 1
                                    for inter = i-nOcclFr:i-1
                                        MarkerCoordsInterp(inter,j:j+2) = MarkerCoordsInterp(i,j:j+2);
                                    end
                                elseif occlfr(end) == nFrames
                                    for inter = i-(nOcclFr-1):i
                                        MarkerCoordsInterp(inter,j:j+2) = MarkerCoordsInterp(i-nOcclFr,j:j+2);
                                    end
                                else
                                    %                                         a)SPLINE interpolation
                                    %                                         time = 1:nFrames;
                                    %                                         Frames = [1:occlfr(1)-1,occlfr(end)+1:nFrames];
                                    %                                         MarkerCoordsValue = [MarkerCoordsInterp(1:occlfr(1)-1,j:j+2);MarkerCoordsInterp(occlfr(end)+1:nFrames,j:j+2)];
                                    %                                         pp_pos = spline(Frames,MarkerCoordsValue');
                                    %                                         OclusedMarker = ppval(pp_pos,time)';
                                    %                                         %                                     plot(Frames,MarkerCoordsValue(:,1)); hold on,
                                    %                                         %                                     plot(time,MarkerCoordsInterpEste(:,1),'r');
                                    %                                         MarkerCoordsInterp(:,j:j+2)= OclusedMarker;
                                    % b) Linear interpolation
                                    for inter = i-nOcclFr:i-1
                                        MarkerCoordsInterp(inter,j:j+2) = MarkerCoordsInterp(inter-1,j:j+2) +  ...
                                            (MarkerCoordsInterp(i,j:j+2) - MarkerCoordsInterp(i-nOcclFr-1,j:j+2))/(nOcclFr+1);
                                    end
                                end
                            end
                            
                            % ---------------------------------
                            % END of Linear Interpolation
                            % ---------------------------------
                            occlfr = [];
                            k = 1;
                        end
                    end
                end
            end            
            
            % ----------------------------------------------------------------------------
            % Print the missing markers information (Frame numbers with missing)
            % ----------------------------------------------------------------------------
            % Monitor the missing markers in order (only if there have been missing frames)
            [rows, columns] = size(History);
            if rows ~= 0
                for i = 1:rows
                    % display info
                    str = sprintf('    Marker "%s":', History{i,1});
                    fprintf(R.ExperLogFileId,'%s\n', str);
                    if(R.Settings.Display == 1 || R.Settings.Display == 2)
                        disp(str);
                    end
                    % print the frame numbers that marker i is missing
                    for j = 2:columns
                        if length(History{i,j}) == 1
                            if  R.Settings.Interpolation.NInterFrames > 0
                                str = sprintf('      Frame %d is missing. Gap filled by linear interpolation.', History{i,j});
                            else
                                str = sprintf('      Frame %d is missing. Gap not filled.', History{i,j});
                            end
                            if(R.Settings.Display == 1  || R.Settings.Display == 2)
                                disp(str);
                            end
                            fprintf(R.ExperLogFileId,'%s\n', str);
                        elseif length(History{i,j}) == 2
                            str1 = sprintf('      Frames %d to %d are missing.', History{i,j});
                            if History{i,j}(2)-History{i,j}(1) <  R.Settings.Interpolation.NInterFrames
                                str2 = sprintf(' Gap filled by linear interpolation.');
                            else
                                str2 = sprintf(' Gap not filled.');
                            end
                            str = [str1,str2];
                            
                            fprintf(R.ExperLogFileId,'%s\n', str);
                            if(R.Settings.Display == 1  || R.Settings.Display == 2)
                                disp(str);
                            end
                        end
                    end
                end
                % -----------------------------------------------------------------------------------------
                % Summary of missing frames: percentage of missing frames for each skin-marker
                % -----------------------------------------------------------------------------------------
                str = sprintf(['\n    Summary of missing frames:\n', ...
                    '      Marker names     Missing frames/Total frames      %% of missing frames']);
                
                fprintf(R.ExperLogFileId,'%s\n', str);
                if(R.Settings.Display == 1 || R.Settings.Display == 2)
                    disp(str);
                end
                % Percentage of missing frames for each marker
                for i = 1:rows
                    nMissingFrames_Markeri = 0;
                    for j = 2:columns
                        if length(History{i,j}) == 1
                            nMissingFrames_Markeri = nMissingFrames_Markeri + 1;
                        elseif length(History{i,j}) == 2
                            nMissFrames = (History{i,j}(2) - History{i,j}(1)) + 1;
                            nMissingFrames_Markeri = nMissingFrames_Markeri + nMissFrames;
                        end
                    end
                    str = sprintf('%18s    %14d/%-12d             %5.1f%%',...
                        History{i,1},nMissingFrames_Markeri,nFrames,100*nMissingFrames_Markeri/nFrames);
                   
                    fprintf(R.ExperLogFileId,'%s\n', str);
                    if(R.Settings.Display == 1 || R.Settings.Display == 2)
                        disp(str);
                    end
                end
            else
                str = sprintf('   There are not any missing marker in the motion file');

                fprintf(R.ExperLogFileId,'%s\n', str);
                if(R.Settings.Display == 1 || R.Settings.Display == 2)
                    disp(str);
                end
            end
            
        end
        function doID(R,ResultsPath)
            if R.ErrorID == 0
                time = tic;
                str = sprintf(['\n  -------------------------------------------------------------------------\n', ...
                    '   Inverse Dynamic Motion Reconstruction (DHID) ', ...
                    '\n  -------------------------------------------------------------------------']);
                dispInfo(R.Settings.Display, str, R.ExperLogFileId);                
                % Inicializo la clase INVERSE_DYNAMIC
                Message = '   Solving inverse dynamic problem ...';
                dispInfo(R.Settings.Display, Message,R.ExperLogFileId);
                ID = INVERSE_DYNAMICS();
                if ~isempty(R.Subject.Segments(1).Distals)
                    for i=1:size(R.Subject.Segments,1)
                        R.Subject.Segments(i).w = [];
                        R.Subject.Segments(i).wdot = [];
                        R.Subject.Segments(i).vdotOrBCS = [];
                        R.Subject.Segments(i).Seg_F_CoM = [];
                        R.Subject.Segments(i).Seg_M_CoM = [];
                    end
                    for i=1:size(R.Subject.Joints,1)
                        R.Subject.Joints{i}.F = [];
                        R.Subject.Joints{i}.M = [];
                        R.Subject.Joints{i}.FGlob = [];
                        R.Subject.Joints{i}.MGlob = [];
                        R.Subject.Joints{i}.Angles = [];
                        R.Subject.Joints{i}.Anglesdot = [];
                        R.Subject.Joints{i}.Angles2dot = [];
                    end
                else
                    % Obtengo la jerarquía del modelo
                    R.Subject.setModelHierarchy();
                end
                R.Subject.calcAnglesData(R.q_t,R.Deltat,R.Settings);
                % Obtengo los segmentosTerminales e iniciales
                InitSegment = R.Subject.getInitSegment();
                EndSegments = R.Subject.getEndSegments();
                % Calculo las fuerzas en cada frame
                NFrames = size(R.q_t,1);
                if strcmp(R.Settings.ExType,'CP')
                    R.calcUpperBodyInertialParameters();
                end
                for i = 1 : NFrames
                    % Calculo las fuerzas y momentos de inercia de cada segmento partiendo del Ground hacia fuera
                    InitSegment.w = [0;0;0];
                    InitSegment.wdot = [0;0;0];
                    InitSegment.vdotOrBCS = [0;0;9.80665];
    %                 InitSegment.vdotOrBCS = [0;0;0];
                    InitSegment.calcInertiaForces(ID,R.q_t(i,:),R.qdot2_t(i,:),i);
                    % Calculo las fuezas y momentos en cada Joint empezando de los extremos
                    for j=1:size(EndSegments,1)
                        EndSegments(j).calcJointForces(ID,R.q_t(i,:),i);
                    end

                end
                R.ResultType = 'ID';
                R.resultsID(ResultsPath);
                [H, MI, S] = second2HMS(toc(time));
                str =   ['   ID Solver Time: ',num2str(H),' hour(s) ',num2str(MI),' minute(s) and ',num2str(S),' second(s).'];
                if(R.Settings.Display == 1 || R.Settings.Display == 2)
                    disp(str);
                end
                fprintf(R.ExperLogFileId,'%s\n', str);
            else
                str = ('  WARNING: Imposible to do ID because Imposible definition of Sensor LCS.');
                disp(str);
                fprintf(R.ExperLogFileId, '%s\n', str);
            end
            
        end
        function ExperVar = doIK(R,GuidedFile,ResultsPath,GraphicPath)
            [Frequency,MarkerNames,MarkerCoords,NFrames] = R.readMotionFile();
            if isempty(GuidedFile)
                GuidedMarkers = {'AllMarkers', []};
                GuidedCoords = [];
                R.Deltat = 1/Frequency;
                MessageLines = ['GuidedVarsFile not defined in Interface file.\n',...
                     '           All model markers found in motion file are used to guide the model.']; 
                dispWarning(MessageLines, 0, R.ExperLogFileId)
            else
                if nargout(GuidedFile) == 2     % marker trajectories defined in a C3D file
                    [GuidedMarkers,GuidedCoords] = feval(GuidedFile, NFrames);
                    R.Deltat = 1/Frequency;
                elseif nargout(GuidedFile) == 3 % marker trajectories defined in file GuideFile
                    [GuidedMarkers,GuidedCoords,Deltat] = feval(GuidedFile, NFrames);
                    R.Deltat =  Deltat;
                elseif nargout(GuidedFile) == 4 % marker trajectories defined in file GuideFile
                    [GuidedMarkers,GuidedCoords,Deltat,R.q0User] = feval(GuidedFile, NFrames);
                    if isnan(Frequency)
                        if isempty(Deltat) || Deltat <= 0
                            error('[The sampling frequency value (',num2str(Deltat),' Hz) define by user must be positive.')
                        else
                            R.Deltat =  Deltat;
                        end
                    else
                        if Deltat == 1/Frequency
                            R.Deltat = 1/Frequency;
                        else
                            error(['Sampling frequency of marker trajectories (',num2str(1/Deltat),' Hz) defined ' ...
                                'by user is different from sampling frequency (',num2str(Frequency),' Hz) in ' ...
                                'motion file. Please change sampling frequency in file ',GuidedFile,'.m']);
                        end
                    end
                end
            end
            
            % check if a C3D file has been defined when marker trajectories are not provided in GuidedFile (*_guidedVars.m) 
            if (size(GuidedMarkers,2) == 2) && isempty(MarkerNames)
                % if conditions:
                %   1st condition is true when markers trajectories are not provided in GuidedFile
                %   2nd condition is true when a C3D file is not defined                
                error('Markers trajectories not define in "motion definition file",\n   then at least one C3D file must be defined in the Main file')
            end
            
            % This functions only are executed once in the experiment
            ExpeVar = 0;
            if isempty(R.z)
                if ~strcmpi(R.Settings.SolverIK,'SODERKVIST')
                    if exist([R.AddCtrPath,'additionalCtr.m'],'file')
                        additionalCtr(R.Subject);
                        PhiVectors = R.Subject.mkAdditionalCtrs();
                        R.Subject.mkVectorPhiq(PhiVectors);
                    else
                        % DO NOTHING. This version is for Biomech lectures only and
                        % additionalCtr.m IS NOT USED
%                         MessageLines = ['Additional constraints will not be added to the model.\n',...
%                             '           Function additionalCtr.m not found in path: ',getPrintPath(R.AddCtrPath)];
%                         dispWarning(MessageLines, 0, R.ExperLogFileId)
                    end
                    R.wFileFillPhi();
                    R.wFileFillPhiqSparse();
                end
                R.addRecontructionMarkers(GuidedMarkers);
                R.fillExperimentValues(GuidedCoords);
%                 R.Subject.setModelHierarchy();
                ExpeVar = 1;
            end
            if ~strcmpi(R.Settings.SolverIK,'SODERKVIST')
                % This function only executed once for subject
                if isempty(R.Par)
                    % Add PosInz to Markers of the Subject
                    NMarkers = size(R.Markers,1);
                    if ExpeVar == 0
                        for i=1:NMarkers
                            R.Subject.Markers(i).PosInz = R.Markers(i).PosInz;
                        end
                    end
                    NPar = 0;
                    NSeg    = size(R.Subject.Segments,1);
                    for i=1:NSeg
                        if R.Subject.Segments(i).Fixed ~= 1 % If is NOT Ground
                        %if ~strcmpi(R.Subject.Segments(i).Name,'Ground')
                            NPoints = size(R.Subject.Segments(i).LocalPoints,1);
                            NMarkers = size(R.Subject.Segments(i).LocalMarkers,1);
                            NVectors = size(R.Subject.Segments(i).LocalVectors,1);
                            for j=1:NPoints
                                R.Par(3*j-2+NPar)= R.Subject.Segments(i).LocalPoints(j).LocCoord(1);
                                R.Par(3*j-1+NPar)= R.Subject.Segments(i).LocalPoints(j).LocCoord(2);
                                R.Par(3*j+NPar)  = R.Subject.Segments(i).LocalPoints(j).LocCoord(3);
                            end
                            for j=1:NMarkers
                                R.Par(3*j-2+3*NPoints+NPar)= R.Subject.Segments(i).LocalMarkers(j).LocCoord(1);
                                R.Par(3*j-1+3*NPoints+NPar)= R.Subject.Segments(i).LocalMarkers(j).LocCoord(2);
                                R.Par(3*j+3*NPoints+NPar)  = R.Subject.Segments(i).LocalMarkers(j).LocCoord(3);
                            end
                            for j=1:NVectors
                                R.Par(3*j-2+3*NPoints+3*NMarkers+NPar)= R.Subject.Segments(i).LocalVectors(j).LocCoord(1);
                                R.Par(3*j-1+3*NPoints+3*NMarkers+NPar)= R.Subject.Segments(i).LocalVectors(j).LocCoord(2);
                                R.Par(3*j+3*NPoints+3*NMarkers+NPar)  = R.Subject.Segments(i).LocalVectors(j).LocCoord(3);
                            end
                            NPar = 3*NPoints + 3*NMarkers + 3*NVectors+NPar;
                        end
                    end
                end
            end
            R.getMarkerGuidedValues(MarkerNames,MarkerCoords,GuidedCoords,GuidedMarkers);
            if ~strcmpi(R.Settings.SolverIK,'SODERKVIST')
                R.initq0(GuidedFile);
            end
            
            SettingsVec = [R.Settings.Results.Ramsis; R.Settings.Results.PAM;R.Settings.Results.InitPosture; ...
                R.Settings.Results.CompPlayback; R.Settings.Results.Position; R.Settings.Results.Velocity; ...
                R.Settings.Results.Acceleration; R.Settings.Results.MarkerError; R.Settings.Results.RawMarkerTrajectory; ...
                R.Settings.Results.AcondMarkerTrajectory;R.Settings.Results.Sensor; R.Settings.Results.JointEffort];
            if sum(SettingsVec) == 1  &&  R.Settings.Results.InitPosture == 1 % only initial approx. is asked for
                dispWarning({'Only results for the initial approximation playback are selected.';....
                    'Therefore, motion reconstruction will not be performed'}, R.Settings.Display, R.ExperLogFileId);
            else % IK calculations must be performed
                if ~strcmpi(R.Settings.SolverIK,'SODERKVIST')
                    R.solverIK();
                else
                    R.solverSoderkvist();
                end
                R.ResultType = 'IK';
                if R.Settings.Results.Velocity == 1 || R.Settings.Results.Acceleration == 1 || strcmpi(R.Settings.Type,'ID')
                    [R.qdot_t,R.qdot2_t] = calcVelAcc(R.q_t,R.Deltat);
                end
            end
            str = sprintf(['\n   --------------------------------------\n', ...
                           '    Results files ', ...
                           '\n   --------------------------------------']);
            if(R.Settings.Display == 1 || R.Settings.Display == 2)
                disp(str);
            end
            fprintf(R.ExperLogFileId, '%s\n', str);
            
            R.resultsIK(ResultsPath,GraphicPath);
            % Return to the EXPERIMENT
            ExperVar.z =  R.z;
            ExperVar.Markers = R.Markers;
            ExperVar.Wm = R.Wm;
            ExperVar.Ws = R.Ws;
            ExperVar.gs_t = R.gs_t;
            ExperVar.PoszInq = R.PoszInq;
            ExperVar.FilePhiName = R.FilePhiName;
            ExperVar.FilePhiqName = R.FilePhiqName;
            ExperVar.Par = R.Par;
        end
        function DataPressure = filtAndsyncroPData(R,RawDataPressure,MotFrames,MotFrequency)
            % La función filtra y sincroniza los valores del pressure map a los tiempos de la captura del movimiento
%             PressureTimeVar = 0;
            PressureTime = 0;
            MotionTime = (0:1/MotFrequency:MotFrames/MotFrequency)';
            DataPressure.Seat = [];
            DataPressure.Back = [];
            % Calcular los intervalos de tiempo del Pressure Map
            for i=1:size(RawDataPressure.Seat,2)-1
                SepIndex = strfind(RawDataPressure.Seat(i).Time,':');
                Min = str2num(RawDataPressure.Seat(i).Time(SepIndex(1)+1:SepIndex(2)-1));
                Sec = str2num(RawDataPressure.Seat(i).Time(SepIndex(2)+1:end-1));
                SecIni = 60*Min + Sec;
                SepIndex = strfind(RawDataPressure.Seat(i+1).Time,':');
                Min = str2num(RawDataPressure.Seat(i+1).Time(SepIndex(1)+1:SepIndex(2)-1));
                Sec = str2num(RawDataPressure.Seat(i+1).Time(SepIndex(2)+1:end-1));
                SecEnd= 60*Min + Sec;
                TimeVar_i = SecEnd - SecIni;
%                 PressureTimeVar = [PressureTimeVar;TimeVar_i];
                PressureTime = [PressureTime;TimeVar_i+PressureTime(i)];
            end
            % Obtener los valores filtrados y sincrónizados de cada casilla del presuremap
            for i=1:size(RawDataPressure.Seat(1).Values,1)
                for j=1:size(RawDataPressure.Seat(1).Values,2)
                    SeatData_k = [];
                    BackData_k = [];
                    for k=1:size(RawDataPressure.Seat,2)
                        SeatData_k = [SeatData_k,RawDataPressure.Seat(k).Values(i,j)];
                        BackData_k = [BackData_k,RawDataPressure.Back(k).Values(i,j)];
                    end
                    DataPressure.Seat = R.filtAndsyncroMotFramePData(SeatData_k,DataPressure.Seat,MotionTime,PressureTime,i,j);
                    DataPressure.Back = R.filtAndsyncroMotFramePData(BackData_k,DataPressure.Back,MotionTime,PressureTime,i,j);
                    
%                     if sum(BackData_k) && ~isempty(R.Settings.Smoothing.Method)
%                         FiltBackData_k = filtTraj((BackData)',{'butter';R.Settings.Smoothing.CutFreq; 44});
%                         for k=1:size(RawDataPressure.Seat,2)
%                             DataPressure.Back(k).Values(i,j)= FiltBackData_k(k);
%                         end
%                     end
                end
            end
            
        end
        function DataPressure_k = filtAndsyncroMotFramePData(R,Data_k,DataPressure_k,MotionTime,PressureTime,Row,Col)
            % Si todos los valores no son cero filtrar y calcular el valor en el tiempos del motion file
            if sum(Data_k)>0
                if ~isempty(R.Settings.Smoothing.Method)
                    FiltData_k = filtTraj((Data_k)',{'butter';R.Settings.Smoothing.CutFreq; 44});
                    for i=1:length(FiltData_k)
                        if FiltData_k(i) < 0.005
                            FiltData_k(i) = 0;
                        end
                    end
                else
                    FiltData_k = Data_k';
                end
                % Calcular los valores del pressure map linealizando para los tiempos del motion file
                MTime = 2;
                Endfor = 0;
                Data_MTime(1) = FiltData_k(1);
                for k=2:length(MotionTime)
                    PTime = k;
                    if Endfor == 1
                        break
                    end
                    while MotionTime(MTime) < PressureTime(PTime) && Endfor == 0
                        % PVal_x = PValIni + (PValFin-PValIni)*((Tx-TIni)/(TFin-TIni))
                        Data_MTime(MTime) = FiltData_k(PTime-1) + (FiltData_k(PTime)-FiltData_k(PTime-1))*((MotionTime(MTime)-PressureTime(PTime-1))/((PressureTime(PTime)-PressureTime(PTime-1))));
                        if MTime < length(MotionTime)
                            MTime = MTime+1;
                        else
                            Endfor = 1;
                        end
                    end
                    
                end
                % Añadir los valores en la forma de estructura del Pressure map
                for k=1:length(MotionTime)
                    DataPressure_k(k).Values(Row,Col)= Data_MTime(k);
                end
                % Si son cero se les da valor cero
            else
                for k=1:length(MotionTime)
                    DataPressure_k(k).Values(Row,Col)= 0;
                end
            end
        end
        function fillExperimentValues(R,GuidedCoords)
            % filz z with guided variables
            % initialize
            NCoords = 1;
            NMarkers = size(R.Markers);
            % fill z with markers
            for i=1:NMarkers
                R.z{NCoords,1}      = R.Markers(i).CoordName{1};
                R.PoszInq (NCoords) = R.Markers(i).PosInq;
                R.Markers(i).PosInz = NCoords;
                R.Wm{NCoords,1}     = R.Markers(i).Wm;
                R.z{NCoords+1,1}    = R.Markers(i).CoordName{2};
                R.PoszInq(NCoords+1)= R.Markers(i).PosInq +1;
                R.Wm{NCoords+1,1}   = R.Markers(i).Wm;
                R.z{NCoords+2,1}    = R.Markers(i).CoordName{3};
                R.PoszInq(NCoords+2)= R.Markers(i).PosInq +2;
                R.Wm{NCoords+2,1}   = R.Markers(i).Wm;
                NCoords = NCoords +3;
            end
            % fill z with guided variables
            NGuidedCoords = size(GuidedCoords,1);
            for i=1:NGuidedCoords
                R.z{NCoords,1} = GuidedCoords{i};
                R.PoszInq (NCoords) = R.getPoszInq(GuidedCoords{i});
                R.Wm{NCoords,1} = GuidedCoords{i,2};
                NCoords = NCoords +1;
            end
            
        end
        function getMarkerGuidedValues(R,MarkerNames,MarkerCoords,GuidedCoords,GuidedMarkers)
            NMarkers  = size(R.Markers,1);
            NGuidedCoords = size(GuidedCoords,1); 
            Marker_t = [];
            MarkerName_i = {};
            GuidedCoords_t = [];
            % fill marker trajectories defined by the user
            for i=1:NMarkers
                MarkerName = R.Markers(i).Name;
                if ~isempty(R.Markers(i).MarkerTrajectory)
                    Markeri_t  = R.Markers(i).MarkerTrajectory;
                else
                    Markeri_t  = R.getCoordMarkerGuided(MarkerNames,MarkerCoords,MarkerName);
                end
                Marker_t =[Marker_t,Markeri_t];
                MarkerName_i = [MarkerName_i;MarkerName];
            end
            % fill gaps & check marker trajectories (and display report)
            Marker_t = R.checkMarkers(MarkerName_i,Marker_t);
            % fill guided coord trajectories defined by the user
            for i=1:NGuidedCoords
                GuidedCoordsi_t = GuidedCoords{i,3};
                GuidedCoords_t = [GuidedCoords_t,GuidedCoordsi_t];
            end
            R.g_t = [Marker_t, GuidedCoords_t]; 
        end
        function Marker_t = getCoordMarkerGuided(R,MarkerNames, MarkCoords, MarkerName)
            % get skin-marker index in cell MarkerNames
            SmarkIndex = R.getIndex(MarkerNames,MarkerName); % it can be a number or an empty array
            
            % sizes
            nFrames = size(MarkCoords,1);
            
            % check if marker is in motion file
            if isempty(SmarkIndex)
                % If marker IS NOT in motion file
                str = sprintf(['  WARNING: marker "',MarkerName,'" is not contained in motion file. Trajectory filed with "NaN,NaN,NaN"']);
                fprintf(R.ExperLogFileId,'%s\n', str);
                str = sprintf(['  WARNING: marker "',MarkerName,'" is not contained in motion file.']);
                disp(str);

                % fill marker trajectory with NaN's
                Marker_t = NaN * zeros(nFrames,3);
                
            else
                % If marker IS in motion file
                
                % get position of marker 'MarkerName' for all frames
                Smark_x_t = MarkCoords(:,3*SmarkIndex-2);
                Smark_y_t = MarkCoords(:,3*SmarkIndex-1);
                Smark_z_t = MarkCoords(:,3*SmarkIndex);
                
                % store data for each coord in a single array (nSamples x 3)
                Marker_t = [Smark_x_t, Smark_y_t, Smark_z_t];
            end
        end
        function IndexElement = getIndex(R,ElementArray, ElementName)
            % Find body properties
            nElements = size(ElementArray, 1);
            IndexElement = [];
            
            for i = 1 : nElements
                if strcmp(ElementName, ElementArray{i, 1}) == 1
                    IndexElement = i;
                end
            end
        end
        function PoszInq = getPoszInq(R,GuidedCoordName)
            Nq = size(R.Subject.q,1);
            PoszInq = [];
            for i=1:Nq
                if strcmpi(GuidedCoordName,R.Subject.q{i})
                    PoszInq = i;
                    break;
                end
            end
            if isempty(PoszInq)
                error('Guided coordinate %s does not belong to the model',GuidedCoordName);
            end
        end
        function initq0(R,GuidedFile)
            % INITTq0 fills q0. There are three ways to fill q0.
            % 1-With the data provided by the user
            % 2-With the motion. Using Söderkvist.
            % 3-With the subject parameters file data.
            
            %----------------------------------------------------------------
            % 1-q0 is filled with dates provided by the user in guided file
            %----------------------------------------------------------------
            if ~isempty(R.q0User)
                Nq0User = size(R.q0User,2);
                Nq = size(R.Subject.q,1);
                if Nq ~= Nq0User
                    error(['The number of variables in q0 is ',num2str(Nq0User),' but the user provides ',num2str(Nq),' in file ',GuidedFile,'.m']);
                end
                for i=1:Nq0User
                    R.q0(i,1) = R.q0User(i);
                end
            
            else
                %----------------------------------------------------------------
                % 2-q0 is filled with motion file
                %----------------------------------------------------------------
                NotInitSegs = R.initq0Soderkvist;
                
                % If some segments cannot be initialized send an error
                % For DHErgo The next 7 lines must be commented
                if ~isempty(NotInitSegs)
                    SegmentName = R.Subject.Segments(NotInitSegs.SegIndex).Name;
                    error([' Segment "',SegmentName,'" -> Initialization error: not possible to estimate segment initial position\n',...
                        '      1) Check if segment has enough markers in the MODEL definition file\n',...
                        '      2) Check if segment has enough markers in the MOTION definition file\n',...
                        ]);
                end

                
                %----------------------------------------------------------------
                % 3-q0 Uninitialized segments are obtained with angles values of 
                %      subject parameter file (ONLY POSSIBLE FOR DHErgo)
                %----------------------------------------------------------------
                % For DHErgo use the following code
                if ~isempty(NotInitSegs)
                    NSegments = size(NotInitSegs);
                    for i=1:NSegments
                        SegIndex = NotInitSegs(i).SegIndex;
                        NVectors = size(R.Subject.Segments(SegIndex).LocalVectors,1);
                        NPoints  = size(R.Subject.Segments(SegIndex).LocalPoints,1);
                        NMarkers = size(R.Subject.Segments(SegIndex).LocalMarkers,1);
                        for j=1:NVectors
                            R.Subject.Segments(SegIndex).LocalVectors(j).Vector.GlobalCoord = R.Subject.Segments(SegIndex).Glob_R_Seg * R.Subject.Segments(SegIndex).LocalVectors(j).LocCoord;
                        end
                        for j=1:NPoints
                            R.Subject.Segments(SegIndex).LocalPoints(j).Point.GlobalCoord = ...
                                    R.Subject.Segments(SegIndex).Glob_R_Seg * R.Subject.Segments(SegIndex).LocalPoints(j).LocCoord + ...
                                    R.Subject.Segments(SegIndex).LocalPoints(1).Point.GlobalCoord;
                        end
                        for j=1:NMarkers
                            R.Subject.Segments(SegIndex).LocalMarkers(j).Point.GlobalCoord = ...
                                    R.Subject.Segments(SegIndex).Glob_R_Seg * R.Subject.Segments(SegIndex).LocalMarkers(j).LocCoord + ...
                                    R.Subject.Segments(SegIndex).LocalPoints(1).Point.GlobalCoord;
                        end
                    end
                    
                end
                %----------------------------------------------------------------
                % 4-q0 is filled with dates of subject paramenter file
                %----------------------------------------------------------------
                % If sub par file has data in local coordinates the global coordinates must been calculated
                if isempty(R.Subject.Segments(2).LocalPoints(1).Point.GlobalCoord)
                    NSegments = size(R.Subject.Segments,1);
                    for i=1:NSegments
                        NVectors = size(R.Subject.Segments(i).LocalVectors,1);
                        NPoints  = size(R.Subject.Segments(i).LocalPoints,1);
                        NMarkers = size(R.Subject.Segments(i).LocalMarkers,1);
                        % Vectors in global coord are directly Glob_R_Seg * Seg_Vec
                        for j=1:NVectors
                            if R.Subject.Segments(i).Fixed ~= 1 % If is NOT Ground
                            %if ~strcmpi(R.Subject.Segments(i).Name,'Ground')
                                R.Subject.Segments(i).LocalVectors(j).Vector.GlobalCoord = R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalVectors(j).LocCoord;
                            end
                            
                        end
                        % Points in global coord are Glob_Pos = Glob_R_Seg * Seg_Pos + Glob_Pos_Origin
                        for j=1:NPoints
%                             WHEN RAMSIS SUB PAR FILE IS IN LOCAL COORDINATES                              
%                             if strcmpi(R.Subject.Segments(i).Name,'Ground')
%                                 R.Subject.Segments(i).LocalPoints(j).Point.GlobalCoord = R.Subject.Segments(i).LocalPoints(j).LocCoord;
%                             elseif ~strcmpi(R.Subject.Segments(i).Name,'Ground')
%                                 R.Subject.Segments(i).LocalPoints(j).Point.GlobalCoord = ...
%                                     R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalPoints(j).LocCoord + ...
%                                     R.Subject.Segments(i).LocalPoints(1).Point.GlobalCoord;
%                             end

                            % In the first segment the origin position is the first marker position
                            if strcmpi(R.Subject.Segments(i).Name,'BEC')||strcmpi(R.Subject.Segments(i).Name,'HM50KR_posture')%||...
                                Marker = 0;
                                for k=1:size(R.Subject.Segments(i).LocalMarkers,1)
                                    PosInqBecMarker = R.Subject.Segments(i).LocalMarkers(k).Point.PosInq;
                                    for l=1:length(R.PoszInq)
                                        if (R.PoszInq(l)==PosInqBecMarker)
                                            Marker = l;
                                            break;
                                        end
                                    end
                                    if Marker > 0
                                        break;
                                    end
                                end
                                if Marker == 0
                                    Marker =1;
                                    warning(['There are not markers in body ',R.Subject.Segments(i).Name, '. We choose the first marker to initial position in q0'])
                                end
                                Origin = [R.g_t(1,Marker);R.g_t(1,Marker+1);R.g_t(1,Marker+2)];
                                R.Subject.Segments(i).LocalPoints(j).Point.GlobalCoord = ...
                                    R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalPoints(j).LocCoord + Origin;
                            % WHEN RAMSIS SUB PAR FILE IS IN LOCAL COORDINATES                              
                            elseif R.Subject.Segments(i).Fixed == 1 % If is Ground
                            %elseif strcmpi(R.Subject.Segments(i).Name,'Ground')                                
                                R.Subject.Segments(i).LocalPoints(j).Point.GlobalCoord = R.Subject.Segments(i).LocalPoints(j).LocCoord;
                            elseif R.Subject.Segments(i).Fixed ~= 1 % If is NOT Ground                                
                            %elseif ~strcmpi(R.Subject.Segments(i).Name,'Ground')
                                R.Subject.Segments(i).LocalPoints(j).Point.GlobalCoord = ...
                                    R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalPoints(j).LocCoord + ...
                                    R.Subject.Segments(i).LocalPoints(1).Point.GlobalCoord;
                            
                            % In PAM model the origin has been located
                            elseif strcmpi(R.Subject.Segments(i).Name,'Left_leg')
                                OrLeftLeg = R.Subject.Segments(i).LocalPoints(1).Point.GlobalCoord -...
                                    R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalPoints(1).LocCoord;
                                if isempty(R.Subject.Segments(i).LocalPoints(j).Point.GlobalCoord)
                                    R.Subject.Segments(i).LocalPoints(j).Point.GlobalCoord = ...
                                        R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalPoints(j).LocCoord + OrLeftLeg;
                                end
                            elseif strcmpi(R.Subject.Segments(i).Name,'Left_calf')
                                OrLeftCalf = R.Subject.Segments(i).LocalPoints(1).Point.GlobalCoord -...
                                    R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalPoints(1).LocCoord;
                                if isempty(R.Subject.Segments(i).LocalPoints(j).Point.GlobalCoord)
                                    R.Subject.Segments(i).LocalPoints(j).Point.GlobalCoord = ...
                                        R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalPoints(j).LocCoord + OrLeftCalf;
                                end
                            elseif strcmpi(R.Subject.Segments(i).Name,'Left_foot')
                                OrLeftFoot = R.Subject.Segments(i).LocalPoints(1).Point.GlobalCoord;
                                if isempty(R.Subject.Segments(i).LocalPoints(j).Point.GlobalCoord)
                                    R.Subject.Segments(i).LocalPoints(j).Point.GlobalCoord = ...
                                        R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalPoints(j).LocCoord + OrLeftFoot;
                                end
                            end
                        end
                        % Markers in global coord are Glob_Pos = Glob_R_Seg * Seg_Pos + Glob_Pos_Origin
                        for j=1:NMarkers
                            %                         if strcmpi(R.Subject.Segments(i).Name,'BEC')
                            if strcmpi(R.Subject.Segments(i).Name,'BEC')||strcmpi(R.Subject.Segments(i).Name,'HM50KR_posture')%||...
%                                     strcmpi(R.Subject.Segments(i).Name,'Left_leg')
                                Origin = [R.g_t(1,Marker);R.g_t(1,Marker+1);R.g_t(1,Marker+2)];
                                R.Subject.Segments(i).LocalMarkers(j).Point.GlobalCoord = ...
                                    R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalMarkers(j).LocCoord + Origin;
                                
                                % RAMSIS en Locales
                            elseif R.Subject.Segments(i).Fixed ~= 1 % If is NOT Ground                                
                            %elseif ~strcmpi(R.Subject.Segments(i).Name,'Ground')
                                R.Subject.Segments(i).LocalMarkers(j).Point.GlobalCoord = ...
                                    R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalMarkers(j).LocCoord + ...
                                    R.Subject.Segments(i).LocalPoints(1).Point.GlobalCoord;
                            elseif strcmpi(R.Subject.Segments(i).Name,'Left_leg')
%                                 OrLeftLeg = R.Subject.Segments(i).LocalPoints(1).Point.GlobalCoord +...
%                                     R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalPoints(2).LocCoord;
                                R.Subject.Segments(i).LocalMarkers(j).Point.GlobalCoord = ...
                                    R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalMarkers(j).LocCoord + OrLeftLeg;
                            elseif strcmpi(R.Subject.Segments(i).Name,'Left_calf')
%                                 OrLeftCalf = R.Subject.Segments(i).LocalPoints(1).Point.GlobalCoord +...
%                                     R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalPoints(2).LocCoord;
                                R.Subject.Segments(i).LocalMarkers(j).Point.GlobalCoord = ...
                                    R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalMarkers(j).LocCoord + OrLeftCalf;
                            elseif strcmpi(R.Subject.Segments(i).Name,'Left_foot')
%                                 OrLeftFoot = R.Subject.Segments(i).LocalPoints(1).Point.GlobalCoord;
                                R.Subject.Segments(i).LocalMarkers(j).Point.GlobalCoord = ...
                                    R.Subject.Segments(i).Glob_R_Seg * R.Subject.Segments(i).LocalMarkers(j).LocCoord + OrLeftFoot;
                            end
                        end
                    end
                end
                NSegments = size(R.Subject.Segments,1);
                NRelPoints = size(R.Subject.PointsRel,1);
                % Inizialize of the Rel points
                for i =1:NSegments
                    if ~isempty (R.Subject.Segments(i).LocalPointsRel)
                        R.Subject.Segments(i).LocalPointsRel(1).Point.GlobalCoord = R.Subject.Segments(i).LocalPoints(2).LocCoord+R.Subject.Segments(i+1).LocalPoints(2).LocCoord;
                        % %                     R.Subject.Segments(i).LocalPointsRel(1).Point.GlobalCoord = R.Subject.Segments(i).LocalPoints(2).LocCoord+...
                        %                         [R.Subject.Segments(i+1).LocalVectors(1).LocCoord, R.Subject.Segments(i+1).LocalVectors(2).LocCoord,R.Subject.Segments(i+1).LocalPoints(2).LocCoord]*R.Subject.Segments(i+1).LocalPoints(2).LocCoord;
                    end
                end
                %             for i = 1: NRelPoints
                %                 NameRelPoint = R.Subject.PointsRel(i).Name(1:end-2);
                %                 PointPos = getVecIndex(NameRelPoint,R.Subject.Points);
                %                 R.Subject.PointsRel(i).GlobalCoord = [0;0;0];
                %             end
                
                
                % --------------------------------
                % fill q0 with problem variables
                % --------------------------------
                R.q0 = [];
                % initialize
                NPoints  = size(R.Subject.Points,1);
                NVectors = size(R.Subject.Vectors,1);
                NAngles  = size(R.Subject.Angles,1);
                NMarkers = size(R.Subject.Markers,1);
                % fill q0 with points
                for i = 1:NPoints
                    % Test if point is not fixed
                    if R.Subject.Points(i).Fixed == 0
                        PosInq = R.Subject.Points(i).PosInq;
                        R.q0(PosInq,1)   = R.Subject.Points(i).GlobalCoord(1);
                        R.q0(PosInq+1,1) = R.Subject.Points(i).GlobalCoord(2);
                        R.q0(PosInq+2,1) = R.Subject.Points(i).GlobalCoord(3);
                    end
                end
                % fill q0 with vectors
                for i = 1:NVectors
                    % Test if vector is not fixed
                    if R.Subject.Vectors(i).Fixed == 0
                        PosInq = R.Subject.Vectors(i).PosInq;
                        R.q0(PosInq,1)   = R.Subject.Vectors(i).GlobalCoord(1);
                        R.q0(PosInq+1,1) = R.Subject.Vectors(i).GlobalCoord(2);
                        R.q0(PosInq+2,1) = R.Subject.Vectors(i).GlobalCoord(3);
                    end
                end
                for i = 1:NRelPoints
                    PosInq = R.Subject.PointsRel(i).PosInq;
                    R.q0(PosInq,1)   = R.Subject.PointsRel(i).GlobalCoord(1);
                    R.q0(PosInq+1,1) = R.Subject.PointsRel(i).GlobalCoord(2);
                    R.q0(PosInq+2,1) = R.Subject.PointsRel(i).GlobalCoord(3);
                end
                % fill q0 with guided angles
                for i = 1:NAngles
                    if(strcmpi(R.Subject.Angles(i).Joint.Type,'SPH'))
                        PosInq = R.Subject.Angles(i).PosInq(1);
                        AngleName1 = R.Subject.Angles(i).Name1;
                        AngleName2 = R.Subject.Angles(i).Name2;
                        AngleName3 = R.Subject.Angles(i).Name3;
                        PosOf_     = findstr('_',AngleName1);
                        AngleType1 = AngleName1(PosOf_+1:end);
                        AngleType2 = AngleName2(PosOf_+1:end);
                        AngleType3 = AngleName3(PosOf_+1:end);
                        R.q0(PosInq,1)   = R.Subject.Angles(i).(AngleType1);
                        R.q0(PosInq+1,1) = R.Subject.Angles(i).(AngleType2);
                        R.q0(PosInq+2,1) = R.Subject.Angles(i).(AngleType3);
                    elseif(strcmpi(R.Subject.Angles(i).Joint.Type,'UNI'))
                        PosInq = R.Subject.Angles(i).PosInq(1);
                        AngleName1 = R.Subject.Angles(i).Name1;
                        AngleName2 = R.Subject.Angles(i).Name2;
                        PosOf_     = findstr('_',AngleName1);
                        AngleType1 = AngleName1(PosOf_+1:end);
                        AngleType2 = AngleName2(PosOf_+1:end);
                        R.q0(PosInq,1)   = R.Subject.Angles(i).(AngleType1);
                        R.q0(PosInq+1,1) = R.Subject.Angles(i).(AngleType2);
                    elseif(strcmpi(R.Subject.Angles(i).Joint.Type,'REV'))
                        PosInq = R.Subject.Angles(i).PosInq;
                        AngleName1 = R.Subject.Angles(i).Name1;
                        PosOf_     = findstr('_',AngleName1);
                        AngleType1    = AngleName1(PosOf_+1:end);
                        if isempty(PosOf_)
                            AngleType1 = 'a1';
                        end
                        R.q0(PosInq,1) = R.Subject.Angles(i).(AngleType1);
                    else
                        error('Incorrect type of Joint')
                    end
                end
                % fill q0 with markers
                for i=1:NMarkers
                    PosInq = R.Subject.Markers(i).PosInq;
                    R.q0(PosInq,1)   = R.Subject.Markers(i).GlobalCoord(1);
                    R.q0(PosInq+1,1) = R.Subject.Markers(i).GlobalCoord(2);
                    R.q0(PosInq+2,1) = R.Subject.Markers(i).GlobalCoord(3);
                end
            end
        end
        function NotInitSegs = initq0Soderkvist(R)
            % INITq0SODERKVIST initicialize q0 with the data of motion files by Soderkvist
            
            % variable initialization
            NSegments = size(R.Subject.Segments,1);
            NotInitSegs = [];
            Lt3Points = [];
            AllSegInitInFrame = 0;
            Frame = 1;
            % Are calculated the segments that are not posible calculate in previous frames. Markers NaN
            while AllSegInitInFrame == 0
                % Calc the optimal R and d for each segment with markers positions
                for i=1:NSegments
                    if isempty(NotInitSegs)
                        SegIndex = i;
                    else
                        SegIndex = NotInitSegs(i).SegIndex;
                    end
                    if R.Subject.Segments(SegIndex).Fixed ~= 1 % If is NOT Ground                                                    
                    %if ~strcmp(R.Subject.Segments(SegIndex).Name,'Ground')
                        NMarkers = size(R.Subject.Segments(SegIndex).LocalMarkers,1);
                        Glob_Pos_Markers = [];
                        Seg_Pos_Markers = [];
                        for j=1:NMarkers
                            PosInzMarker = R.Subject.Segments(SegIndex).LocalMarkers(j).Point.PosInz;
                            Glob_Pos_Marker_i = R.g_t(Frame,PosInzMarker:PosInzMarker+2);
                            % Check if experimental marker have a non NaN value
                            if ~isnan(Glob_Pos_Marker_i(1))
                                Glob_Pos_Markers = [Glob_Pos_Markers,Glob_Pos_Marker_i'];
                                Seg_Pos_Markers = [Seg_Pos_Markers,R.Subject.Segments(SegIndex).LocalMarkers(j).LocCoord];
                            end
                        end
                        % Three marker are necesary to calc optimal pose.
                        if size(Glob_Pos_Markers,2) > 2
                            [Glob_R_Seg, Glob_Pos_Or] = R.Subject.Segments(SegIndex).calOptPose(Glob_Pos_Markers, Seg_Pos_Markers);
                            R.Subject.Segments(SegIndex).calcGlobWithOptPose(Glob_R_Seg,Glob_Pos_Or);
                        % The Segment Index and the markes value are stored.
                        else
                            Lt3Points_i.SegIndex = SegIndex;
                            Lt3Points_i.Glob_Pos_Markers = Glob_Pos_Markers;
                            Lt3Points_i.Seg_Pos_Markers  = Seg_Pos_Markers;
                            Lt3Points_i.OptPose = 0;
                            Lt3Points = [Lt3Points;Lt3Points_i];
                        end
                    end
                end
                % If markesr are not enough to calculate opt pose segments points will be used
                NSegments = size(Lt3Points,1);
                if NSegments == 0
                    % Inicializated correctly
                    AllSegInitMarkers = 1;
                    AllSegInitInFrame = 1;
                    NotInitSegs = [];
                else
                    AllSegInitMarkers = 0;
                end
                % Are calculated the segments with joint points
                while AllSegInitMarkers == 0 
                    for i=1:NSegments
                        SegIndex = Lt3Points(i).SegIndex;
                        Glob_Pos_Points = Lt3Points(i).Glob_Pos_Markers;
                        Seg_Pos_Points  = Lt3Points(i).Seg_Pos_Markers;
                        
                        NPoints  = size(R.Subject.Segments(SegIndex).LocalPoints,1);
                        NMarkers = size(Glob_Pos_Points,2);
                        % is necessary that the sum of Markers and Points is three
                        if NPoints + NMarkers > 2
                            for j=1:NPoints
                                % The point value is valid when the other segment has been calculated
                                if (R.Subject.Segments(SegIndex).LocalPoints(j).Point.OptPose == 1 || ...
                                   R.Subject.Segments(SegIndex).LocalPoints(j).Point.Fixed == 1) && ...
                                   ~strcmpi(R.Subject.Segments(SegIndex).LocalPoints(j).Name,'GBRK')    
                                    Glob_Pos_Point_j = R.Subject.Segments(SegIndex).LocalPoints(j).Point.GlobalCoord;
                                    Seg_Pos_Pointj = R.Subject.Segments(SegIndex).LocalPoints(j).LocCoord;
                                    Glob_Pos_Points = [Glob_Pos_Points,Glob_Pos_Point_j];
                                    Seg_Pos_Points  = [Seg_Pos_Points,Seg_Pos_Pointj];
                                end
                            end
                            % Three points are necesary to calculate optimal pose
                            if size(Glob_Pos_Points,2) > 2
                                [Glob_R_Seg, Glob_Pos_Or] = R.Subject.Segments(SegIndex).calOptPose(Glob_Pos_Points, Seg_Pos_Points);
                                R.Subject.Segments(SegIndex).calcGlobWithOptPose(Glob_R_Seg,Glob_Pos_Or);
                                Lt3Points(i).OptPose = 1;
                            end
                        end
                    end
                    AllSegInitMarkers = 1;
                    NewLt3Points = [];
                    % Check if segments has not been calculated
                    for i=1:NSegments
                        if Lt3Points(i).OptPose == 0
                            AllSegInitMarkers = 0;
                            NewLt3Points_i.SegIndex = Lt3Points(i).SegIndex;
                            NewLt3Points_i.Glob_Pos_Markers = Lt3Points(i).Glob_Pos_Markers;
                            NewLt3Points_i.Seg_Pos_Markers  = Lt3Points(i).Seg_Pos_Markers;
                            NewLt3Points_i.OptPose = Lt3Points(i).OptPose;
                            NewLt3Points = [NewLt3Points;NewLt3Points_i];
                        end
                    end
                    % Check if is posible to calculate more segments with other segments points
                    if NSegments == size(NewLt3Points,1)
                        % Not possible to calculate more segments
                        AllSegInitMarkers = 1;
                        AllSegInitInFrame = 0;
                        NSegments = size(NewLt3Points,1);
                        NotInitSegs = NewLt3Points;
                        NewLt3Points = [];
                        Lt3Points = [];
                    else
                        if isempty(NewLt3Points)
                            % All segments correctly inicializated with other segment points
                            AllSegInitInFrame = 1;
                            NotInitSegs = [];
                        else
                            % Some segments correcly inicializated with other segments points check again other segments
                            NSegments = size(NewLt3Points,1);
                            Lt3Points = NewLt3Points;
                        end
                    end
                    
                end
                % If for the i Frame it is not posible check with the next frame
                Frame = Frame +1;
                % If frame is greater than a value a warning is showed.
                if Frame > (R.Settings.Interpolation.NInterFrames+5)
                    AllSegInitInFrame = 1;
                    str = ['  WARNING: Initialization not possible using the motion file. Not enough data in the first ',num2str(R.Settings.Interpolation.NInterFrames+5),' frames.'];
                    disp(str);
                    fprintf(R.ExperLogFileId, '%s\n', str);
                end
                    
            end
        end
        function plotPressureMap(R,DataPressureMap)
%             for i=1:40:size(DataPressureMap.Seat,2)
%                 figure('Name',[R.ExperimentName,'_SeatF',num2str(i)],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%                 bar3(DataPressureMap.Seat(i).Values)
%                 figure('Name',[R.ExperimentName,'_SeatBackF',num2str(i)],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%                 bar3(DataPressureMap.Back(i).Values)
%             end
            SeatData = [];
            BackData = [];
            t=[];
            for i=1:size(DataPressureMap.Seat,2)
                t=[t,i];
                SeatData = [SeatData,DataPressureMap.Seat(i).Values(5,25)];
                BackData = [BackData,DataPressureMap.Back(i).Values(22,25)];
            end
            figure('Name',[R.ExperimentName,'_Seat5_25'],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
            plot(t,SeatData)
            figure('Name',[R.ExperimentName,'_SeatBack22_25'],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
            plot(t,BackData)

        end
        function [Frequency, MarkerNames, MarkerCoords, DataCSV] = readCSV(R)
            % READCSV reads a motion capture file with .CSV format
            %   [Frequency, MarkerNames, MarkerCoords] = readCSV()
            %   Outputs:
            %     + Frequency is the sampling frequency in Hz.
            %     + MarkerNames is a cell (nMarker x 1) with the names of the skin-markers. Each
            %       element of the cell is a string.
            %     + MarkerCoords is a double array(nFrames x 3*nMarker) with the 3D coordinates of
            %       skin-markers. The coordinates x,y,z of skin-marker i are in position
            %       MarkerCoords(3*i-2,3*i-1,3*i)
            
            % Path and File
            InputPath = R.Motion.Path;
            InputFile = R.Motion.File;
            % Initialize variables
            DataCSV = [];
            % Initialize counters
            i = 1;  % This counter represents the read row of the file
            p = 1;  % This counter represents the row of the matrix where the numeric information is inserted
            r = 1;  % This counter represents the rows of the analog
            Analog = 0; % This counter represents when start the analog part
            
            % ----------------------------------------------------------------------------
            % Read the markers 3D coord. from the 'InputFile'
            % ----------------------------------------------------------------------------
            % Open the input file for reading
            FidInput = fopen([InputPath,InputFile]);
            if FidInput == -1
                error(['The file ',InputFile,' can not be opened in path (',InputPath,'). ']);
            end
            % 1) There are not heading lines. First line starts with keyword TRAJECTORIES or TRAJECTOIRES.
            % 2) When a skin-markers is missing the format is ",NaN,NaN,NaN"
            Frequency_LineIndex    = 2;
            MarkerNames_LineIndex  = 3;
            MarkerCoords_LineIndex = 5;
            AnalogFrecuency_LineIndex = 0;
            AnalogNames_LineIndex = 0;
            AnalogUnits_LineIndex = 0;
            AnalogData_LineIndex = 0;
            

            % While the reading file hasn't reached the end of file keep calculating
            while feof(FidInput) == 0
                
                % Get the whole line
                Line = fgets(FidInput);
                
                % In line 'Frequency_LineIndex' get the Frequency
                if i == Frequency_LineIndex
                    % Different formats are found: 50.000Hz  "50.000;Hz"  50,Hz
                    % As normally a frequency of 50Hz is used I consider this value fixed!!!!!
                    k = strfind(Line, ',');
                    Frequency = str2num(Line(1:(k-1)));
%                     Frequency = 100;
                    % In line 'MarkerNames_LineIndex' get the specific number of the measured markers
                elseif i == MarkerNames_LineIndex
                    MarkerNames = {};
                    nMarker = 1;
                    % delimiters: comma (ASCII 44), new line (ASCII 10 - \n), carriage return (ASCII 13 - \r)
                    [MarkerName_i, Remainder] = strtok(Line, [44 10 13]);
                    while isempty(MarkerName_i) == 0
                        MarkerNames{nMarker,1} = upper(MarkerName_i);
                        nMarker = nMarker + 1;
                        [MarkerName_i, Remainder] = strtok(Remainder, [44 10 13]);
                    end
                    
                    % From line MarkerCoords_LineIndex till the end get the 3D coordinates of the measured markers
                elseif i >= MarkerCoords_LineIndex && Analog == 0
                if 	strcmpi(deblank(Line),'ANALOG') || strcmp(deblank(Line),'FORCES')
                    AnalogFrecuency_LineIndex = i+1;
                    AnalogNames_LineIndex = i+2;
                    AnalogUnits_LineIndex = i+3;
                    AnalogData_LineIndex = i+4;
                    Analog=1;
                end
                    k = strfind(Line, ',,');    % Returns the positions of the ',,' characters. These commas represent an occluded coordinate
                    % The format of NaN is not unique in CSV files. There are at least two version: Nan and NaN.
                    % Here al possible combinations are transformed to NaN in order to be compatible with Matlab.
                    Line = strrep(Line,'nan','NaN');
                    Line = strrep(Line,'Nan','NaN');
                    Line = strrep(Line,'NAn','NaN');
                    Line = strrep(Line,'NAN','NaN');
                    Line = strrep(Line,'naN','NaN');
                    Line = strrep(Line,'nAN','NaN');
                    
                    % If the line isn't empty fill the MarkerCoord matrix
                    Line2num = str2num(Line);
                    if norm(Line2num) ~= 0
                        MarkerCoords(p,:) = Line2num(2:end);   % The first column isn't of interest
                        p = p + 1;
                    end
                elseif i == AnalogFrecuency_LineIndex
                    k = strfind(Line, ',');
                    DataCSV.Analog.Frequency = str2num(Line(1:(k-1)));
                elseif i == AnalogNames_LineIndex
                    % delimiters: comma (ASCII 44), new line (ASCII 10 - \n), carriage return (ASCII 13 - \r)
                    [AnalogName_first, Remainder] = strtok(Line, [44 10 13]); % The first name is the keyword Sample#
                    [AnalogName_i, Remainder] = strtok(Remainder, [44 10 13]);
                    nAnalog = 1;
                    while isempty(AnalogName_i) == 0
                        AnalogNames{nAnalog,1} = AnalogName_i;
                        nAnalog = nAnalog +1;
                        [AnalogName_i, Remainder] = strtok(Remainder, [44 10 13]);
                    end
                    
                elseif i == AnalogUnits_LineIndex
                    nAnalog = nAnalog - 1;
                    [AnalogUnits_first, Remainder] = strtok(Line, [44 10 13]); % The first name is the keyword Units:
                    [AnalogUnits_i, Remainder] = strtok(Remainder, [44 10 13]);
                    for j=1:nAnalog
                         DataCSV.Analog.(AnalogNames{j}).Units = AnalogUnits_i;
                         [AnalogUnits_i, Remainder] = strtok(Remainder, [44 10 13]);
                    end
                elseif i >= AnalogData_LineIndex && Analog == 1
                    Line2num = str2num(Line);
                    if ~isempty(Line2num)
                        for j=1:nAnalog
                            DataCSV.Analog.(AnalogNames{j}).Values(r) = Line2num(j+1);
                        end
                        r = r + 1;
                    end
                end
                i = i + 1;
            end
            
            % Close the input file
            fclose(FidInput);
            
            % mm to m
            MarkerCoords = MarkerCoords/1000; % mm to m
            % remove the front of :
            for i=1:size(MarkerNames,1)
                Index = findstr(MarkerNames{i},':');
                if ~isempty(Index)
                    MarkerNames{i} = MarkerNames{i}(Index+1:end);
                end
            end
            if exist('DataCSV','var')==0 && ~strcmpi(R.Settings.Type,'IK')
                R.ErrorID = 1;
                str = (['  WARNING: The motion file ',InputFile, ' does not have forces section.']);
                disp(str);
                fprintf(R.ExperLogFileId, '%s\n', str);
            end
        end
        function [Frequency, MarkerNames, MarkerCoords] = readMAT(R)
            % READMAT read a .MAT file in which data from a .C3D has been stored
            % Path and File
            InputPath = R.Motion.Path;
            InputFile = R.Motion.File; % Es el C3D. Coger solo el filename y añadirle .mat

            AceptedExtension{1,1} = 'mat';
            checkFileAndPath(InputPath,InputFile,AceptedExtension);            

            [Filename, ~] = getFilenameAndExt(InputFile);
            
            C3D_DATA = load([InputPath,Filename,'.mat']);            
            Frequency = C3D_DATA.Frequency;
            MarkerNames = C3D_DATA.MarkerNames;
            MarkerCoords = C3D_DATA.MarkerCoords;
        end        
        function [Frequency, MarkerNames, MarkerCoords] = readC3D(R)
            % READC3D read a motion capture file in .C3D format using C3Dserver (Shriners Hospital, 2006)
            %   [Frequency,MarkerNames, MarkerCoords] = readC3D()
            %   Outputs:
            %     + MarkerNames is a cell (nMarkers x 1) with the names of the skin-markers. Each
            %       element of the cell is a string.
            %     + MarkerCoords is a double array(nFrames x 3*nMarkers) with the 3D coordinates of
            %       skin-markers. The coordinates x,y,z of skin-marker i are in position
            %       MarkerCoords(3*i-2,3*i-1,3*i)
            %     + Frequency is a double array (1x1) with the frequency in Hertz at which the motion
            %       was recorded.
            %     + Markers is a structure with nMarkers fields. The name of each field is the name of
            %       the marker. For example, for marker 'M4' it coordinates are Markers.M4, this is an
            %       array (nFramesx3) with M4 x,y,z coordinates for each frame
            
            % Path and File
            InputPath = R.Motion.Path;
            InputFile = R.Motion.File;
            
            % Activate C3Dserver as a COM object
            SERVER = c3dserver();
            
            % Load file C3D in C3Dserver
            openc3d(SERVER, 0, [InputPath, InputFile]);
            
            % Read marker coordinates
            Markers = get3dtargets(SERVER);
            
            % Read Video Frame Rate (VFR)
            VFR_Index = SERVER.GetParameterIndex('POINT', 'RATE');
            Frequency = double(SERVER.GetParameterValue(VFR_Index, 0));
            
            % close file and C3Dserver
            closec3d(SERVER);
            
            % Units of marker coordinates
            Units = Markers.units;
            % Remove field 'units', only markers are left
            Markers = rmfield(Markers, 'units');
            
            % Extraer nombres de los markers
            MarkerNames = fieldnames(Markers);
            nMarkers = length(MarkerNames);
            
            % Convert marker coordinates units to metres if necessary.
            if strcmpi(Units,'mm')
                for i=1:nMarkers
                    Markers.(MarkerNames{i}) = (1/1000) * Markers.(MarkerNames{i});
                end
                
            elseif ~strcmpi(Units,'m')
                error('Marker coordinate units are neither metres [m] nor milimetres [mm]. Appropriated conversion should be implemented');
            end
            
            % Crear matriz MarkerCoords array(nFrames x 3*nMarkers)
            MarkerCoords = [];
            for i=1:nMarkers
                MarkerCoords = [MarkerCoords, Markers.(MarkerNames{i})];
            end
            
            % Cut out of trajectories. In the 1st frame must be at least one marker (not all may be lost)
            % We also calculat the first frame where they appear all markers
            Index = 0;
            IndexAllVisible = 0;
            nFrames  = size(MarkerCoords,1);
            i = 1;
            while IndexAllVisible == 0  &&  i<=nFrames
                if ((sum(isnan(MarkerCoords(i,:))) < 3*nMarkers - 3) && (Index==0)) % At least 1 visible Marker
                    Index = i;
                elseif (sum(isnan(MarkerCoords(i,:))) == 0)                         % All markers visibles
                    IndexAllVisible = i;
                else
                    i = i + 1;
                end
            end
%             if IndexAllVisible == 0
%                 dispWarning('A frame where all markers are visibles does not exist');
%             end
            MarkerCoords = MarkerCoords(Index:end,:);
            MarkerCoords = double(MarkerCoords);
            for i=1:nMarkers
                Markers.(MarkerNames{i}) = double(Markers.(MarkerNames{i})(Index:end,:));
            end
            
        end
        function [Frequency, MarkerNames, MarkerCoords] = readMOT(R)
            % READMOT reads a motion capture file with .Mot format
            %   [Frequency, MarkerNames, MarkerCoords] = readMOT()
            %   Outputs:
            %     + Frequency is the sampling frequency in Hz.
            %     + MarkerNames is a cell (nMarker x 1) with the names of the skin-markers. Each
            %       element of the cell is a string.
            %     + MarkerCoords is a double array(nFrames x 3*nMarker) with the 3D coordinates of
            %       skin-markers. The coordinates x,y,z of skin-marker i are in position
            %       MarkerCoords(3*i-2,3*i-1,3*i)
            
            % Path and File
            InputPath = R.Motion.Path;
            InputFile = R.Motion.File;
            
            % Initialize counters
            i = 1;  % This counter represents the read row of the file
            p = 1;  % This counter represents the row of the matrix where the numeric information is inserted
            
            % ----------------------------------------------------------------------------
            % Read the markers 3D coord. from the 'InputFile'
            % ----------------------------------------------------------------------------
            % Open the input file for reading
            FidInput = fopen([InputPath,InputFile]);
            if FidInput == -1
                error(['The file ',InputFile,' can not be opened in path (',InputPath,'). ']);
            end
            % 1) There are not heading lines. First line starts with keyword TRAJECTORIES or TRAJECTOIRES.
            % 2) When a skin-markers is missing the format is ",NaN,NaN,NaN"
            Frequency_LineIndex    = 2;
            MarkerNames_LineIndex  = 3;
            MarkerCoords_LineIndex = 4;
            
            % Close the input file
            fclose(FidInput);
            
            FidInput = fopen([InputPath,InputFile]);
            % While the reading file hasn't reached the end of file keep calculating
            while feof(FidInput) == 0
                
                % Get the whole line
                Line = fgets(FidInput);
                
                % In line 'Frequency_LineIndex' get the Frequency
                if i == Frequency_LineIndex
                    % Different formats are found: 50.000Hz  "50.000;Hz"  50,Hz
                    k = strfind(Line, ',');
                    Frequency = str2num(Line(1:(k-1)));
                % In line 'MarkerNames_LineIndex' get the specific number of the measured markers
                elseif i == MarkerNames_LineIndex
                    MarkerNames = {};
                    nMarker = 1;
                    % delimiters: comma (ASCII 44), new line (ASCII 10 - \n), carriage return (ASCII 13 - \r)
                    [MarkerName_i, Remainder] = strtok(Line, [44 10 13]);
                    while isempty(MarkerName_i) == 0
                        MarkerNames{nMarker,1} = upper(MarkerName_i);
                        nMarker = nMarker + 1;
                        [MarkerName_i, Remainder] = strtok(Remainder, [44 10 13]);
                    end
                    
                % From line MarkerCoords_LineIndex till the end get the 3D coordinates of the measured markers
                elseif i >= MarkerCoords_LineIndex
                    k = strfind(Line, ',,');    % Returns the positions of the ',,' characters. These commas represent an occluded coordinate
                    % The format of NaN is not unique in CSV files. There are at least two version: Nan and NaN.
                    % Here al possible combinations are transformed to NaN in order to be compatible with Matlab.
                    Line = strrep(Line,'nan','NaN');
                    Line = strrep(Line,'Nan','NaN');
                    Line = strrep(Line,'NAn','NaN');
                    Line = strrep(Line,'NAN','NaN');
                    Line = strrep(Line,'naN','NaN');
                    Line = strrep(Line,'nAN','NaN');
                    
                    % If the line isn't empty fill the MarkerCoord matrix
                    Line2num = str2num(Line);
                    if norm(Line2num) ~= 0
                        MarkerCoords(p,:) = Line2num(2:end);   % The first column isn't of interest
                        p = p + 1;
                    end
                end
                i = i + 1;
            end
            
            % Close the input file
            fclose(FidInput);
            
            % mm to m
            MarkerCoords = MarkerCoords/1000; % mm to m
            
        end
        function readForces(R,Path,File)
             % Path and File
            InputPath = Path;
            InputFile = File;
            AceptedExtension{1,1} = 'for';
            checkFileAndPath(InputPath,InputFile,AceptedExtension);
%             % Initialize counters
            NLine = 1;  % This counter represents the read row of the file
           
            % ----------------------------------------------------------------------------
            % Read the forces from the 'InputFile'
            % ----------------------------------------------------------------------------
            % Open the input file for reading
            FidInput = fopen([InputPath,InputFile]);
            if FidInput == -1
                error(['The file ',InputFile,' can not be opened in path (',InputPath,'). ']);
            end
            % While the reading file hasn't reached the end of file keep calculating
            while feof(FidInput) == 0
                
                % Get the whole line
                Line = fgets(FidInput);
                if NLine == 1
                    Frequency = str2num(Line);
                elseif NLine == 2
                    NSegments = str2num(Line);
                    CSegment = 1;
                    Segments = [];
                elseif NLine > 2 && NLine <= (2+NSegments)
                    k = strfind(Line, ',');
                    if strcmp(Line(1:k(1)-1),'Pelvis')
                        SegIndex = getVecIndex('BEC',R.Subject.Segments);
                        if SegIndex == 0
                            SegIndex = getVecIndex('HM50KR_posture',R.Subject.Segments);
                        end
                        if SegIndex == 0
                            error('There is not Pelvis segment in the model')
                        end
                    elseif strcmp(Line(1:k-1),'LeftThigh')
                        SegIndex = getVecIndex('OSL',R.Subject.Segments);
                        if SegIndex == 0
                            SegIndex = getVecIndex('Left_leg',R.Subject.Segments);
                        end
                        if SegIndex == 0
                            error('There is not LeftThigh segment in the model')
                        end
                    elseif strcmp(Line(1:k-1),'LeftShank')
                        SegIndex = getVecIndex('USL',R.Subject.Segments);
                        if SegIndex == 0
                            SegIndex = getVecIndex('Left_calf',R.Subject.Segments);
                        end
                        if SegIndex == 0
                            error('There is not LeftShank segment in the model')
                        end
                    elseif strcmp(Line(1:k-1),'LeftFoot')
                        SegIndex = getVecIndex('FUL',R.Subject.Segments);
                        if SegIndex == 0
                            SegIndex = getVecIndex('Left_foot',R.Subject.Segments);
                        end
                        if SegIndex == 0
                            error('There is not LeftFoot segment in the model')
                        end
                    end
                    Seg = R.Subject.Segments(SegIndex);
                    Sys = str2num(Line(k(1)+1:k(2)-1));
                    if Sys~=0 && Sys~=1
                        error (['The keyword for the system ',Sys,' is incorrect keyword in the file ',InputPath,InputFile]);
                    end
                    Seg.F_Ext.Sys = Sys;
                    Segments = [Segments;Seg];
                    CSegment = CSegment+1;
                elseif NLine > (3+ NSegments)
                    k = strfind(Line,',');
                    Frame = str2num(Line(1:k(1)-1));
                    CSegment = str2num(Line(k(1)+1:k(2)-1));
                    PosX = str2num(Line(k(2)+1:k(3)-1))/1000;
                    PosY = str2num(Line(k(3)+1:k(4)-1))/1000;
                    PosZ = str2num(Line(k(4)+1:k(5)-1))/1000;
                    FX = str2num(Line(k(5)+1:k(6)-1));
                    FY = str2num(Line(k(6)+1:k(7)-1));
                    FZ = str2num(Line(k(7)+1:end));
                    Segments(CSegment).F_Ext.Value(Frame,:) = [FX,FY,FZ];
                    Segments(CSegment).F_Ext.Pos(Frame,:) = [PosX,PosY,PosZ];
                end
                NLine = NLine+1;
            end
            
        end
        function [Frequency, MarkerNames, MarkerCoords, NFrames] = readMotionFile(R)
            if isempty(R.Motion)
                Frequency = NaN;
                NFrames = NaN;
                MarkerNames = [];
                MarkerCoords = [];
                %str = sprintf('  WARNING: the motion file is not defined. Marker trajectories should be define in file where guided variables are defined');
                %disp(str);
                %fprintf(R.ExperLogFileId, '%s\n', str);
            elseif strcmpi(R.Motion.File(end-2:end),'CSV')
                str = sprintf(['  Reading the motion file:',R.Motion.File,'...']);
                if(R.Settings.Display == 1 || R.Settings.Display == 2), disp(str); end
                fprintf(R.ExperLogFileId, '%s\n', str);
                % read marker trajectories
                [Frequency, MarkerNames, MarkerCoords, DataCSV] = R.readCSV();
                R.RawMarkerCoords = MarkerCoords;
                R.RawMarkerNames = MarkerNames;
                % smoothing
                if strcmpi(R.Settings.Smoothing.Method, 'butter')
                    CutFreq = R.Settings.Smoothing.CutFreq;
                    str = ['  Filtering marker trajectories... Method: butter, Cut-off freq: ',num2str(CutFreq)];
                    if(R.Settings.Display == 1 || R.Settings.Display == 2), disp(str); end
                    fprintf(R.ExperLogFileId, '%s\n', str);
                    MarkerCoords = filtTrajNaN(MarkerCoords, {'butter'; CutFreq; Frequency});                    
                    
                elseif isempty(R.Settings.Smoothing.Method)
                    str = '  No filtering applied to marker trajectories.';
                    if(R.Settings.Display == 1 || R.Settings.Display == 2), disp(str); end
                    fprintf(R.ExperLogFileId, '%s\n', str);
                    
                else
                    str = ['  Filtering method "',R.Settings.Smoothing.Method,'" is unknown'];
                    fprintf(R.ExperLogFileId, '%s\n', str);
                    error(str);
                    
                end                
                if ~isempty(DataCSV)
                    ForceNames = fieldnames(DataCSV.Analog);
                    for i=2:size(ForceNames,1)
                        DataCSV.Analog.(ForceNames{i}).ValuesNoFiltrados = DataCSV.Analog.(ForceNames{i}).Values;
                        if ~isempty(R.Settings.Smoothing.Method)
%                             DataCSV.Analog.(ForceNames{i}).Values = filtTraj((DataCSV.Analog.(ForceNames{i}).Values)',{'butter';1; DataCSV.Analog.Frequency});
                              DataCSV.Analog.(ForceNames{i}).Values = filtTraj((DataCSV.Analog.(ForceNames{i}).Values)',{'butter';R.Settings.Smoothing.CutFreq; DataCSV.Analog.Frequency});

                        end
                        for j=1:size(DataCSV.Analog.(ForceNames{i}).ValuesNoFiltrados,2)
                            if (DataCSV.Analog.(ForceNames{i}).ValuesNoFiltrados(1,j) == 0)
                                DataCSV.Analog.(ForceNames{i}).Values(j,1) = 0;
                            end
                        end
                    end
                    R.storeForceValues(DataCSV,MarkerNames,MarkerCoords,Frequency);
                end
                DataCSV.MarkTraj.MarkerNames = MarkerNames;
                DataCSV.MarkTraj.MarkerCoords = MarkerCoords;
                DataCSV.MarkTraj.Frequency = Frequency;
                R.DataCSV = DataCSV;
                if exist(R.PressurePath,'dir')
                    R.Settings.ExType = 'CP';
%                     time = tic;
                    RawDataPressure = R.readPressureFile();
%                     [H, MI, S] = second2HMS(toc(time))
                    if ~isempty(RawDataPressure)
                        R.plotPressureMap(RawDataPressure);
                        MotFrames = (size(MarkerCoords,1)-1);
                        DataPressure = R.filtAndsyncroPData(RawDataPressure,MotFrames,Frequency);
                        % get the zero line in Seat data
                        NCol = size(DataPressure.Seat(1).Values,2);
                        NRow = size(DataPressure.Seat(1).Values,1);
                        Fr = length(DataPressure.Seat);
                        for i =1:Fr
                            CBreak = 0;
                            for j=1:NRow
                                if CBreak == 1
                                    break;
                                end
                                SumRow = sum(DataPressure.Seat(i).Values(j,:));
                                if SumRow > 0
                                    for k=1:NCol-1
                                        SumColk = DataPressure.Seat(i).Values(j,k) + DataPressure.Seat(i).Values(j,k+1);
                                        if SumColk>0 && SumColk == DataPressure.Seat(i).Values(j,k) && sum(DataPressure.Seat(i).Values(j,k+1:NCol)) > 0 && sum(DataPressure.Seat(i).Values(j:NRow,k+1)) == 0
                                            DataPressure.Seat(i).ZeroRow = j;
                                            DataPressure.Seat(i).ZeroCol = k+1;
                                            CBreak =1;
                                            break;
                                        end
                                    end
                                end
                            end
                        end
%                         PepF = [];
%                         for i=1:486
%                             PepF = [PepF;DataPressure.Seat(i).Values(21,8)];
%                         end
                        DataPressure.Units = RawDataPressure.Units;
                        R.plotPressureMap(DataPressure);
                        R.storePressureData(DataPressure,MarkerNames,MarkerCoords,Frequency);
                        R.DataCSV.PMap = DataPressure;
                    end
                end
%                 figure('Name','NoFilt F2z','NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%                 % figure settings
%                 grid on
%                 plot(tiempo,DataCSV.Analog.Fz2.ValuesNoFiltrados);
%                 figure('Name','Filt  F2z','NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%                 plot(tiempo,DataCSV.Analog.Fz2.Values);
%                 for i=1:size(MarkerNames,1)
%                     figure('Name',['Filt ',MarkerNames{i},'_x'],'NumberTitle','off','WindowStyle','docked'), axes('FontSize',12), hold on
%                     % figure settings
%                     grid on
%                     plot(tiempo,MarkerCoords(:,(3*i-2)));
%                 end
%                 
%                 MarkerCoords = R.checkMarkers(MarkerNames,MarkerCoords);                
                NFrames = size(MarkerCoords,1);

            elseif strcmpi(R.Motion.File(end-2:end),'mot')||strcmpi(R.Motion.File(end-2:end),'v12')
                str = sprintf(['  Reading the motion file:',R.Motion.File,'...']);
                if(R.Settings.Display == 1 || R.Settings.Display == 2)
                    disp(str);
                end
                fprintf(R.ExperLogFileId, '%s\n', str);
                [Frequency, MarkerNames, MarkerCoords] = R.readMOT();
                R.RawMarkerCoords = MarkerCoords;
                R.RawMarkerNames = MarkerNames;
                MarkerCoords = filtTrajNaN(MarkerCoords, {'butter';R.Settings.Smoothing.CutFreq; Frequency});
                NFrames = size(MarkerCoords,1);

            elseif strcmpi(R.Motion.File(end-2:end),'c3d') || strcmpi(R.Motion.File(end-2:end),'mat')
                str = sprintf(['  Reading the motion file:',R.Motion.File,'...']);
                if(R.Settings.Display == 1 || R.Settings.Display == 2)
                    disp(str);
                end
                fprintf(R.ExperLogFileId, '%s\n', str);
                
                
                [~, Extension] = getFilenameAndExt(R.Motion.File);
                
                if strcmpi(Extension,'c3d')
                    % Option 1: read marker trajectories file from C3D FILE
                    [Frequency, MarkerNames, MarkerCoords] = R.readC3D();
                
                elseif strcmpi(Extension,'mat')
                    % Option 2: read marker trajectories file from MAT FILE
                    [Frequency, MarkerNames, MarkerCoords] = R.readMAT();
                end
                
                
                R.RawMarkerCoords = MarkerCoords;
                R.RawMarkerNames = MarkerNames;
                % smoothing
                if strcmpi(R.Settings.Smoothing.Method, 'butter')
                    CutFreq = R.Settings.Smoothing.CutFreq;
                    str = ['  Filtering marker trajectories... Method: butter, Cut-off freq: ',num2str(CutFreq)];
                    if(R.Settings.Display == 1 || R.Settings.Display == 2), disp(str); end
                    fprintf(R.ExperLogFileId, '%s\n', str);                    
                    MarkerCoords = filtTrajNaN(MarkerCoords, {'butter'; CutFreq; Frequency});                    
                    
                elseif isempty(R.Settings.Smoothing.Method)
                    str = '  No filtering applied to marker trajectories.';
                    if(R.Settings.Display == 1 || R.Settings.Display == 2), disp(str); end
                    fprintf(R.ExperLogFileId, '%s\n', str);
                    
                else
                    str = ['  Filtering method "',R.Settings.Smoothing.Method,'" is unknown.'];
                    fprintf(R.ExperLogFileId, '%s\n', str);
                    error(str);
                    
                end
                NFrames = size(MarkerCoords,1);
                
            end
        end
        function DataPressure = readPressureFile(R)
            % Path and File
            InputPath = R.PressurePath;
            InputFile = R.Motion.File; % The name of Pressure maps files must be the same as motion file
            % Open the input file for reading
            FidInput = fopen([InputPath,InputFile]);
            if FidInput == -1
                str = (['  WARNING:The pressure file ',InputFile,' can not be opened in path (',InputPath,'). ']);
                disp(str);
                fprintf(R.ExperLogFileId, '%s\n', str);
                DataPressure = [];
%                 error(['The pressure file ',InputFile,' can not be opened in path (',InputPath,'). ']);
                return;
            end
            % Inizialize variables
            DataPressure = [];
            Fr=0;   % frame counter
%             i = 1;  % This counter represents the read row of the file
            r = 0;  % Counter of row for each sensor
            c = 0;  % Counter of column for each sensor
            TypeSensor = 0; % 1-Seat sensor, 2-Back rest sensor
            % While the reading file hasn't reached the end of file keep calculating
            while feof(FidInput) == 0
                % Get the whole line
                Line = fgets(FidInput);
                k = strfind(Line, ',');
                LineInfo = Line(1:k-1);
                if strcmpi(LineInfo,'Frame:')
                    Fr = Fr +1;
                elseif strcmpi(LineInfo,'Units:')
                    DataPressure.Units = deblank(Line(k+1:end));
                elseif strcmpi(LineInfo,'Time:')
                    DataPressure.Back(Fr).Time = deblank(Line(k+1:end));
                    DataPressure.Seat(Fr).Time = deblank(Line(k+1:end));
                elseif strcmpi(LineInfo,'Sensor:')
                    TypeSensor = str2num(deblank(Line(k+1:end)));
                    r = 0;
                    c = 0;
%                 elseif strcmpi(LineInfo,'COP Row:')
%                     if TypeSensor == 2
%                         DataPressure.Back(Fr).CoPRow = str2num(deblank(Line(k+1:end)));
%                     end
%                 elseif strcmpi(LineInfo,'COP Column:')
%                     if TypeSensor == 2
%                         DataPressure.Back(Fr).CoPCol = str2num(deblank(Line(k+1:end)));
%                     end
                elseif strcmpi(LineInfo,'Sensels:')
                    r = r+1;
                    c = c+1;
                    NRow = size(k,2);
                    for i=1:NRow-1
                        if TypeSensor == 1
                            DataPressure.Seat(Fr).Values(r,c) = str2num(Line(k(i):k(i+1)));
                        elseif TypeSensor == 2
                            DataPressure.Back(Fr).Values(r,c) = str2num(Line(k(i):k(i+1)));
                        else 
                            error('In pressure map only two type of sensor can be defined')
                        end
                        c = c+1;
                    end
                    if TypeSensor == 1
                        DataPressure.Seat(Fr).Values(r,c) = str2num(Line(k(NRow):end));
                    elseif TypeSensor == 2
                        DataPressure.Back(Fr).Values(r,c) = str2num(Line(k(NRow):end));
                    end
                    c = 0;
                end
            end
%             NCol = size(DataPressure.Seat(1).Values,2);
%             for i =1:Fr
%                 CBreak = 0;
%                 for j=1:NRow
%                     if CBreak == 1
%                         break;
%                     end
%                     SumRow = sum(DataPressure.Seat(i).Values(j,:));
%                     if SumRow > 0
%                         for k=1:NCol-1
%                             SumColk = DataPressure.Seat(i).Values(j,k) + DataPressure.Seat(i).Values(j,k+1);
%                             if SumColk>0 && SumColk == DataPressure.Seat(i).Values(j,k) && sum(DataPressure.Seat(i).Values(j,k+1:NCol)) > 0 && sum(DataPressure.Seat(i).Values(j:NRow,k+1)) == 0
%                                 DataPressure.Seat(i).ZeroRow = j;
%                                 DataPressure.Seat(i).ZeroCol = k+1;
%                                 CBreak =1;
%                                 break;
%                             end
%                         end
%                     end
%                 end
%             end
            
        end 
        function resultsIK(R,ResultsPath,GraphicPath)
            % ResultName
            if isempty(R.Motion)
                ResultName = [R.ExperimentName,'_',R.Subject.ModelName];
            else
                ResultName = [R.ExperimentName,'_',R.Subject.ModelName,'_',R.Motion.File(1:(end-4))];
            end
            % Result for the init posture
            if R.Settings.Results.InitPosture == 1
                str = '    Writing initial aproximation in COMPAMM format:';
                fprintf(R.ExperLogFileId, '%s\n', str);
                if(R.Settings.Display == 1 || R.Settings.Display == 2)
                    disp(str);
                end
                %Results_to = RESULTS(R.q0',R.z,R.g_t,R.qdot_t,R.qdot2_t, R.Subject,ResultsPath,R.Deltat,R.Markers,R.ExperLogFileId,R.Settings,R.RawMarkerCoords);
                Results_to = RESULTS(R,R.q0',ResultsPath);
                Results_to.writeFilePlayback([ResultName,'_t0'],GraphicPath,R.MarkerSensor);
                Results_to.writeFileMatrices([ResultName,'_t0'],R.MarkerSensor);
                %                 Results_to.writeFileMatrices([ResultName,'_t0']);
            end
            % Kinematic results of the reconstruction
            if strcmpi(R.ResultType,'IK') 
                Results = RESULTS(R,R.q_t,ResultsPath);
                Results.getSensorData();
                                
                if R.Settings.Results.CompPlayback == 1
                    str = '    Writing reconstructed motion in COMPAMM format:';
                    fprintf(R.ExperLogFileId, '%s\n', str);
                    if(R.Settings.Display == 1 || R.Settings.Display == 2)
                        disp(str);
                    end
                    Results.writeFilePlayback(ResultName,GraphicPath,R.MarkerSensor);
                    Results.writeFileMatrices(ResultName,R.MarkerSensor);
                    %Results.writeFileMatrices(ResultName);
                end
                if (R.Settings.Results.Position == 1 || R.Settings.Results.Velocity == 1 || R.Settings.Results.Acceleration == 1 || ...
                        R.Settings.Results.MarkerError == 1 || R.Settings.Results.RawMarkerTrajectory == 1 || R.Settings.Results.Sensor == 1 || ...
                        R.Settings.Results.CentreOfMass == 1 || R.Settings.Results.AcondMarkerTrajectory == 1)
                    Results.result_kin(ResultName);
                end
                if strcmpi(R.Subject.ModelType,'Ramsis')&& R.Settings.Results.Ramsis == 1
                    Results.result_Ramsis(ResultName);
                elseif strcmpi(R.Subject.ModelType,'PAM')&& R.Settings.Results.PAM == 1
                    Results.result_PAM(ResultName);
                elseif (strcmpi(R.Subject.ModelType,'Ramsis')&& R.Settings.Results.Ramsis == 0) || ...
                       (strcmpi(R.Subject.ModelType,'PAM')&& R.Settings.Results.PAM == 0) 
                    str = ['  WARNING: Model type is ',R.Subject.ModelType,' BUT IK results in ',R.Subject.ModelType,' format not asked for'];
                    fprintf(R.ExperLogFileId, '%s\n', str);
                    disp(str);
                end
                %save([ResultsPath,ResultName,'_IK'],'R');
            end
            
        end
        function resultsID(R,ResultsPath)
            % ResultName
            if isempty(R.Motion)
                ResultName = [R.ExperimentName,'_',R.Subject.ModelName];
            else
                ResultName = [R.ExperimentName,'_',R.Subject.ModelName,'_',R.Motion.File(1:(end-4))];
            end
            % Dynamic results of the reconstruction
            if strcmpi(R.ResultType,'ID')
                Results = RESULTS(R,R.q_t,ResultsPath);
                if R.Settings.Results.JointEffort == 1
                    Results.resultID(ResultName);
                end
                if strcmpi(R.Subject.ModelType,'Ramsis')&& R.Settings.Results.Ramsis == 1
                    Results.getSensorData();
                    Results.result_RamsisDyn(ResultName);
                elseif (strcmpi(R.Subject.ModelType,'Ramsis')&& R.Settings.Results.Ramsis == 0)
                    str = ['  WARNING: Model type is ',R.Subject.ModelType,' BUT ID results in ',R.Subject.ModelType,' format not asked for'];
                    fprintf(R.ExperLogFileId, '%s\n', str);
                    disp(str);
                end
                if strcmp(R.Settings.ExType,'CP')
                    Results.resultsForIDOptimization(ResultName);
                end
            end
        end
        function solverIK(R)
            time = tic;
            str = sprintf(['\n  -------------------------------------------------------------------------\n', ...
                           '   Inverse Kinematic Motion Reconstruction (DHIK) ', ...
                           '\n  -------------------------------------------------------------------------']);
            if(R.Settings.Display == 1 || R.Settings.Display == 2)
                disp(str);
            end
            fprintf(R.ExperLogFileId, '%s\n', str);
            Solver = SOLVER(R,R.ExperLogFileId,R.Settings);
            Solver.(R.SolverType);
            [H, MI, S] = second2HMS(toc(time));
            str =   ['   IK Solver Time: ',num2str(H),' hour(s) ',num2str(MI),' minute(s) and ',num2str(S),' second(s).'];
            if(R.Settings.Display == 1 || R.Settings.Display == 2)
                disp(str);
            end
            fprintf(R.ExperLogFileId,'%s\n', str);
        end
        function solverSoderkvist(R)
            time = tic;
            str = sprintf(['\n  -------------------------------------------------------------------------\n', ...
                           '   Inverse Kinematic Motion Reconstruction (DHIK) ', ...
                           '\n  -------------------------------------------------------------------------\n', ...
                           '   Solving problem with Soderkvist Method ...']);
            if(R.Settings.Display == 1 || R.Settings.Display == 2)
                disp(str);
            end
            fprintf(R.ExperLogFileId, '%s\n', str);
            NFrames = size(R.g_t,1);
            NSegments = size(R.Subject.Segments,1);
            for Fr =1:NFrames
                for i=1:NSegments
                    if R.Subject.Segments(i).Fixed ~= 1 % If is NOT Ground                                                    
                    %if ~strcmp(R.Subject.Segments(i).Name,'Ground')
                        NMarkers = size(R.Subject.Segments(i).LocalMarkers,1);
                        Glob_Pos_Markers = [];
                        Seg_Pos_Markers = [];
                        for j=1:NMarkers
                            PosInzMarker = R.Subject.Segments(i).LocalMarkers(j).Point.PosInz;
                            Glob_Pos_Marker_i = R.g_t(Fr,PosInzMarker:PosInzMarker+2);
                            % Check if experimental marker have a non NaN value
                            if ~isnan(Glob_Pos_Marker_i(1))
                                Glob_Pos_Markers = [Glob_Pos_Markers,Glob_Pos_Marker_i'];
                                Seg_Pos_Markers = [Seg_Pos_Markers,R.Subject.Segments(i).LocalMarkers(j).LocCoord];
                            end
                        end
                        % Three marker are necesary to calc optimal pose.
                        if size(Glob_Pos_Markers,2) > 2
%                             disp(['Frame:',num2str(Fr),' Segment:',num2str(i)]);
                            [Glob_R_Seg, Glob_Pos_Or] = R.Subject.Segments(i).calOptPose(Glob_Pos_Markers, Seg_Pos_Markers);
                            R.q_t = R.Subject.Segments(i).calcq_tWithOptPose(Glob_R_Seg,Glob_Pos_Or,Fr,R.q_t);
                        else
                            % If there are 2 markers or less, Soderkvist can not calculate the segment pose (R,d)
                            % Then, we set R, d to NaN so that other functions can detect this fact
                            Glob_R_Seg = NaN*eye(3);
                            Glob_Pos_Or = [NaN;NaN;NaN];
                            R.q_t = R.Subject.Segments(i).calcq_tWithOptPose(Glob_R_Seg,Glob_Pos_Or,Fr,R.q_t);
                        end
                    end
                end
            end
            [H, MI, S] = second2HMS(toc(time));
            str =   ['   IK Solver Time: ',num2str(H),' hour(s) ',num2str(MI),' minute(s) and ',num2str(S),' second(s).'];
            if(R.Settings.Display == 1 || R.Settings.Display == 2)
                disp(str);
            end
            fprintf(R.ExperLogFileId,'%s\n', str);
        end
        function storeEffortValues(R,DataCSV,MarkerNames,MarkerCoords,Frequency)
            ForceDom = [];
            
            % EffortFile = get from Settings
            % EffortPath = get from Settings
            
            FileExist = exist([EffortPath,EffortFile], 'file');
            %FileExist = exist([R.ExperimentPath,R.Subject.ModelName,'.forces'],'file');
            if FileExist ~=2
                str = sprintf(['  WARNING: There is not sensor force definition file: ',R.Subject.ModelName,'.forces\n'...
                               '           in path ',getPrintPath(R.ExperimentPath)]);
                disp(str);
                fprintf(R.ExperLogFileId,'%s\n', str);
            else
                % Read the force file in xml format
                ForceDom = xmlread([R.ExperimentPath,R.Subject.ModelName,'.forces']);
            end
            if ~isempty(ForceDom)
                
%                 SensorFile =  get from EffortFile in element SensorFile
                
                SensorDom = xmlread([R.ExperimentPath,'Sensor.sdf']);
                Sensor = SensorDom.getElementsByTagName('SENSOR');
                for i = 1:Sensor.getLength
                    SensorListItem = Sensor.item(i-1);
                    SensorName = char(SensorListItem.getAttribute('Name'));
                    Sen.(SensorName)= [];
                    if strcmp(SensorName,'Sensor1')
                        Sen.Sensor1.Marker_X.Label        = char(SensorListItem.getElementsByTagName('MARKER_X_DIR_POS').item(0).getAttribute('Label'));
                        Sen.Sensor1.Marker_X.LocCoord     = str2num(char(SensorListItem.getElementsByTagName('MARKER_X_DIR_POS').item(0).getAttribute('LocCoord')));
                        Sen.Sensor1.Marker_Y_Neg.Label    = char(SensorListItem.getElementsByTagName('MARKER_Y_DIR_NEG').item(0).getAttribute('Label'));
                        Sen.Sensor1.Marker_Y_Neg.LocCoord = str2num(char(SensorListItem.getElementsByTagName('MARKER_Y_DIR_NEG').item(0).getAttribute('LocCoord')));
                        Sen.Sensor1.Marker_Y_Pos.Label    = char(SensorListItem.getElementsByTagName('MARKER_Y_DIR_POS').item(0).getAttribute('Label'));
                        Sen.Sensor1.Marker_Y_Pos.LocCoord = str2num(char(SensorListItem.getElementsByTagName('MARKER_Y_DIR_POS').item(0).getAttribute('LocCoord')));
                        Sen.Sensor1.Marker_Z.Label        = char(SensorListItem.getElementsByTagName('MARKER_Z_DIR_POS').item(0).getAttribute('Label'));
                        Sen.Sensor1.Marker_Z.LocCoord     = str2num(char(SensorListItem.getElementsByTagName('MARKER_Z_DIR_POS').item(0).getAttribute('LocCoord')));
                        Sen.Sensor1.Marker_Rot1.Label     = char(SensorListItem.getElementsByTagName('MARKER_DIR_ROT1').item(0).getAttribute('Label'));
                        Sen.Sensor1.Marker_Rot1.LocCoord  = str2num(char(SensorListItem.getElementsByTagName('MARKER_DIR_ROT1').item(0).getAttribute('LocCoord')));
                        Sen.Sensor1.Marker_Rot2.Label     = char(SensorListItem.getElementsByTagName('MARKER_DIR_ROT2').item(0).getAttribute('Label'));
                        Sen.Sensor1.Marker_Rot2.LocCoord  = str2num(char(SensorListItem.getElementsByTagName('MARKER_DIR_ROT2').item(0).getAttribute('LocCoord')));
                        Sen.Sensor1.Marker_X_Axis_Neg.Label = char(SensorListItem.getElementsByTagName('MARKER_X_AXIS_NEG').item(0).getAttribute('Label'));
                        Sen.Sensor1.Marker_Y_Axis_Neg.Label = char(SensorListItem.getElementsByTagName('MARKER_Y_AXIS_NEG').item(0).getAttribute('Label'));
                        Sen.Sensor1.DifMarkers.LocCoord   =  str2num(char(SensorListItem.getElementsByTagName('DIF_AXIS').item(0).getAttribute('LocCoord')));
                        Sen.Sensor1.Label_Fx = char(SensorListItem.getAttribute('LABEL_FORCE_X'));
                        Sen.Sensor1.Label_Fy = char(SensorListItem.getAttribute('LABEL_FORCE_Y'));
                        Sen.Sensor1.Label_Fz = char(SensorListItem.getAttribute('LABEL_FORCE_Z'));
                    end
                end
                
                % read EffortFile
                AppliedForce = ForceDom.getElementsByTagName('AppliedForce');
                for i = 1:AppliedForce.getLength
                    ApliedForceListItem = AppliedForce.item(i-1);
                    SegmentName = char(ApliedForceListItem.getAttribute('Segment'));
                    SensorName  = char(ApliedForceListItem.getAttribute('Sensor'));
                    R.calcSegmentForces(DataCSV,MarkerNames,MarkerCoords,Sen.(SensorName),SensorName,SegmentName,Frequency);
                    
                end
            end
        end
        function storeForceValues(R,DataCSV,MarkerNames,MarkerCoords,Frequency)
            ForceDom = [];
            FileExist = exist([R.ExperimentPath,R.Subject.ModelName,'.forces'],'file');
            if FileExist ~=2
                str = sprintf(['  WARNING: There is not sensor force definition file: ',R.Subject.ModelName,'.forces\n'...
                               '           in path ',getPrintPath(R.ExperimentPath)]);
                disp(str);
                fprintf(R.ExperLogFileId,'%s\n', str);
            else
                % Read the force file in xml format
                ForceDom = xmlread([R.ExperimentPath,R.Subject.ModelName,'.forces']);
            end
            if ~isempty(ForceDom)
                SensorDom = xmlread([R.ExperimentPath,'Sensor.sdf']);
                Sensor = SensorDom.getElementsByTagName('SENSOR');
                for i = 1:Sensor.getLength
                    SensorListItem = Sensor.item(i-1);
                    SensorName = char(SensorListItem.getAttribute('Name'));
                    Sen.(SensorName)= [];
                    if strcmp(SensorName,'Sensor1')
                        Sen.Sensor1.Marker_X.Label        = char(SensorListItem.getElementsByTagName('MARKER_X_DIR_POS').item(0).getAttribute('Label'));
                        Sen.Sensor1.Marker_X.LocCoord     = str2num(char(SensorListItem.getElementsByTagName('MARKER_X_DIR_POS').item(0).getAttribute('LocCoord')));
                        Sen.Sensor1.Marker_Y_Neg.Label    = char(SensorListItem.getElementsByTagName('MARKER_Y_DIR_NEG').item(0).getAttribute('Label'));
                        Sen.Sensor1.Marker_Y_Neg.LocCoord = str2num(char(SensorListItem.getElementsByTagName('MARKER_Y_DIR_NEG').item(0).getAttribute('LocCoord')));
                        Sen.Sensor1.Marker_Y_Pos.Label    = char(SensorListItem.getElementsByTagName('MARKER_Y_DIR_POS').item(0).getAttribute('Label'));
                        Sen.Sensor1.Marker_Y_Pos.LocCoord = str2num(char(SensorListItem.getElementsByTagName('MARKER_Y_DIR_POS').item(0).getAttribute('LocCoord')));
                        Sen.Sensor1.Marker_Z.Label        = char(SensorListItem.getElementsByTagName('MARKER_Z_DIR_POS').item(0).getAttribute('Label'));
                        Sen.Sensor1.Marker_Z.LocCoord     = str2num(char(SensorListItem.getElementsByTagName('MARKER_Z_DIR_POS').item(0).getAttribute('LocCoord')));
                        Sen.Sensor1.Marker_Rot1.Label     = char(SensorListItem.getElementsByTagName('MARKER_DIR_ROT1').item(0).getAttribute('Label'));
                        Sen.Sensor1.Marker_Rot1.LocCoord  = str2num(char(SensorListItem.getElementsByTagName('MARKER_DIR_ROT1').item(0).getAttribute('LocCoord')));
                        Sen.Sensor1.Marker_Rot2.Label     = char(SensorListItem.getElementsByTagName('MARKER_DIR_ROT2').item(0).getAttribute('Label'));
                        Sen.Sensor1.Marker_Rot2.LocCoord  = str2num(char(SensorListItem.getElementsByTagName('MARKER_DIR_ROT2').item(0).getAttribute('LocCoord')));
                        Sen.Sensor1.Marker_X_Axis_Neg.Label = char(SensorListItem.getElementsByTagName('MARKER_X_AXIS_NEG').item(0).getAttribute('Label'));
                        Sen.Sensor1.Marker_Y_Axis_Neg.Label = char(SensorListItem.getElementsByTagName('MARKER_Y_AXIS_NEG').item(0).getAttribute('Label'));
                        Sen.Sensor1.DifMarkers.LocCoord   =  str2num(char(SensorListItem.getElementsByTagName('DIF_AXIS').item(0).getAttribute('LocCoord')));
                        Sen.Sensor1.Label_Fx = char(SensorListItem.getAttribute('LABEL_FORCE_X'));
                        Sen.Sensor1.Label_Fy = char(SensorListItem.getAttribute('LABEL_FORCE_Y'));
                        Sen.Sensor1.Label_Fz = char(SensorListItem.getAttribute('LABEL_FORCE_Z'));
                    end
                end
                AppliedForce = ForceDom.getElementsByTagName('AppliedForce');
                for i = 1:AppliedForce.getLength
                    ApliedForceListItem = AppliedForce.item(i-1);
                    SegmentName = char(ApliedForceListItem.getAttribute('Segment'));
                    SensorName  = char(ApliedForceListItem.getAttribute('Sensor'));
                    R.calcSegmentForces(DataCSV,MarkerNames,MarkerCoords,Sen.(SensorName),SensorName,SegmentName,Frequency);
                    
                end
            end
        end
        function storePressureData(R,DataPressure,MarkerNames,MarkerCoords,Frequency)
            ForceDom = [];
            FileExist = exist([R.ExperimentPath,R.Subject.ModelName,'.effort'],'file');
            if FileExist ~=2
                str = sprintf(['  WARNING: There is not sensor force definition file: ',R.Subject.ModelName,'.effort\n'...
                               '           in path ',getPrintPath(R.ExperimentPath)]);
                disp(str);
                fprintf(R.ExperLogFileId,'%s\n', str);
            else
                % Read the force file in xml format
                ForceDom = xmlread([R.ExperimentPath,R.Subject.ModelName,'.effort']);
                AppliedForce = ForceDom.getElementsByTagName('AppliedEffort');
                % For each AppliedForce we need calc this force in the segment
                for i = 1:AppliedForce.getLength
                    ApliedForceListItem = AppliedForce.item(i-1);
                    SegmentName = char(ApliedForceListItem.getAttribute('Segment'));
                    SensorName  = char(ApliedForceListItem.getAttribute('SensorName'));
                    SensorFile = char(ApliedForceListItem.getAttribute('SensorFile'));
                    SensorDom = xmlread([R.ExperimentPath,SensorFile]);
                    Sensor = SensorDom.getElementsByTagName('SENSOR');
                    for j = 1:Sensor.getLength
                        SensorListItem = Sensor.item(j-1);
                        SensorNamej = char(SensorListItem.getAttribute('Name'));
                        % Sensor Information
                        if strcmp(SensorName,SensorNamej)
                            AplPointItem = SensorListItem.getElementsByTagName('APLICATION_POINT').item(0);
                            Sen.(SensorName).AplPoint.Name  = char(AplPointItem.getAttribute('Name'));
                            Sen.(SensorName).AplPoint.Coord = str2num(char(AplPointItem.getAttribute('GlobCoord')));
                            DirList = SensorListItem.getElementsByTagName('DIR');
                            for k=1:DirList.getLength
                                DirListItem  = DirList.item(k-1);
                                DirName = char(DirListItem.getAttribute('Name'));
                                P1Name = char(DirListItem.getAttribute('Point1')); if ~isempty(P1Name) Sen.(SensorName).(DirName).P1Name = P1Name; end
                                P2Name = char(DirListItem.getAttribute('Point2')); if ~isempty(P2Name) Sen.(SensorName).(DirName).P2Name = P2Name; end
                                P3Name = char(DirListItem.getAttribute('Point3')); if ~isempty(P3Name) Sen.(SensorName).(DirName).P3Name = P3Name; end
                                P4Name = char(DirListItem.getAttribute('Point4')); if ~isempty(P4Name) Sen.(SensorName).(DirName).P4Name = P4Name; end
                                AngleY = str2num(char(DirListItem.getAttribute('AngleY'))); if ~isempty(AngleY) Sen.(SensorName).(DirName).AngleY = AngleY; end
                            end
                        end
                    end
                    R.calcSegmentForces(DataPressure,MarkerNames,MarkerCoords,Sen.(SensorName),SensorName,SegmentName,Frequency);
                end
            end
            
        end
        function wFileFillPhi(R)
            % define function name for phi
            Filename = [R.Subject.ModelName,'_',R.ExperimentName,'_fillphi'];
            FilenameExt = [Filename,'.m'];
            R.FilePhiName = Filename;
            
            % sizes
            NSeg    = size(R.Subject.Segments,1);
            nCoords = size(R.Subject.q,1);
            NEqs    = size(R.Subject.Phi,1);
            Time = clock;
            
            % open file
            fid = fopen([R.ModelEqPath,FilenameExt],'w');
            fprintf(fid,['function Phi = ',Filename,'(q,Par)\n\n']);
            fprintf(fid, '%% Author  : CEIT\n');
            fprintf(fid,['%% Date  : ',date,'\n']);
            fprintf(fid,['%% Time  : ',num2str(Time(4)),':',num2str(Time(5)),'\n']);
            fprintf(fid,['%% Model : ',R.Subject.ModelName,'\n']);
            fprintf(fid, '%% Version: 2.0 CEIT\n\n');
            
            % write local coordinates of points and markers
            NPar = 0;
            for i=1:NSeg
                if R.Subject.Segments(i).Fixed ~= 1 % If is NOT Ground                                                    
                %if ~strcmpi(R.Subject.Segments(i).Name,'Ground')
                    NPoints = size(R.Subject.Segments(i).LocalPoints,1);
                    NMarkers = size(R.Subject.Segments(i).LocalMarkers,1);
                    NVectors = size(R.Subject.Segments(i).LocalVectors,1);
                    %                 for j=1:NPoints
                    %                     fprintf(fid,[R.Subject.Segments(i).LocalPoints(j).CoordName{1}, ...
                    %                         ' = Subject.Segments(',num2str(i),').LocalPoints(',num2str(j),').LocCoord(1);\n']);
                    %                     fprintf(fid,[R.Subject.Segments(i).LocalPoints(j).CoordName{2},...
                    %                         ' = Subject.Segments(',num2str(i),').LocalPoints(',num2str(j),').LocCoord(2);\n']);
                    %                     fprintf(fid,[R.Subject.Segments(i).LocalPoints(j).CoordName{3},...
                    %                         ' = Subject.Segments(',num2str(i),').LocalPoints(',num2str(j),').LocCoord(3);\n']);
                    %                 end
                    %                 for j=1:NMarkers
                    %                     fprintf(fid,[R.Subject.Segments(i).LocalMarkers(j).CoordName{1}, ...
                    %                         ' = Subject.Segments(',num2str(i),').LocalMarkers(',num2str(j),').LocCoord(1);\n']);
                    %                     fprintf(fid,[R.Subject.Segments(i).LocalMarkers(j).CoordName{2},...
                    %                         ' = Subject.Segments(',num2str(i),').LocalMarkers(',num2str(j),').LocCoord(2);\n']);
                    %                     fprintf(fid,[R.Subject.Segments(i).LocalMarkers(j).CoordName{3},...
                    %                         ' = Subject.Segments(',num2str(i),').LocalMarkers(',num2str(j),').LocCoord(3);\n']);
                    %                 end
                    for j=1:NPoints
                        fprintf(fid,[R.Subject.Segments(i).LocalPoints(j).CoordName{1}, ...
                            ' = Par(',num2str(3*j-2+NPar),');\n']);
                        fprintf(fid,[R.Subject.Segments(i).LocalPoints(j).CoordName{2},...
                            ' = Par(',num2str(3*j-1+NPar),');\n']);
                        fprintf(fid,[R.Subject.Segments(i).LocalPoints(j).CoordName{3},...
                            ' = Par(',num2str(3*j+NPar),');\n']);
                    end
                    for j=1:NMarkers
                        fprintf(fid,[R.Subject.Segments(i).LocalMarkers(j).CoordName{1}, ...
                            ' = Par(',num2str(3*j-2+3*NPoints+NPar),');\n']);
                        fprintf(fid,[R.Subject.Segments(i).LocalMarkers(j).CoordName{2},...
                            ' = Par(',num2str(3*j-1+3*NPoints+NPar),');\n']);
                        fprintf(fid,[R.Subject.Segments(i).LocalMarkers(j).CoordName{3},...
                            ' = Par(',num2str(3*j+3*NPoints+NPar),');\n']);
                    end
                    for j=1:NVectors
                        fprintf(fid,[R.Subject.Segments(i).LocalVectors(j).CoordName{1}, ...
                            ' = Par(',num2str(3*j-2+3*NPoints+3*NMarkers+NPar),');\n']);
                        fprintf(fid,[R.Subject.Segments(i).LocalVectors(j).CoordName{2}, ...
                            ' = Par(',num2str(3*j-1+3*NPoints+3*NMarkers+NPar),');\n']);
                        fprintf(fid,[R.Subject.Segments(i).LocalVectors(j).CoordName{3}, ...
                            ' = Par(',num2str(3*j+3*NPoints+3*NMarkers+NPar),');\n']);
                    end
                    NPar = 3*NPoints + 3*NMarkers + 3*NVectors+NPar;    
                end
            end
            fprintf(fid,'\n');


            % write generalized coordinates
            fprintf(fid,'\n');
            for i=1:nCoords
                fprintf(fid,[char(R.Subject.q(i)),' = q(', num2str(i),');\n'] );
            end
            fprintf(fid,'\n');
            
            % write constant variables of the FIXED Body
            % get FIXED body (Ground) data
            for i=1:NSeg
                if R.Subject.Segments(i).Fixed == 1 % If is Fixed Body (Ground)
                    FixedBodyIndex = i;
                end
            end            
            
            NPoints = size(R.Subject.Segments(FixedBodyIndex).LocalPoints,1);
            NVectors = size(R.Subject.Segments(FixedBodyIndex).LocalVectors,1);
       
            for i = 1:NPoints
                fprintf(fid,[R.Subject.Segments(FixedBodyIndex).LocalPoints(i).Point.CoordName{1},' = ', num2str(R.Subject.Segments(FixedBodyIndex).LocalPoints(i).LocCoord(1)),';\n']);
                fprintf(fid,[R.Subject.Segments(FixedBodyIndex).LocalPoints(i).Point.CoordName{2},' = ', num2str(R.Subject.Segments(FixedBodyIndex).LocalPoints(i).LocCoord(2)),';\n']);
                fprintf(fid,[R.Subject.Segments(FixedBodyIndex).LocalPoints(i).Point.CoordName{3},' = ', num2str(R.Subject.Segments(FixedBodyIndex).LocalPoints(i).LocCoord(3)),';\n']);
            end
            
            for i = 1:NVectors
                fprintf(fid,[R.Subject.Segments(FixedBodyIndex).LocalVectors(i).Vector.CoordName{1},' = ', num2str(R.Subject.Segments(FixedBodyIndex).LocalVectors(i).LocCoord(1)),';\n']);
                fprintf(fid,[R.Subject.Segments(FixedBodyIndex).LocalVectors(i).Vector.CoordName{2},' = ', num2str(R.Subject.Segments(FixedBodyIndex).LocalVectors(i).LocCoord(2)),';\n']);
                fprintf(fid,[R.Subject.Segments(FixedBodyIndex).LocalVectors(i).Vector.CoordName{3},' = ', num2str(R.Subject.Segments(FixedBodyIndex).LocalVectors(i).LocCoord(3)),';\n']);                         
            end
            fprintf(fid,'\n');
            
            % write equations
            fprintf(fid,['Phi = zeros(',num2str(NEqs),',1);\n']);
            for i=1:NEqs
                fprintf(fid,['Phi(',num2str(i),') = ',R.Subject.Phi{i,1},';\n']);
            end
            fprintf(fid,'\n');

            % close file
            fclose(fid);

        end
        function wFileFillPhiqSparse(R)
            % define function name for phi
            Filename = [R.Subject.ModelName,'_',R.ExperimentName,'_fillphiq'];
            FilenameExt = [Filename,'.m'];
            R.FilePhiqName = Filename;
            % sizes
            NSeg    = size(R.Subject.Segments,1);
            nCoords = size(R.Subject.q,1);
            NEqs    = size(R.Subject.Phi,1);
            Time    = clock;
            % open file
            fid = fopen([R.ModelEqPath,FilenameExt],'w');
            % write function heading 
            fprintf(fid,['function Phiq = ',Filename,'(q,Par)\n\n']);
            fprintf(fid, '%% Author  : CEIT\n');
            fprintf(fid,['%% Date  : ',date,'\n']);
            fprintf(fid,['%% Time  : ',num2str(Time(4)),':',num2str(Time(5)),'\n']);
            fprintf(fid,['%% Model : ',R.Subject.ModelName,'\n']);
            fprintf(fid, '%% Version: 2.0 CEIT\n\n');
            
            % write local coordinates
            NPar = 0;
            for i=1:NSeg
                if R.Subject.Segments(i).Fixed ~= 1 % If is NOT Ground                                                    
                %if ~strcmpi(R.Subject.Segments(i).Name,'Ground')
                    NPoints = size(R.Subject.Segments(i).LocalPoints,1);
                    NMarkers = size(R.Subject.Segments(i).LocalMarkers,1);
                    NVectors = size(R.Subject.Segments(i).LocalVectors,1);
                    %                 for j=1:NPoints
                    %                     fprintf(fid,[R.Subject.Segments(i).LocalPoints(j).CoordName{1}, ...
                    %                         ' = Subject.Segments(',num2str(i),').LocalPoints(',num2str(j),').LocCoord(1);\n']);
                    %                     fprintf(fid,[R.Subject.Segments(i).LocalPoints(j).CoordName{2},...
                    %                         ' = Subject.Segments(',num2str(i),').LocalPoints(',num2str(j),').LocCoord(2);\n']);
                    %                     fprintf(fid,[R.Subject.Segments(i).LocalPoints(j).CoordName{3},...
                    %                         ' = Subject.Segments(',num2str(i),').LocalPoints(',num2str(j),').LocCoord(3);\n']);
                    %                 end
                    %                 for j=1:NMarkers
                    %                     fprintf(fid,[R.Subject.Segments(i).LocalMarkers(j).CoordName{1}, ...
                    %                         ' = Subject.Segments(',num2str(i),').LocalMarkers(',num2str(j),').LocCoord(1);\n']);
                    %                     fprintf(fid,[R.Subject.Segments(i).LocalMarkers(j).CoordName{2},...
                    %                         ' = Subject.Segments(',num2str(i),').LocalMarkers(',num2str(j),').LocCoord(2);\n']);
                    %                     fprintf(fid,[R.Subject.Segments(i).LocalMarkers(j).CoordName{3},...
                    %                         ' = Subject.Segments(',num2str(i),').LocalMarkers(',num2str(j),').LocCoord(3);\n']);
                    %                 end
                    for j=1:NPoints
                        fprintf(fid,[R.Subject.Segments(i).LocalPoints(j).CoordName{1}, ...
                            ' = Par(',num2str(3*j-2+NPar),');\n']);
                        fprintf(fid,[R.Subject.Segments(i).LocalPoints(j).CoordName{2},...
                            ' = Par(',num2str(3*j-1+NPar),');\n']);
                        fprintf(fid,[R.Subject.Segments(i).LocalPoints(j).CoordName{3},...
                            ' = Par(',num2str(3*j+NPar),');\n']);
                    end
                    for j=1:NMarkers
                        fprintf(fid,[R.Subject.Segments(i).LocalMarkers(j).CoordName{1}, ...
                            ' = Par(',num2str(3*j-2+3*NPoints+NPar),');\n']);
                        fprintf(fid,[R.Subject.Segments(i).LocalMarkers(j).CoordName{2},...
                            ' = Par(',num2str(3*j-1+3*NPoints+NPar),');\n']);
                        fprintf(fid,[R.Subject.Segments(i).LocalMarkers(j).CoordName{3},...
                            ' = Par(',num2str(3*j+3*NPoints+NPar),');\n']);
                    end
                    for j=1:NVectors
                        fprintf(fid,[R.Subject.Segments(i).LocalVectors(j).CoordName{1}, ...
                            ' = Par(',num2str(3*j-2+3*NPoints+3*NMarkers+NPar),');\n']);
                        fprintf(fid,[R.Subject.Segments(i).LocalVectors(j).CoordName{2}, ...
                            ' = Par(',num2str(3*j-1+3*NPoints+3*NMarkers+NPar),');\n']);
                        fprintf(fid,[R.Subject.Segments(i).LocalVectors(j).CoordName{3}, ...
                            ' = Par(',num2str(3*j+3*NPoints+3*NMarkers+NPar),');\n']);
                    end
                    NPar = 3*NPoints + 3*NMarkers + 3*NVectors+NPar;   
                end
            end
            fprintf(fid,'\n');
            
            % write generalized coordinates
            for i=1:nCoords
                fprintf(fid,[char(R.Subject.q(i)),' = q(', num2str(i),');\n'] );
            end
            
            % write constant variables of the FIXED Body
            % get FIXED body (Ground) data
            for i=1:NSeg
                if R.Subject.Segments(i).Fixed == 1 % If is Fixed Body (Ground)
                    FixedBodyIndex = i;
                end
            end
            
            
            
            NPoints = size(R.Subject.Segments(FixedBodyIndex).LocalPoints,1);
            NVectors = size(R.Subject.Segments(FixedBodyIndex).LocalVectors,1);
       
            for i = 1:NPoints
                fprintf(fid,[R.Subject.Segments(FixedBodyIndex).LocalPoints(i).Point.CoordName{1},' = ', num2str(R.Subject.Segments(FixedBodyIndex).LocalPoints(i).LocCoord(1)),';\n']);
                fprintf(fid,[R.Subject.Segments(FixedBodyIndex).LocalPoints(i).Point.CoordName{2},' = ', num2str(R.Subject.Segments(FixedBodyIndex).LocalPoints(i).LocCoord(2)),';\n']);
                fprintf(fid,[R.Subject.Segments(FixedBodyIndex).LocalPoints(i).Point.CoordName{3},' = ', num2str(R.Subject.Segments(FixedBodyIndex).LocalPoints(i).LocCoord(3)),';\n']);
            end
            
            for i = 1:NVectors
                fprintf(fid,[R.Subject.Segments(FixedBodyIndex).LocalVectors(i).Vector.CoordName{1},' = ', num2str(R.Subject.Segments(FixedBodyIndex).LocalVectors(i).LocCoord(1)),';\n']);
                fprintf(fid,[R.Subject.Segments(FixedBodyIndex).LocalVectors(i).Vector.CoordName{2},' = ', num2str(R.Subject.Segments(FixedBodyIndex).LocalVectors(i).LocCoord(2)),';\n']);
                fprintf(fid,[R.Subject.Segments(FixedBodyIndex).LocalVectors(i).Vector.CoordName{3},' = ', num2str(R.Subject.Segments(FixedBodyIndex).LocalVectors(i).LocCoord(3)),';\n']);                         
            end
            
            % write weighting factors of driver coordinates exactly satisfied
            % write value of driven coord. in driver constraints exactly satisfied
            % write sparse jacobian matrix (Phiq)
            % write vector of row indexes
            nonZeroElements = length(R.Subject.Phiq.rows);
            fprintf(fid,'i = [ ');
            for i=1:nonZeroElements
                fprintf(fid,'%d ',R.Subject.Phiq.rows(i));
            end
            fprintf(fid,'];\n\n');
            % write vector of col indexes
            fprintf(fid,'j = [ ');
            for i=1:nonZeroElements
                fprintf(fid,'%d ',R.Subject.Phiq.cols(i));
            end
            fprintf(fid,'];\n\n');
            % write non zero matrix elements
            fprintf(fid,'s = [ ');
            for i=1:nonZeroElements
                fprintf(fid,[R.Subject.Phiq.s{i}, '\n']);
            end
            fprintf(fid,'];\n\n');
            %define de sparse Matlab variable
            fprintf(fid,'Phiq = sparse(i,j,s,%d,%d);\n\n',NEqs,nCoords);
            
            
            % close file ===================================================================
            fclose(fid);
        end
        
        
       
    end
    
end
