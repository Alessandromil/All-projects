classdef RESULTS < handle
    %RESULTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        q_t             % Value of eachgeneralized coordinate at each sample time.      double [NFrames x NVars]
        qdot_t          % The first time derivative of the generalized coordinates.     double [NFrames x NVars]
        qdot2_t         % The 2nd time derivative of the generalized coordinates.       double [NFrames x NVars]
        Subject         % Model with parameters                                         SUBJECT
        DataCSV         % Contain all the data of the CSV
                        % DataCSV.MarkerTraj contain MarkerNames,MarkerCoords,Frecuency
                        % DataCSV.Analog contain Frecuency, and structure with Forces, Torques and EMG
                        %                Each structure have units and values. Example: DataCSV.Analog.Fx1
                        % DataCSV.PMap contain Units, Seat, and Back, Seat and Back contain values
        z
        g_t
        ResultsPath     % Results File path                                             char
        Deltat
        Markers
        MarkerNames
        SensorNames
        SensorValues    % The values of the sensors in deg                              double [NFrames x (NSensorsx3)]
        ForcesNames     % The names of the Forces JointName + Fx/Fy/Fz/Mx/My/Mz         cell{6xNJoints x1}
        ExperLogFileId
        Display
        RawMarkerCoords % The raw coordinates of all the markers in the motion file
        RawMarkerNames  % The names of all the markers in motion file
        Settings        % Struct that contains settings for motion reconstruction
                        % Settings.Type. Inverse Kinematics or Inverse Dynamics. Possible options: IK, ID for IK+ID
                        % Settings.Display. Define ammount of feedback to user: 0(minimum) 1(standard) 2(complete)
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
                        % Settings.Results.AcondMarkerTrajectory. Aconditioned (smoothed & gaps filled) marker trajectories in Compamm format (*.amt)
                        % Settings.Results.Sensor. Variables measured by sensors in Compamm format (*.sen)
                        % Settings.Results.JointEffort. Forces & torques at all the joints in the model. Includes reactions & motor efforts
                        % Settings.Interpolation.Method. Method for interpolating gaps in marker trajectories. Options:`'linear'
                        % Settings.Interpolation.NInterFrames. Threshold (in frames) for interpolating missing markers. Options: any positive integer
                        % Settings.Smoothing.Method. Smoothing method for marker trajectories. Options: 'butter' (butterworth filter)
                        % Settings.Smoothing.CutFreq. Cutt-off frequency in Hz for butterworth filter
    end
    
    methods
        function R = RESULTS(Rec,q_t,ResultsPath)
            R.q_t = q_t; % This variable could be q0 or q_t depending of the calling to the constructor
            R.qdot_t = Rec.qdot_t;
            R.qdot2_t = Rec.qdot2_t;
            R.Subject = Rec.Subject;
            R.z = Rec.z;
            R.g_t = Rec.g_t;
            R.ResultsPath = ResultsPath;
            R.Deltat = Rec.Deltat;
            R.Markers = Rec.Markers;
            R.ExperLogFileId = Rec.ExperLogFileId;
            R.Settings = Rec.Settings;
            R.RawMarkerCoords = Rec.RawMarkerCoords;
            R.RawMarkerNames  = Rec.RawMarkerNames;
            R.DataCSV = Rec.DataCSV;
        end
        function DataDist = calcMarkerError(R)
            % sizes
            [NSamples, NCoords] = size(R.q_t);
            NMarkersExp = size(R.MarkerNames,2);
            % get number of experimetally measured skin-markers and their index in z (or g_t)
%             [MarkerExpIndex, NMarkersExp] = getSmarkIndex(R.z, R.MarkerNames);
            PosInq = [];
            PosInz = [];
            for i=1:NMarkersExp
                MarkerIndex = getVecIndex(R.MarkerNames{i},R.Subject.Markers);
                PosInq = [PosInq;R.Subject.Markers(MarkerIndex).PosInq];
                PosInz = [PosInz;R.Subject.Markers(MarkerIndex).PosInz];
%                 PosInz = [PosInz;R.Markers(MarkerIndex).PosInz];%%%%% Revisar bien
            end
            DataDist  = zeros(NSamples, NMarkersExp);
            % Only one for marker must be show the warning
            Warningj = [];
            for i = 1:NSamples
                for j = 1:NMarkersExp
                    ExperMarker = R.g_t(i,PosInz(j):PosInz(j)+2);
                    SegMarker = R.q_t(i,PosInq(j):PosInq(j)+2);
                    DataDist(i,j) = norm(ExperMarker-SegMarker);
%                     if isnan(DataDist(i,j)) %Para matlab no hace falta ke sea cero/ Para compam si
%                         DataDist(i,j) = 0;
%                     end
                    if DataDist(i,j)>0.10 
                        IsInWarningj = 0;
                        for k=1:size(Warningj)
                            if Warningj(k)== j
                                IsInWarningj = 1;
                                break;
                            end
                        end
                        if IsInWarningj == 0
                            str=['  WARNING: Marker ',R.Subject.q{PosInq(j)}(1:end-1),' could be swapped or have a large error.'];
                            disp(str);
                            fprintf(R.ExperLogFileId, '%s\n', str);
                            Warningj = [j;Warningj];
                        end
                    end
                end
            end
        end
        function CheckValue = checkSensorName(R,RamsisName)
            NSensors = size(R.Subject.Sensors,1);
            CheckValue = 0;
            for i=2:NSensors
                if strcmpi(RamsisName,R.Subject.Sensors{i}.Name)
                    CheckValue = 1;
                end
            end
        end
        function getSensorData(R)
            NSensors = size(R.Subject.Sensors,1);
            j =1;
            for i=1:NSensors
                Sensor = R.Subject.Sensors{i}.getSensorData(R.q_t);
                R.SensorNames{j} = Sensor.Name1;
                R.SensorNames{j+1} = Sensor.Name2;
                R.SensorNames{j+2} = Sensor.Name3;
                j = j+3;
                R.SensorValues = [R.SensorValues,Sensor.Val];
            end
        end
        function CoM = getCoM(R)
            NSegments = size(R.Subject.Segments,1);
            NFrames = size(R.q_t,1);
            for i=1:NSegments
                if R.Subject.Segments(i).Fixed ~= 1 % If is not Fixed(Ground)                    
                    CoMj = [];
                    for j=1:NFrames
                        [Glob_R_Seg,Glob_Pos_OrSeg] = R.Subject.Segments(i).getRd(R.q_t(j,:));
                        CoMj = [CoMj;(Glob_R_Seg*R.Subject.Segments(i).CoM.LocCoord + Glob_Pos_OrSeg)'];
                    end
                    CoM(:,(3*(i-1)-2):3*(i-1)) = CoMj;
                end
            end
        end
        function resultID(R,ResultName)
            
            SaveSubject = R.Subject;
            save([R.ResultsPath,ResultName,'_ID'],'SaveSubject');
            % add force to the file .ind
            fid = fopen([R.ResultsPath,ResultName,'.ind'],'a');
            fprintf(fid,[ResultName,'.jft\r\n']);
            fclose(fid);
            NJoints = size(R.Subject.Joints,1);
            ForceValues = [];
            for i=1:NJoints
                F{1,1} = [R.Subject.Joints{i}.Name,'_Fx'];
                F{2,1} = [R.Subject.Joints{i}.Name,'_Fy'];
                F{3,1} = [R.Subject.Joints{i}.Name,'_Fz'];
                F{4,1} = [R.Subject.Joints{i}.Name,'_Mx'];
                F{5,1} = [R.Subject.Joints{i}.Name,'_My'];
                F{6,1} = [R.Subject.Joints{i}.Name,'_Mz'];
                F{7,1} = [R.Subject.Joints{i}.Name,'_GlobFx'];
                F{8,1} = [R.Subject.Joints{i}.Name,'_GlobFy'];
                F{9,1} = [R.Subject.Joints{i}.Name,'_GlobFz'];
                F{10,1} = [R.Subject.Joints{i}.Name,'_GlobMx'];
                F{11,1} = [R.Subject.Joints{i}.Name,'_GlobMy'];
                F{12,1} = [R.Subject.Joints{i}.Name,'_GlobMz'];
                R.ForcesNames = [R.ForcesNames;F];
                Forces  = R.Subject.Joints{i}.F;
                Moments = R.Subject.Joints{i}.M;
                GlobForces = R.Subject.Joints{i}.FGlob;
                GlobMoments = R.Subject.Joints{i}.MGlob;
                ForceValuei = [Forces,Moments,GlobForces,GlobMoments];
                ForceValues = [ForceValues,ForceValuei];
            end
            % Print force file (*.jft)
            str = sprintf(['    Writing output data in COMPAMM format:', ...
                '\n      ',ResultName,'.jft']);
            fprintf(R.ExperLogFileId, '%s\n', str);
            if(R.Settings.Display == 1 || R.Settings.Display == 2)
                disp(str);
            end
            R.printCompFile([ResultName,'.jft'], ForceValues, R.ForcesNames, 'FORCE');
            
            
        end
        function resultsForIDOptimization(R,ResultName)
            Data4IDOptimization.DataCSV = R.DataCSV;
            DataExtEffort = [];
            DataJointsEffort = [];
            for i=1:length(R.Subject.Segments)
                Data = [];
                if ~isempty(R.Subject.Segments(i).F_Ext)
                    Data.Name = R.Subject.Segments(i).Name;
                    Data.F_Ext = R.Subject.Segments(i).F_Ext;
                end
                    
                if ~isempty((R.Subject.Segments(i).M_Ext))
                    Data.M_Ext = R.Subject.Segments(i).M_Ext;
                end
                DataExtEffort = [DataExtEffort;Data];
            end
            for i=1:length(R.Subject.Joints)
                DataJ.Name = R.Subject.Joints{i}.Name;
                DataJ.F = R.Subject.Joints{i}.F;
                DataJ.M = R.Subject.Joints{i}.M;
                DataJ.FGlob = R.Subject.Joints{i}.FGlob;
                DataJ.MGlob = R.Subject.Joints{i}.MGlob;
                DataJointsEffort = [DataJointsEffort;DataJ];
            end
            Data4IDOptimization.DataExtEffort = DataExtEffort;
            Data4IDOptimization.DataJointsEffort = DataJointsEffort;
            save([R.ResultsPath,ResultName,'_4IDOptimization'],'Data4IDOptimization');
        end
        function result_kin(R,ResultName)
            % RESULTS_KIN writes the kinematic results of the model in COMPAMM file
            % format. The files generated can be open with COMPGRAPH. Three or five
            % files are generated depending on the number of input arguments
            %
            %   This function call be used in three different ways. Depending on
            %   the number of input arguments the files generated are different.
            %    - 8 Input arguments: three files are generated with the position,
            %      velocity and acceleration of the the generalized coordinates.
            %      Files .pos .vel and .acc
            %
            %      result_kin(Filename, Path, Deltat, Body, q, q_t, qdot_t, qdot2_t)
            %
            %    - 11 Input arguments: five files are generated, three with the
            %      position, velocity and acceleration of the the generalized coordinates.
            %      (files .pos .vel and .acc), a fourth file (.exp) with the coordinates of
            %      the measured skin-markers and a fifth file (.dis) with the distance between
            %      each measured skin-marker and the skin-marker rigidly attached to each
            %      body. This distance can be considered as a measure of the error.
            %
            %      result_kin(Filename, Path, Deltat, Body, q, q_t, qdot_t, qdot2_t, z, g_t)
            %
            %   Inputs:
            %     + Filename is the name of the the result files - string
            %     + Path is the path for the result files - string
            %     + Deltat is the sample time of the motion capture file - scalar [seconds]
            %     + Body is a cell with the description of each body of the model
            %       See function MKMODEL for details
            %     + q is the symbolic vector of generalized coordinates
            %     + q_t is a double array(nFrames x nVars) with the value of each
            %       generalized coordinate at each sample time. For the sample time
            %       i, q_t(i,:) containts the values of the generalized coordinates.
            %       The order of the variables is the same as in q.
            %     + qdot_t is a double array(nFrames x nVars) with the first time
            %       derivative of the generalized coordinates for each sample time.
            %     + qdot2_t is a double array(nFrames x nVars) with the 2nd time
            %       derivative of the generalized coordinates for each sample time.
            %     + z is a symbolic array (nInputs x 1) that containts the generalized
            %       coordinates of the model chosen as inputs
            %     + g_t is a double array(nFrames x nInputs) that containts the value for
            %       the inputs for each sample time. See INPMODELIK or SIMULAIK for details.
            %     + SmarkNames is a cell with the name of the skin-markers.
            %   Outputs:
            %     NONE
        
            
            % set variables
            Filename = ResultName;
            Path     = R.ResultsPath;

            q        = R.Subject.q;
            q_t      = R.q_t;
            qdot_t   = R.qdot_t;
            qdot2_t  = R.qdot2_t;


%             % sizes
            [nSamples, nCoords] = size(q_t);
%             
            %check sizes
%             if nCoords ~= size(q,1) || nCoords ~= size(qdot_t,2) || nCoords ~= size(qdot2_t,2)
%                 error('q, q_t, qdot_t and qdot2_t MUST have the same length!!');
%             end
            
            % Calculate extensions
            IndExt = 'ind';  
            PosExt = 'pos';  
            VelExt = 'vel';  
            AccExt = 'acc';  
            SenExt = 'sen';  
            CoMExt = 'com';  
            AcondTrajExt = 'amt';  
            DisExt = 'dis';  
            RawTrajExt = 'rmt';

            % filenames
            FilenamePos = [Filename,'.',PosExt];
            FilenameVel = [Filename,'.',VelExt];
            FilenameAcc = [Filename,'.',AccExt];
            FilenameSen = [Filename,'.',SenExt];
            FilenameCoM = [Filename,'.',CoMExt];
            FilenameAcondTraj = [Filename,'.',AcondTrajExt];
            FilenameErr = [Filename,'.',DisExt];
            FilenameRawTraj = [Filename,'.',RawTrajExt];
            
            % display info
            
            str = sprintf('    Writing output data in COMPAMM format:');
            fprintf(R.ExperLogFileId, '%s\n', str);
            if(R.Settings.Display == 1 || R.Settings.Display == 2)
                disp(str);
            end
            % create file .ind
            %             fid = fopen([Path,Filename,'.',IndExt],'a');
            fid = fopen([Path,Filename,'.',IndExt],'w');
            if R.Settings.Results.Position == 1
                str = sprintf(['      ',FilenamePos]);
                fprintf(R.ExperLogFileId, '%s\n', str);
                if(R.Settings.Display == 1 || R.Settings.Display == 2)
                    disp(str);
                end
                fprintf(fid,[FilenamePos,'\r\n']);
                % Print position file (*.pos)
                R.printCompFile(FilenamePos, q_t, R.Subject.q, 'POSITION');
            end
            if R.Settings.Results.Velocity == 1
                str = sprintf(['      ',FilenameVel]);
                fprintf(R.ExperLogFileId, '%s\n', str);
                if(R.Settings.Display == 1 || R.Settings.Display == 2)
                    disp(str);
                end
                fprintf(fid,[FilenameVel,'\r\n']);
                % Print velocity file (*.vel)
                R.printCompFile(FilenameVel, qdot_t, R.Subject.q, 'VELOCITY');
            end
            if R.Settings.Results.Acceleration == 1
                str = sprintf(['      ',FilenameAcc]);
                fprintf(R.ExperLogFileId, '%s\n', str);
                if(R.Settings.Display == 1 || R.Settings.Display == 2)
                    disp(str);
                end
                fprintf(fid,[FilenameAcc,'\r\n']);
                % Print acceleration file (*.acc)
                R.printCompFile(FilenameAcc, qdot2_t, R.Subject.q, 'ACCELERATION');
            end
            if R.Settings.Results.Sensor == 1 && ~isempty(R.Subject.Sensors)
                str = sprintf(['      ',FilenameSen]);
                fprintf(R.ExperLogFileId, '%s\n', str);
                if(R.Settings.Display == 1 || R.Settings.Display == 2)
                    disp(str);
                end
                fprintf(fid,[FilenameSen,'\r\n']);
                % Print sensor file (*.sen)
                R.printCompFile(FilenameSen,R.SensorValues, R.SensorNames, 'SENSOR')
% %               % Solo para guardar los sensores
%                 NAngles = size(Data,2);
%                 j=1;
%                 for i=1:NAngles
%                     if(abs(Data(11,i))>0.01) 
% %                         AngleName = ['ang',num2str(j)];
%                         AngleName = ['ang',Var{i}];
%                         Sensors.(AngleName)= Data(:,i);
%                         j=j+1;
%                     end
%                 end
%                 save([Path,Filename(1:end-4),'_SensorsNat'],'Sensors');
            end
            if R.Settings.Results.CentreOfMass == 1
                str = sprintf(['      ',FilenameCoM]);
                fprintf(R.ExperLogFileId, '%s\n', str);
                if(R.Settings.Display == 1 || R.Settings.Display == 2)
                    disp(str);
                end
                fprintf(fid,[FilenameCoM,'\r\n']);
                % Print position of CoM (*.com)
                CoM = R.getCoM();
                NSegments = size(R.Subject.Segments,1);
                CoMNames = cell(3*NSegments,1);
                for i=1:NSegments
                    if R.Subject.Segments(i).Fixed ~= 1 % If is not Fixed(Ground)
                        CoMNames{3*(i-1)-2}= [R.Subject.Segments(i).Name,'x'];
                        CoMNames{3*(i-1)-1}= [R.Subject.Segments(i).Name,'y'];
                        CoMNames{3*(i-1)}  = [R.Subject.Segments(i).Name,'z'];
                    end
                end
                R.printCompFile(FilenameCoM, CoM, CoMNames, 'CoM');
            end
            if R.Settings.Results.AcondMarkerTrajectory == 1
                str = sprintf(['      ',FilenameAcondTraj]);
                fprintf(R.ExperLogFileId, '%s\n', str);
                if(R.Settings.Display == 1 || R.Settings.Display == 2)
                    disp(str);
                end
                fprintf(fid,[FilenameAcondTraj,'\r\n']);
                % Print Experimental markers trajectories (*.traj)
                R.printCompFile(FilenameAcondTraj, R.g_t, R.z, 'ACONDITIONED_MARKER_POSITION');
            end
            if R.Settings.Results.RawMarkerTrajectory == 1
                str = sprintf(['      ',FilenameRawTraj]);
                fprintf(R.ExperLogFileId, '%s\n', str);
                if(R.Settings.Display == 1 || R.Settings.Display == 2)
                    disp(str);
                end
                fprintf(fid,[FilenameRawTraj,'\r\n']);
                % Print Experimental markers trajectories (*.traj)
                NRawMarkers = size(R.RawMarkerNames,1);
                for i=1:NRawMarkers
                    RawMarkerCoordNames{3*i-2} = [R.RawMarkerNames{i},'x'];
                    RawMarkerCoordNames{3*i-1} = [R.RawMarkerNames{i},'y'];
                    RawMarkerCoordNames{3*i}   = [R.RawMarkerNames{i},'z'];
                end
                R.printCompFile(FilenameRawTraj, R.RawMarkerCoords, RawMarkerCoordNames, 'RAW_MARKER_POSITION');
            end
            if R.Settings.Results.MarkerError == 1
                str = sprintf(['      ',FilenameErr]);
                fprintf(R.ExperLogFileId, '%s\n', str);
                if(R.Settings.Display == 1 || R.Settings.Display == 2)
                    disp(str);
                end
                fprintf(fid,[FilenameErr,'\r\n']);
                MarkerErrorDist = R.calcMarkerError();
                % Print Error data file. (Type *.dis)
                R.printCompFile(FilenameErr, MarkerErrorDist, R.MarkerNames, 'DIST_ExpMARK_&_SegMARK');
            end
            fclose(fid);
            
        end
        function result_Ramsis(R,ResultName)
            % set variables
            Filename = ResultName;
            Path     = R.ResultsPath;
            Segments  = R.Subject.Segments;
            Sensor    = R.Subject.Sensors;
            q        = R.Subject.q;
            q_t      = R.q_t;
             % constants
            rad2deg = 180/pi;
            
            % sizes
            nSensors = size(Sensor,1);
            [nSamples, nCoords] = size(q_t);
            
            %check sizes
            if nCoords ~= length(q)
                error('q and q_t MUST have the same length!!');
            end
            % Extension
            KfrExt = 'kfr';  
            % filenames
            FilenameExt = [Filename,'.',KfrExt];
            
            
            % display info
            str = sprintf(['    Writing reconstructed motion in RAMSIS format:\n      ',FilenameExt]);
            fprintf(R.ExperLogFileId, '%s\n', str);
            if(R.Settings.Display == 1 || R.Settings.Display == 2)
                disp(str);
            end
            % Check if the Sensors are for the RAMSIS model:
            %   There must be a single root element called GHZ and
            %   this sensor must be the only sensor of type 'TRS'.
            counter = 0;
            IndexTRSinSensor = 0;
            for i=1:nSensors
                SensorType = Sensor{i}.Type;
                if strcmp(SensorType,'TRS')
                    counter = counter + 1;
                    IndexTRSinSensor = i; % This variable is only used if there is only one TRS sensor
                end
            end
            if counter > 1
                str = sprintf(['function Result_rmr: the cell "Sensor" contains more than one translational sensor\n'...
                    '"Sensor" does not complay with kfr format. An output kfr file cannot be generated']);
                warning(str);
                return
            end
            if IndexTRSinSensor ~= 1
                str = sprintf(['function Result_rmr: the translational sensor must be place in the first position of\n'...
                    'the cell "Sensor". It does not complay with kfr format and an kfr file cannot be generated']);
                warning(str);
                return
            end
            
            
            % calculate the sensor data for every frame
%             [SensorVars, SensorData] = fillSensorData(Body, RefBody, Sensor, q, q_t);
  
            
            % open file
            fid = fopen([Path,FilenameExt],'w');
            
            % Print header file
            fprintf(fid,'ANIMATION_V1.4\r\n');
            fprintf(fid,'20582 20601 20561 21283\r\n');
            fprintf(fid,['0 ',num2str(1/R.Deltat),' 1\r\n']);
            fprintf(fid,[num2str(nSamples),'\r\n']);

            
            % Print data for each step
            for i = 1 : nSamples
                % Print position of the root element GHZ in mm
                fprintf(fid,'%8.2f', 1000*R.SensorValues(i, [1:3])); fprintf(fid,'\r\n');
                % Print the number of rotational sensors
                fprintf(fid,'56\r\n');
                for j=1:nSensors
                    if strcmp(Sensor{j}.Type,'SPH')
                        %Print joint name
                        fprintf(fid,'%-4s', Sensor{j}.Name);
                        % Print joint angles in degrees for each sensor (in the same line)
                        RotOrder = Sensor{j}.RotSeq;
                        if strcmp(RotOrder,'321') % change order to xyz
                            %                         fprintf(fid,' %7.2f', R.SensorValues(i,[3*j, 3*j-1, 3*j-2]) );
                            %                         if strcmp(Sensor{j}.Name,'GSPL')
                            %                             a = R.SensorValues(i,3*j);
                            %                             b = R.SensorValues(i,3*j-1);
                            %                             c = R.SensorValues(i,3*j-2);
                            %                             c = c +90;
                            %                             fprintf(fid,' %7.2f', a,b,c );
                            %                         else
                            fprintf(fid,' %7.2f', R.SensorValues(i,[3*j, 3*j-1, 3*j-2]) );
                            %                         end
                        elseif strcmp(RotOrder,'123')
                            fprintf(fid,' %7.2f', R.SensorValues(i,[3*j-2:3*j]) );
                        end
                        fprintf(fid,'\r\n');
                    end
                end
                
                % Print value of the joints
                % Check if the value is used in the reconstruction if not use the default value
                if ~R.checkSensorName('GHZ') fprintf(fid,'GHZ     0.00    0.00  -25.40\r\n'); end
                if ~R.checkSensorName('GHUL')fprintf(fid,'GHUL  -12.40   11.10   82.10\r\n'); end
                if ~R.checkSensorName('GKNL')fprintf(fid,'GKNL    0.00   63.10    0.00\r\n'); end
                if ~R.checkSensorName('GSPL')fprintf(fid,'GSPL    0.00    0.00   72.00\r\n'); end
                if ~R.checkSensorName('GLK') fprintf(fid,'GLK     0.00    0.00  -14.69\r\n'); end
                if ~R.checkSensorName('GLL') fprintf(fid,'GLL     0.00    0.00   16.60\r\n'); end
                if ~R.checkSensorName('GBL') fprintf(fid,'GBL     0.00    0.00    7.90\r\n'); end
                if ~R.checkSensorName('GBB') fprintf(fid,'GBB     0.00    0.00    9.60\r\n'); end
                if ~R.checkSensorName('GHB') fprintf(fid,'GHB     0.00    0.00   12.50\r\n'); end
                if ~R.checkSensorName('GHH') fprintf(fid,'GHH     0.00    0.00   -5.70\r\n'); end
                if ~R.checkSensorName('GKH') fprintf(fid,'GKH     0.00    0.00   -1.40\r\n'); end
                if ~R.checkSensorName('GFBL')fprintf(fid,'GFBL    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GHUR')fprintf(fid,'GHUR   12.40  -11.10   82.10\r\n'); end
                if ~R.checkSensorName('GKNR')fprintf(fid,'GKNR    0.00   63.10    0.00\r\n'); end
                if ~R.checkSensorName('GSPR')fprintf(fid,'GSPR    0.00    0.00   72.00\r\n'); end
                if ~R.checkSensorName('GFBR')fprintf(fid,'GFBR    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GBRK')fprintf(fid,'GBRK    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GSBL')fprintf(fid,'GSBL    0.00   -7.80    5.70\r\n'); end
                if ~R.checkSensorName('GSL') fprintf(fid,'GSL   -66.90  -33.00   75.60\r\n'); end
                if ~R.checkSensorName('GELL')fprintf(fid,'GELL    4.60  -54.00    0.00\r\n'); end
                if ~R.checkSensorName('GHAL')fprintf(fid,'GHAL    0.00    6.60    8.50\r\n'); end
                if ~R.checkSensorName('GSBR')fprintf(fid,'GSBR    0.00    7.80    5.70\r\n'); end
                if ~R.checkSensorName('GSR') fprintf(fid,'GSR    66.90   33.00   75.60\r\n'); end
                if ~R.checkSensorName('GELR')fprintf(fid,'GELR   -4.60  -54.00    0.00\r\n'); end
                if ~R.checkSensorName('GHAR')fprintf(fid,'GHAR    0.00   -6.60    8.50\r\n'); end
                
                % Print fixed joint values not included in the RAMSIS-REALMAN model
                fprintf(fid,'GAUM    0.00    0.00   -5.00\r\n');
                fprintf(fid,'GD1R   30.00   45.00   30.00\r\n');
                fprintf(fid,'GZ1R    0.00   10.00   -5.00\r\n');
                fprintf(fid,'GM1R    0.00   10.00   -5.00\r\n');
                fprintf(fid,'GR1R    0.00   10.00   -5.00\r\n');
                fprintf(fid,'GK1R    0.00   10.00   -5.00\r\n');
                fprintf(fid,'GD2R    0.00   10.00    0.00\r\n');
                fprintf(fid,'GD3R    0.00    0.00    0.00\r\n');
                fprintf(fid,'GZ2R    0.00   10.00    0.00\r\n');
                fprintf(fid,'GZ3R    0.00    5.00    0.00\r\n');
                fprintf(fid,'GM2R    0.00   10.00    0.00\r\n');
                fprintf(fid,'GM3R    0.00    5.00    0.00\r\n');
                fprintf(fid,'GR2R    0.00   10.00    0.00\r\n');
                fprintf(fid,'GR3R    0.00    5.00    0.00\r\n');
                fprintf(fid,'GK2R    0.00   10.00    0.00\r\n');
                fprintf(fid,'GK3R    0.00    5.00    0.00\r\n');
                fprintf(fid,'GD1L  -30.00   45.00  -30.00\r\n');
                fprintf(fid,'GZ1L    0.00  -10.00   -5.00\r\n');
                fprintf(fid,'GM1L    0.00  -10.00   -5.00\r\n');
                fprintf(fid,'GR1L    0.00  -10.00   -5.00\r\n');
                fprintf(fid,'GK1L    0.00  -10.00   -5.00\r\n');
                fprintf(fid,'GD2L    0.00   10.00    0.00\r\n');
                fprintf(fid,'GD3L    0.00    0.00    0.00\r\n');
                fprintf(fid,'GZ2L    0.00  -10.00    0.00\r\n');
                fprintf(fid,'GZ3L    0.00   -5.00    0.00\r\n');
                fprintf(fid,'GM2L    0.00  -10.00    0.00\r\n');
                fprintf(fid,'GM3L    0.00   -5.00    0.00\r\n');
                fprintf(fid,'GR2L    0.00  -10.00    0.00\r\n');
                fprintf(fid,'GR3L    0.00   -5.00    0.00\r\n');
                fprintf(fid,'GK2L    0.00  -10.00    0.00\r\n');
                fprintf(fid,'GK3L    0.00  -5.00    0.00\r\n');
                % Print sight line length
                fprintf(fid,'1000.00\r\n');
                % Print number of intermediate frames
                fprintf(fid,'0\r\n');
                
            end
            
            % close file
            fclose(fid);
        end
        function result_RamsisDyn(R,ResultName)
            % set variables
            Filename = ResultName;
            Path     = R.ResultsPath;
            Segments  = R.Subject.Segments;
            Sensor    = R.Subject.Sensors;
            q        = R.Subject.q;
            q_t      = R.q_t;
             % constants
            rad2deg = 180/pi;
            
            % sizes
            nSensors = size(Sensor,1);
            [nSamples, nCoords] = size(q_t);
            
            %check sizes
            if nCoords ~= length(q)
                error('q and q_t MUST have the same length!!');
            end
            % Extension
            KfrExt = 'dfr';  
            % filenames
            FilenameExt = [Filename,'.',KfrExt];
            
            
            % display info
            str = sprintf(['   Writing reconstructed motion in RAMSIS format:\n     ',FilenameExt]);
            fprintf(R.ExperLogFileId, '%s\n', str);
            if(R.Settings.Display == 1 || R.Settings.Display == 2)
                disp(str);
            end
            % Check if the Sensors are for the RAMSIS model:
            %   There must be a single root element called GHZ and
            %   this sensor must be the only sensor of type 'TRS'.
            counter = 0;
            for i=1:nSensors
                SensorType = Sensor{i}.Type;
                if strcmp(SensorType,'TRS')
                    counter = counter + 1;
                    IndexTRSinSensor = i; % This variable is only used if there is only one TRS sensor
                end
            end
            if counter > 1
                str = sprintf(['function Result_rmr: the cell "Sensor" contains more than one translational sensor\n'...
                    '"Sensor" does not complay with kfr format. An output kfr file cannot be generated']);
                warning(str);
                return
            end
            if IndexTRSinSensor ~= 1
                str = sprintf(['function Result_rmr: the translational sensor must be place in the first position of\n'...
                    'the cell "Sensor". It does not complay with kfr format and an kfr file cannot be generated']);
                warning(str);
                return
            end
            
            
            % calculate the sensor data for every frame
%             [SensorVars, SensorData] = fillSensorData(Body, RefBody, Sensor, q, q_t);
  
            
            % open file
            fid = fopen([Path,FilenameExt],'w');
            
            % Print header file
            fprintf(fid,'ANIMATION_V1.5\r\n');
            fprintf(fid,'20582 20601 20561 21283\r\n');
            fprintf(fid,['0 ',num2str(1/R.Deltat),' 1\r\n']);
            fprintf(fid,[num2str(nSamples),'\r\n']);

            % Print data for each step
            for i = 1 : nSamples
                % Print position of the root element GHZ in mm
                 % Get global forces and moments
                 % En posición GHZ no hay que poner F y M
%                 Seg1 = R.Subject.Sensors{2}.Segment1;
%                 Glob_R_Seg1 = getRd(Seg1,R.q_t(i,:));
%                 JointU = R.Subject.Sensors{2}.Perm1x;
%                 JointV = R.Subject.Sensors{2}.Perm1y;
%                 JointW = cross(JointU,JointV);
%                 Glob_Forces = (Glob_R_Seg1*[JointU,JointV,JointW])*R.Subject.Joints{1}.F(i,(1:3))';
%                 Glob_Moments =(Glob_R_Seg1*[JointU,JointV,JointW])*R.Subject.Joints{1}.M(i,(1:3))';
                GHZPos = 1000*R.SensorValues(i, [1:3]);
                fprintf(fid,'%8.2f', GHZPos); fprintf(fid,'\r\n');
                % Print the number of rotational sensors
                fprintf(fid,'56\r\n');
                for j=2:nSensors
                    %Print joint name
                    fprintf(fid,'%-4s', Sensor{j}.Name);
                    % Get global forces and moments
%                     Seg1 = R.Subject.Sensors{j}.Segment1;
%                     Glob_R_Seg1 = getRd(Seg1,R.q_t(i,:));
%                     JointU = R.Subject.Sensors{j}.Perm1x;
%                     JointV = R.Subject.Sensors{j}.Perm1y;
%                     JointW = cross(JointU,JointV);
%                     Glob_Forces = (Glob_R_Seg1*[JointU,JointV,JointW])*R.Subject.Joints{j-1}.F(i,(1:3))';
%                     Glob_Moments =(Glob_R_Seg1*[JointU,JointV,JointW])*R.Subject.Joints{j-1}.M(i,(1:3))';
                    % Directly the F y M in the Global frame
                    JointIndex = getVecIndex(Sensor{j}.Name,R.Subject.Joints);
                    Glob_Forces  = R.Subject.Joints{JointIndex}.FGlob(i,(1:3));
                    Glob_Moments = R.Subject.Joints{JointIndex}.MGlob(i,(1:3));
                    % Print joint angles in degrees for each sensor (in the same line)
                    RotOrder = Sensor{j}.RotSeq;
                    if strcmp(RotOrder,'321') % change order to xyz
                        %                         fprintf(fid,' %7.2f', R.SensorValues(i,[3*j, 3*j-1, 3*j-2]) );
                        %                         if strcmp(Sensor{j}.Name,'GSPL')
                        %                             a = R.SensorValues(i,3*j);
                        %                             b = R.SensorValues(i,3*j-1);
                        %                             c = R.SensorValues(i,3*j-2);
                        %                             c = c +90;
                        %                             fprintf(fid,' %7.2f', a,b,c );
                        %                         else
                        AngleForces = [R.SensorValues(i,[3*j, 3*j-1, 3*j-2]),Glob_Forces ,Glob_Moments];
                        fprintf(fid,' %7.2f', AngleForces );
                        %                         end
                    elseif strcmp(RotOrder,'123')
                        AngleForces = [R.SensorValues(i,[3*j-2:3*j]),Glob_Forces ,Glob_Moments];
                        fprintf(fid,' %7.2f', AngleForces);
                    end
                    fprintf(fid,'\r\n');
                end
%                 for j=2:nSensors
%                     % Print joint name
%                     fprintf(fid,'%-4s', Sensor{j,1});
%                     % Print joint angles in degrees for each sensor (in the same line)
%                     % The angle order is always X Y Z independently of the joint rotation order
%                     RotOrder = Sensor{j,4};
%                     if strcmp(RotOrder,'321') % change order to xyz
%                         fprintf(fid,' %7.2f', rad2deg*SensorData(i,[3*j, 3*j-1, 3*j-2]) );
%                     elseif strcmp(RotOrder,'123')
%                         fprintf(fid,' %7.2f', rad2deg*SensorData(i,[3*j-2:3*j]) );
%                     end
%                     fprintf(fid,'\r\n');
%                 end
                
%                 fprintf(fid,[nombreVarQueNoEsta{i},'    0.00    0.00    0.00\r\n']);
                
                % Print value of the joints
                % Check if the value is used in the reconstruction if not use the default value
                if ~R.checkSensorName('GHZ') fprintf(fid,'GHZ     0.00    0.00  -25.40    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GHUL')fprintf(fid,'GHUL  -12.40   11.10   82.10    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GKNL')fprintf(fid,'GKNL    0.00   63.10    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GSPL')fprintf(fid,'GSPL    0.00    0.00   72.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GLK') fprintf(fid,'GLK     0.00    0.00  -14.69    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GLL') fprintf(fid,'GLL     0.00    0.00   16.60    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GBL') fprintf(fid,'GBL     0.00    0.00    7.90    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GBB') fprintf(fid,'GBB     0.00    0.00    9.60    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GHB') fprintf(fid,'GHB     0.00    0.00   12.50    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GHH') fprintf(fid,'GHH     0.00    0.00   -5.70    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GKH') fprintf(fid,'GKH     0.00    0.00   -1.40    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GFBL')fprintf(fid,'GFBL    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GHUR')fprintf(fid,'GHUR   12.40  -11.10   82.10    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GKNR')fprintf(fid,'GKNR    0.00   63.10    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GSPR')fprintf(fid,'GSPR    0.00    0.00   72.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GFBR')fprintf(fid,'GFBR    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GBRK')fprintf(fid,'GBRK    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GSBL')fprintf(fid,'GSBL    0.00   -7.80    5.70    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GSL') fprintf(fid,'GSL   -66.90  -33.00   75.60    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GELL')fprintf(fid,'GELL    4.60  -54.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end
                if ~R.checkSensorName('GHAL')fprintf(fid,'GHAL    0.00    6.60    8.50    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end    
                if ~R.checkSensorName('GSBR')fprintf(fid,'GSBR    0.00    7.80    5.70    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end    
                if ~R.checkSensorName('GSR') fprintf(fid,'GSR    66.90   33.00   75.60    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end    
                if ~R.checkSensorName('GELR')fprintf(fid,'GELR   -4.60  -54.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end    
                if ~R.checkSensorName('GHAR')fprintf(fid,'GHAR    0.00   -6.60    8.50    0.00    0.00    0.00    0.00    0.00    0.00\r\n'); end 
                % Print fixed joint values not included in the RAMSIS-REALMAN model
                fprintf(fid,'GAUM    0.00    0.00   -5.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GD1R   30.00   45.00   30.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GZ1R    0.00   10.00   -5.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GM1R    0.00   10.00   -5.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GR1R    0.00   10.00   -5.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GK1R    0.00   10.00   -5.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GD2R    0.00   10.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GD3R    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GZ2R    0.00   10.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GZ3R    0.00    5.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GM2R    0.00   10.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GM3R    0.00    5.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GR2R    0.00   10.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GR3R    0.00    5.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GK2R    0.00   10.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GK3R    0.00    5.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GD1L  -30.00   45.00  -30.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GZ1L    0.00  -10.00   -5.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GM1L    0.00  -10.00   -5.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GR1L    0.00  -10.00   -5.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GK1L    0.00  -10.00   -5.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GD2L    0.00   10.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GD3L    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GZ2L    0.00  -10.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GZ3L    0.00   -5.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GM2L    0.00  -10.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GM3L    0.00   -5.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GR2L    0.00  -10.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GR3L    0.00   -5.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GK2L    0.00  -10.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                fprintf(fid,'GK3L    0.00   -5.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00\r\n');
                % Print sight line length 
                fprintf(fid,'1000.00\r\n');
                % Print number of intermediate frames
                fprintf(fid,'0\r\n');
            end
            
            % close file
            fclose(fid);
        end
        function result_PAM(R,ResultName)
            % display info
            str = sprintf(['    Writing reconstructed motion in PAM format:\n      ',ResultName,'.xml']);
            fprintf(R.ExperLogFileId, '%s\n', str);
            if(R.Settings.Display == 1 || R.Settings.Display == 2)
                disp(str);
            end
            [NSamples, NCoords] = size(R.q_t);
            NSegments = size(R.Subject.Segments);
            % Create a  XML document.
            docNode = com.mathworks.xml.XMLUtils.createDocument('PAMModelPosture');
            docRootNode = docNode.getDocumentElement;
            SystemElement = docNode.createElement('System');
            docRootNode.appendChild(SystemElement);
            NameAtt = docNode.createAttribute('Name');
            TypeAtt = docNode.createAttribute('Type');
            UniAtt = docNode.createAttribute('UnitSystem');
            SystemElement.setAttributeNode(NameAtt);
            SystemElement.setAttributeNode(TypeAtt);
            SystemElement.setAttributeNode(UniAtt);
            NameAtt.setValue('HM50KR posture');
            TypeAtt.setValue('Human');
            UniAtt.setValue('mm Kg ms Deg');
            for i=1:NSamples
                PostureElement = SystemElement.appendChild(docNode.createElement('Posture'));
                TimeAttribute = docNode.createAttribute('Time');
                TimeAttribute.setValue(num2str(1000*R.Deltat*(i-1)));
                PostureElement.setAttributeNode(TimeAttribute);
                for j=1:NSegments
                    if R.Subject.Segments(i).Fixed ~= 1 % If is not Fixed(Ground)
                        
                        SegmentElement = PostureElement.appendChild(docNode.createElement('Segment'));
                        NameAttribute = docNode.createAttribute('Name');
                        SegmentElement.setAttributeNode(NameAttribute);
                        SegName = regexprep(R.Subject.Segments(j).Name,'_',' ');
                        NameAttribute.setValue(SegName)
                        JointElement = SegmentElement.appendChild(docNode.createElement('Joint'));
                        JNameAttr = docNode.createAttribute('Mnemonic');  JointElement.setAttributeNode(JNameAttr);
                        GRefElement = JointElement.appendChild(docNode.createElement('GlobalReferenceCoordSystem'));
                        GMovElement = JointElement.appendChild(docNode.createElement('GlobalMovingCoordSystem'));
                        URAttr = docNode.createAttribute('U');          GRefElement.setAttributeNode(URAttr);
                        VRAttr = docNode.createAttribute('V');          GRefElement.setAttributeNode(VRAttr);
                        WRAttr = docNode.createAttribute('W');          GRefElement.setAttributeNode(WRAttr);
                        OrRAttr = docNode.createAttribute('Origin');    GRefElement.setAttributeNode(OrRAttr);
                        UMAttr = docNode.createAttribute('U');          GMovElement.setAttributeNode(UMAttr);
                        VMAttr = docNode.createAttribute('V');          GMovElement.setAttributeNode(VMAttr);
                        WMAttr = docNode.createAttribute('W');          GMovElement.setAttributeNode(WMAttr);
                        OrMAttr = docNode.createAttribute('Origin');    GMovElement.setAttributeNode(OrMAttr);
                        NSensors = size(R.Subject.Sensors,1);
                        for k=1:NSensors
                            if ~strcmpi(R.Subject.Sensors{k}.Type,'TRS')
                                if strcmpi(R.Subject.Segments(j).Name,R.Subject.Sensors{k}.Segment2.Name)
                                    JointName = R.Subject.Sensors{k}.Name;
                                    Seg1 = R.Subject.Sensors{k}.Segment1;
                                    Glob_R_Seg1 = getRd(Seg1,R.q_t(i,:));
                                    JointU = R.Subject.Sensors{k}.Perm1x;
                                    JointV = R.Subject.Sensors{k}.Perm1y;
                                    JointW = cross(JointU,JointV);
                                    URValue = num2str((Glob_R_Seg1*JointU)','%-10.6f');
                                    VRValue = num2str((Glob_R_Seg1*JointV)','%-10.6f');
                                    WRValue = num2str((Glob_R_Seg1*JointW)','%-10.6f');
                                    Seg2 = R.Subject.Sensors{k}.Segment2;
                                    Glob_R_Seg2 = getRd(Seg2,R.q_t(i,:));
                                    JointU = R.Subject.Sensors{k}.Perm2x;
                                    JointV = R.Subject.Sensors{k}.Perm2y;
                                    JointW = cross(JointU,JointV);
                                    UMValue = num2str((Glob_R_Seg2*JointU)','%-10.6f');
                                    VMValue = num2str((Glob_R_Seg2*JointV)','%-10.6f');
                                    WMValue = num2str((Glob_R_Seg2*JointW)','%-10.6f');
                                end
                            end
                        end
                        
                        PointIndex = getVecIndex(JointName,R.Subject.Points);
                        PointIndexInq = R.Subject.Points(PointIndex).PosInq;
                        OrgValue = 1000*R.q_t(i,(PointIndexInq:PointIndexInq+2));
                        OrgValue1 = num2str(OrgValue(1),'%-6.2f');
                        OrgValue2 = num2str(OrgValue(2),'%-6.2f');
                        OrgValue3 = num2str(OrgValue(3),'%-6.2f');
                        OrgRValue = [OrgValue1,' ',OrgValue2,' ',OrgValue3];
                        OrgMValue = OrgRValue;
                        if strcmpi(R.Subject.Segments(j).Name,'HM50KR_posture')
                            OrgRValue = '0.00 0.00 0.00';
                        end
                        
                        JNameAttr.setValue(JointName);
                        URAttr.setValue(URValue);
                        VRAttr.setValue(VRValue);
                        WRAttr.setValue(WRValue);
                        OrRAttr.setValue(OrgRValue);
                        UMAttr.setValue(UMValue);
                        VMAttr.setValue(VMValue);
                        WMAttr.setValue(WMValue);
                        OrMAttr.setValue(OrgMValue);
                        if ~isempty(R.Subject.Segments(j).F_Ext)
                            if ~isempty(R.Subject.Segments(j).F_Ext.Value)
                                ForceElement = SegmentElement.appendChild(docNode.createElement('ExternForce'));
                                MnemonicAttr = docNode.createAttribute('Mnemonic');  ForceElement.setAttributeNode(MnemonicAttr);
                                MnemonicAttr.setValue('Pedal');
                                GlobPosAttr = docNode.createAttribute('GlobalPosition');      ForceElement.setAttributeNode(GlobPosAttr);
                                GlobVecAttr = docNode.createAttribute('GlobalVector');      ForceElement.setAttributeNode(GlobVecAttr);
                                %                             FyAttr = docNode.createAttribute('Fy');         ForceElement.setAttributeNode(FyAttr);
                                %                             FzAttr = docNode.createAttribute('Fz');         ForceElement.setAttributeNode(FzAttr);
                                %                             PosxAttr = docNode.createAttribute('Posx');     ForceElement.setAttributeNode(PosxAttr);
                                %                             PosyAttr = docNode.createAttribute('Posy');     ForceElement.setAttributeNode(PosyAttr);
                                %                             PoszAttr = docNode.createAttribute('Posz');     ForceElement.setAttributeNode(PoszAttr);
                                Fx = num2str(R.Subject.Segments(j).F_Ext.Value(i,1)*0.001,'%-7.5f');
                                Fy = num2str(R.Subject.Segments(j).F_Ext.Value(i,2)*0.001,'%-7.5f');
                                Fz = num2str(R.Subject.Segments(j).F_Ext.Value(i,3)*0.001,'%-7.5f');
                                Posx = num2str(R.Subject.Segments(j).F_Ext.Pos(i,1)*1000,'%-6.2f');
                                Posy = num2str(R.Subject.Segments(j).F_Ext.Pos(i,2)*1000,'%-6.2f');
                                Posz = num2str(R.Subject.Segments(j).F_Ext.Pos(i,3)*1000,'%-6.2f');
                                FGlob = [Fx,' ',Fy,' ',Fz];
                                PosGlob = [Posx,' ',Posy,' ',Posz];
                                GlobPosAttr.setValue(PosGlob);
                                GlobVecAttr.setValue(FGlob);
                            end
                        end
                    end
                end
            end
            
            % Save the sample XML document.
            xmlFileName = [R.ResultsPath,ResultName,'.xml'];
            xmlwrite(xmlFileName,docNode);
%             edit(xmlFileName);
        end
        function writeFileMatrices(R,ResultName,SensorMarkers)%###
                       
            % set variables
            Filename = ResultName;
            Path     = R.ResultsPath;
            Deltat   = R.Deltat;
            Segments = R.Subject.Segments;
%             RefBody  = varargin{5};
            q        = R.Subject.q;
            q_t      = R.q_t;
            z        = R.z;
            g_t      = R.g_t;
            MarkerNames = R.MarkerNames; %###
            
            % get number of measured skin-markers and their index in z (or g_t)
            [MarkerIndex, nMarkers] = getSmarkIndex(z, MarkerNames);
            
            % sizes
            nSegments = size(Segments,1);
            [nSamples, nCoords] = size(q_t);
            nMatrices = nSegments + 2 * nMarkers;
            nObjects  = nMatrices; % only one graphic object for each body is allowed.
            
            % filename
            FilenameExt = [Filename,'.sim'];
            % display info
            str = ['      ',FilenameExt,];
            fprintf(R.ExperLogFileId, '%s\n', str);
            if(R.Settings.Display == 1 || R.Settings.Display == 2)
                disp(str);
            end
            
            %----------------------------------------------------------------------------------------------
            % Results file with body matrices (.sim)
            %----------------------------------------------------------------------------------------------
            % open file
            fid = fopen([Path,FilenameExt],'w');
            
            %--------------------------------------------------
            % File header
            %--------------------------------------------------
            if isempty(SensorMarkers)
                fprintf(fid,[num2str(nMatrices),'  ',num2str(nObjects),'\r\n']);
            else
                fprintf(fid,[num2str(nMatrices+9),'  ',num2str(nObjects+9),'\r\n']); %###
            end
            for i = 1 : nSegments
                fprintf(fid,[Segments(i).Name,' ',num2str(i-1),'\r\n']);
            end
            for i = 1 : nMarkers %Auxiliar bodies
                %Reorder q Markers in function o z Markers
                j=1;
                zMarker=z{MarkerIndex(i)};
                while j<nCoords
                    qMarker=q{j};
                    if strcmpi(zMarker, qMarker)==1
                        MarkerIndexq (1,i) = j;
                        break
                    else
                        j=j+1;
                    end
                end
                %---------------------
                Pname=q{MarkerIndexq(i)}; Pname=Pname(1:end-1);
                fprintf(fid,['Aux_',Pname,' ',num2str(nSegments + i-1),'\r\n']);
            end
            for i = 1 : nMarkers
                DrivenCoordName = z{MarkerIndex(i)};  DrivenCoordName = DrivenCoordName(1:end-1);
                fprintf(fid,['MeasuredMarker_', DrivenCoordName,' ',num2str(nSegments + nMarkers + i - 1),'\r\n']);
            end
            if ~isempty(SensorMarkers)
                % For the force position ###
                fprintf(fid,['ForceOrigin ',num2str(nSegments + 2*nMarkers),'\r\n']);
                fprintf(fid,['MarkerXP ',num2str(nSegments + 2*nMarkers+1),'\r\n']);
                fprintf(fid,['MarkerYF ',num2str(nSegments + 2*nMarkers+2),'\r\n']);
                fprintf(fid,['MarkerYB ',num2str(nSegments + 2*nMarkers+3),'\r\n']);
                fprintf(fid,['MarkerZP ',num2str(nSegments + 2*nMarkers+4),'\r\n']);
                fprintf(fid,['MarkerFU ',num2str(nSegments + 2*nMarkers+5),'\r\n']);
                fprintf(fid,['MarkerBU ',num2str(nSegments + 2*nMarkers+6),'\r\n']);
                fprintf(fid,['MarkerXM ',num2str(nSegments + 2*nMarkers+7),'\r\n']);
                fprintf(fid,['MarkerYM ',num2str(nSegments + 2*nMarkers+8),'\r\n']);
            end
            %--------------------------------------------------
            % Print matrices for each step
            %--------------------------------------------------
            
            % Initialize variables to save the T matrixes <=== Ignacio
            T_Segments = cell(nSegments,nSamples);
            T_ModelMarkers = cell(nMarkers,nSamples);
            T_MeasuredMarkers = cell(nMarkers,nSamples);
            Time = cell(1,nSamples);
            maxCoord=[0 0;0 0;0 0]; %Variable to save Maximum and Minimum Coordinate to set it as the axis limit <===== Ignacio
            first=1; %Variable to save first d vector in the maxCord vector, then we compare it with the next one <===== Ignacio
            maxd=zeros(nSegments,1); %Variable to save max d vector coordinate from the T_Segments transformation matrix <=== Ignacio
            
            % get data
            for i = 1 : nSamples
                % Print time
                fprintf(fid,'%10.4f /Time\r\n', Deltat*(i-1));
                % Get time of each step <===== Ignacio
                Time{1,i}=(Deltat*(i-1));
                % Print each body transformation matrix
                for j=1:nSegments
                    [R , Or] = Segments(j).getRd(q_t(i,:));
                    T = [R,Or;zeros(1,3),1];
                    fprintf(fid,'%10.6f %10.6f %10.6f %10.6f\r\n',T');
                    % get T for 3D view in Matlab of segments <===== Ignacio
                    tmp=max(abs(Or));
                    if tmp>maxd(j,1)
                        maxd(j,1)=tmp;
                    end
                    T_Segments{j,i} = T;  
                end
              
                % Print each auxiliar body transformation matrix for MODEL-MARKERS
                for j=1:nMarkers
                    Or = g_t(i,[MarkerIndex(j), MarkerIndex(j)+1, MarkerIndex(j)+2])';
                    SegMarker =   q_t(i,MarkerIndexq(j):MarkerIndexq(j)+2);
                    if isnan(Or(1)) || isnan(SegMarker(1))
                        % This means the the marker is missing. Then to have a realistic
                        % representation the marker is located at the origin and scaled with
                        % a very small factor. Then the marker "dissappears" for the appropriated frames
                        Or = [0; 0; 0];
                        R  = 0.0001*eye(3);                        
                        DataDist(i,j) = 0;
                        
                    else
                        R=zeros(3,3);
                        Distance=eye(3,3);
                        ExperMarker = g_t(i,MarkerIndex(j):MarkerIndex(j)+2);
                        SegMarker =   q_t(i,MarkerIndexq(j):MarkerIndexq(j)+2);
                        DataDist(i,j) = norm(ExperMarker-SegMarker);
                        R(:,1) = q_t(i,[MarkerIndexq(j), MarkerIndexq(j)+1, MarkerIndexq(j)+2])'- g_t(i,[MarkerIndex(j), MarkerIndex(j)+1, MarkerIndex(j)+2])'; 
                        Distance(:,1)=R(:,1);
                        if norm(R(:,1)) == 0  % Measured- and model-marker at the same point. Error = 0!!
                            R = zeros(3); % This is redundant. It is already zeros(3)
                        else
                            R(:,1) = R(:,1)/norm(R(:,1));
                            R(:,2) = anyNormal(R(:,1));
                            R(:,3) = cross(R(:,1), R(:,2));
                            Distance(:,1)=R'*Distance(:,1);
                            R=R*Distance;
                        end
                    end
                    T = [R,Or;zeros(1,3),1];
                    fprintf(fid,'%10.6f %10.6f %10.6f %10.6f\r\n',T');
                    % get T for 3D view in Matlab of model markers <===== Ignacio
                    T_ModelMarkers{j,i} = T;  
                end
                
                % Print each MEASURED marker transformation matrix
                for j=1:nMarkers
                    R = eye(3);
                    d = g_t(i,[MarkerIndex(j), MarkerIndex(j)+1, MarkerIndex(j)+2])';
                    if isnan(d(1))
                        % This means the the marker is missing. Then to have a realistic
                        % representation the marker is located at the origin and scaled with
                        % a very small factor. Then the marker "dissappears" for the appropriated frames
                        d = [0; 0; 0];
                        R = 0.0001*R;
                    end
                    T = [R,d;zeros(1,3),1];
                    fprintf(fid,'%10.6f %10.6f %10.6f %10.6f\r\n',T');
                    % get Maximum and Minimum Coordinate to set it as the axis limit <===== Ignacio
                    if first==1 %Only the first time
                        maxCoord(1,1)=d(1,1);
                        maxCoord(1,2)=d(1,1);
                        maxCoord(2,1)=d(2,1);
                        maxCoord(2,2)=d(2,1);
                        maxCoord(3,1)=d(3,1);
                        maxCoord(3,2)=d(3,1);
                        first=0;
                    end
                    if d(1,1)<maxCoord(1,1)
                        maxCoord(1,1)=d(1,1);
                    end
                    if d(1,1)>maxCoord(1,2)
                        maxCoord(1,2)=d(1,1);
                    end
                    if d(2,1)<maxCoord(2,1)
                        maxCoord(2,1)=d(2,1);
                    end
                    if d(2,1)>maxCoord(2,2)
                        maxCoord(2,2)=d(2,1);
                    end
                    if d(3,1)<maxCoord(3,1)
                        maxCoord(3,1)=d(3,1);
                    end
                    if d(3,1)>maxCoord(3,2)
                        maxCoord(3,2)=d(3,1);
                    end
                    % get T for 3D view in Matlab of model markers <===== Ignacio
                    T_MeasuredMarkers{j,i} = T; 
                end
                % Print the origin of the force %###
                if ~isempty(SensorMarkers)
                    R = eye(3);
                    d = Segments(5).F_Ext.Pos(i,[1,2,3])';
                    T = [R,d;zeros(1,3),1];
                    fprintf(fid,'%10.6f %10.6f %10.6f %10.6f\r\n',T');
                    % Print the markers %###
                    R = eye(3);
                    for j=1:size(SensorMarkers,2)/3
                        if isnan(SensorMarkers(i,3*j-2))
                            d = [0; 0; 0];
                            R = 0.0001*R;
                        else
                            d = SensorMarkers(i,[3*j-2,3*j-1,3*j])';
                        end
                        T = [R,d;zeros(1,3),1];
                        fprintf(fid,'%10.6f %10.6f %10.6f %10.6f\r\n',T');
                    end
                end
            end
            % close file
            fclose(fid);
            
            % save data in .MAT 3D view in Matlab <===== Ignacio
            save([Path,Filename],'T_Segments','T_ModelMarkers','T_MeasuredMarkers','Time','maxCoord','maxd','-append');
            
        end
        function writeFilePlayback(R,Filename,GraphicsPath,SensorMarkers)
            
            %inizialize
            PlaybackPath = R.ResultsPath;
%             GraphicsPath = R.GraphicPath;
            % sizes
            NSegments = size(R.Subject.Segments,1);
            % calculate extension
            PbExt  = 'pb';   % default extension if SOLVER and TOL are not defined
            SimExt = 'sim';  % default extension if SOLVER and TOL are not defined
            
            %             if ~isempty(SOLVER) &&  ~isempty(TOL)
            %     SolExt = getExtSol(SOLVER);
            %     TolExt = getExtTol(TOL);
            %     PbExt  = [TolExt,SolExt,'.',PbExt];
            %     SimExt = [TolExt,SolExt,'.',SimExt];
            % % else
            %     error(['wfile_playback:: both parameters, SOLVER and TOL, must be defined'])
            % filenames
            PlaybackFileExt    = [Filename,'.',PbExt];
            MatricesFileExt    = [Filename,'.',SimExt];
            % get relative path for the graphic objects.
            RelPath = getRelPath(PlaybackPath, GraphicsPath);
            RelPath = getPrintPath(RelPath);
            % display info
            str = sprintf(['      ',PlaybackFileExt]);
            fprintf(R.ExperLogFileId, '%s\n', str);
            if(R.Settings.Display == 1 || R.Settings.Display == 2)
                disp(str);
            end
            
            %Save R and GraphicsPath (.stl files path) in .MAT 3D view in Matlab <===== Ignacio
            save([R.ResultsPath,Filename],'R','GraphicsPath','RelPath');            
            
            
            %--------------------------------------
            % Parameters
            %--------------------------------------
            LmarkRadius = 0.007;
            JointRadius = 0.015;
            MarkerRadius = 0.015;
            WireRadius  = 0.008;
            WireColor   = 'gray';
            MuscleRadius = 0.006;
            AuxWireRadius = 0.004;
            MarkerName_ScaleX = 0.15;
            MarkerName_ScaleY = 0.15;
            MarkerName_ScaleZ = 0.15;
            MarkerName_TransX = 0;
            MarkerName_TransY = 0;
            MarkerName_TransZ = MarkerRadius + 0.005;
            %MarkerName_RotX = -90*3.1416/180; 
            %MarkerName_RotY = -90*3.1416/180;
            MarkerName_RotX = 0; 
            MarkerName_RotY = 0;
            MarkerName_RotZ = 0;
            ModelMarkers=cell(3,1); %Variable to save Model Markers information <===== Ignacio
            ModelMarkers{2,1}=MarkerRadius;
            ModelMarkers{3,1}='b';
            MeasuredMarkers=cell(3,1); %Variable to save Measured Markers information <===== Ignacio
            MeasuredMarkers{2,1}=MarkerRadius;
            MeasuredMarkers{3,1}='g';
            Joints=cell(3,1); %Variable to save Joints information <===== Ignacio
            Joints{2,1}=JointRadius;
            Joints{3,1}='r';
            %----------------------------------------------------------------------------------------------
            % Main Playback file
            %----------------------------------------------------------------------------------------------
            % open file
            fid = fopen([PlaybackPath,PlaybackFileExt],'w');
            % Header
            fprintf(fid,'! ----------------------------------------------------\n');
            fprintf(fid,'! Playback File generated automaticaly - CEIT 2003.\n');
            fprintf(fid,'! ----------------------------------------------------\n');
            fprintf(fid,'Simulation Playback\n');
            fprintf(fid,['ObjectMotion=myMotion File=',MatricesFileExt,'\n\n']);
            
            % Dimension of the graphic objects
            fprintf(fid,'! ----------------------------------------------------\n');
            fprintf(fid,'! Dimensions of the graphic objects\n');
            fprintf(fid,'! ----------------------------------------------------\n');
            fprintf(fid,'! Sphere radius of the bony palpable landmarks\n');
            fprintf(fid,['Parameter = LmarkRadius   Value = ',num2str(LmarkRadius),';','\n']);
            fprintf(fid,'! Sphere radius of each joint centre\n');
            fprintf(fid,['Parameter = JointRadius   Value = ',num2str(JointRadius),';','\n']);
            fprintf(fid,'! Sphere radius of the skin-markers\n');
            fprintf(fid,['Parameter = MarkerRadius   Value = ',num2str(MarkerRadius),';','\n']);
%             fprintf(fid,'! Radius of the cylinder that represents each body segment\n');
%             fprintf(fid,['Parameter = WireRadius   Value = ',num2str(WireRadius),';','\n']);
            fprintf(fid,'!Radius of the cylinder that connects each experimental skin-marker and body skin-marker\n');
            fprintf(fid,['Parameter = AuxWireRadius   Value = ',num2str(AuxWireRadius),';','\n']);
            fprintf(fid,'!Radius of the cylinder that represents each muscles of the model\n');
            fprintf(fid,['Parameter = MuscleRadius   Value = ',num2str(MuscleRadius),';','\n\n']);
            
            % graphics file
            fprintf(fid,'! ----------------------------------------------------\n');
            fprintf(fid,['!Graphics files (bodies, joints, landmarks)','\n']);
            fprintf(fid,'! ----------------------------------------------------\n\n');
            fprintf(fid,'! ----------------------------------------------------\n');
            fprintf(fid,'!Graphic file - Body graphics \n');
            fprintf(fid,'!----------------------------------------------------\n\n');

            WireframeLocalCoords=cell(NSegments,1); %Cell for Wireframe local coordinates <===== Ignacio
            Segments=strings(NSegments,1); %Cell for Segments name <===== Ignacio
            for i=1:NSegments
                BodyName  = R.Subject.Segments(i).Name;
                Segments(i,1)=BodyName; %<===== Ignacio
                P         = R.Subject.Segments(i).LocalPoints;
                M         = R.Subject.Segments(i).LocalMarkers;
                if ~isempty(R.Subject.Segments(i).WireFrame)
                    DrawSeq   = R.Subject.Segments(i).WireFrame.DrawSeq;
                    if ~isempty(R.Subject.Segments(i).WireFrame.Color)
                        WireColor = R.Subject.Segments(i).WireFrame.Color;
                    end
                    if ~isempty(R.Subject.Segments(i).WireFrame.Radius)
                        WireRadius = str2num(R.Subject.Segments(i).WireFrame.Radius);
                    end
                else
                    DrawSeq = [];
                end
                GraphFile = R.Subject.Segments(i).Graphic;
                % Wireframe joining the body points defined in DrawSeq
                % If DrawSeq is empty or has only length 1 nothing is ploted
                nPoints = length(DrawSeq);
                WireframeLocalCoords{i,1}=cell(nPoints,1); % <===== Ignacio
                for j = 1 : (nPoints - 1) % if nPoints > 1
                    Pname_1 = DrawSeq{j};   
                    Pname_2 = DrawSeq{j+1};
                    Draw1Index = getLocVecIndex(Pname_1,P);
                    if Draw1Index == 0
                       Draw1Index = getLocVecIndex(Pname_1,M); 
                       locP1 = R.Subject.Segments(i).LocalMarkers(Draw1Index).LocCoord;
                       WireframeLocalCoords{i,1}{j,1}=locP1; % <===== Ignacio
                    else
                        locP1 = R.Subject.Segments(i).LocalPoints(Draw1Index).LocCoord;
                        WireframeLocalCoords{i,1}{j,1}=locP1; % <===== Ignacio
                    end
                    Draw2Index = getLocVecIndex(Pname_2,P);
                    if Draw2Index == 0
                       Draw2Index = getLocVecIndex(Pname_2,M); 
                       locP2 = R.Subject.Segments(i).LocalMarkers(Draw2Index).LocCoord;
                       WireframeLocalCoords{i,1}{j+1,1}=locP2; % <===== Ignacio
                    else
                       locP2 = R.Subject.Segments(i).LocalPoints(Draw2Index).LocCoord;
                       WireframeLocalCoords{i,1}{j+1,1}=locP2; % <===== Ignacio
                    end
                    fprintf(fid,['! Wireframe graphics between ',Pname_1,'-',Pname_2,' in Body ',BodyName,'\n']);
                    fprintf(fid,['Graphic=',BodyName,'_Wire',Pname_1,Pname_2,' CYLINDER x0=',num2str(locP1(1)), ...
                        '  y0=',num2str(locP1(2)), '  z0=',num2str(locP1(3)), '  x1=',num2str(locP2(1)),...
                        '  y1=',num2str(locP2(2)), '  z1=',num2str(locP2(3)), '  Radius=' ,num2str(WireRadius),';  Facets=40  Material=',WireColor,'\n']);
                    fprintf(fid,['Object=',BodyName,'   Entity=',BodyName,'_Wire',Pname_1,Pname_2,'   Motion=myMotion\n\n']);
                end
                save([R.ResultsPath,Filename],'WireframeLocalCoords','Segments','-append'); %Save Wireframe Local Coordinates to .mat file <===== Ignacio
                if iscell(GraphFile)  &&  ~isempty(GraphFile)
                    nGraphFiles = size(GraphFile,1);
                    fprintf(fid,['! ',BodyName,' graphics\n']);
                    for j = 1 : nGraphFiles
                        % get data from cell
                        GraphFilename = GraphFile{j,1};
                        GraphTrans    = GraphFile{j,2};
                        GraphRot      = GraphFile{j,3};
                        GraphScale    = GraphFile{j,4};
                        
                        % Text strings
                        GraphFile_String  = sprintf(['Graphic=',BodyName,'_GraphFile',num2str(j),'  File=',RelPath,GraphFilename]);
                        if ~isempty(GraphTrans)
                            GraphTrans_String = sprintf(['   Tx=',num2str(GraphTrans(1)),'  Ty=',num2str(GraphTrans(2)),'  Tz=',num2str(GraphTrans(3))]);
                            GraphFile_String = [GraphFile_String, GraphTrans_String];
                        end
                        if ~isempty(GraphRot)
                            GraphRot_String   = sprintf(['   Rx=',num2str(GraphRot(1)),'  Ry=',num2str(GraphRot(2)),'  Rz=',num2str(GraphRot(3))]);
                            GraphFile_String = [GraphFile_String, GraphRot_String];
                        end
                        if ~isempty(GraphScale)
                            GraphScale_String = sprintf(['   Sx=',num2str(GraphScale(1)),'  Sy=',num2str(GraphScale(2)),'  Sz=',num2str(GraphScale(3))]);
                            GraphFile_String = [GraphFile_String, GraphScale_String];
                        end
                        GraphFile_String = getPrintPath(GraphFile_String);
                        fprintf(fid,[GraphFile_String,'\n']);
                        fprintf(fid,['Object=',BodyName,'   Entity = ',BodyName,'_GraphFile',num2str(j),'   Motion = myMotion\n\n']);
                    end
                end  
            end
            fprintf(fid,'! ----------------------------------------------------\n');
            fprintf(fid,'!Graphic file - Joints \n');
            fprintf(fid,'!----------------------------------------------------\n\n');
%             fprintf(fid,['include ',LmarkGraphFileExt,'\n']);

            
            NJoints = size(R.Subject.Joints,1);
            Joints{1,1}=cell(NSegments,1);  %Cell for Joints local coordinates <===== Ignacio
            for i=1:NJoints
                SegName = R.Subject.Joints{i}.Seg2;
                if strcmpi(R.Subject.Joints{i}.Type,'Float')
                    SegIndex = getVecIndex(SegName,R.Subject.Segments);
                    PName   = R.Subject.Segments(SegIndex).LocalPoints(1).Name;
                    locP    = R.Subject.Segments(SegIndex).LocalPoints(1).LocCoord;
                else
                    PName   = R.Subject.Joints{i}.Point2.Point.Name;
                    locP    = R.Subject.Joints{i}.Point2.LocCoord;
                end
                for j=1:NSegments %Loop to save each joint linked to the corresponding segment <===== Ignacio
                    if strcmp(SegName,Segments(j,1))==1
                        if isempty(Joints{1,1}{j,1})
                            Joints{1,1}{j,1}=cell(1,1);
                            Joints{1,1}{j,1}{1,1}=locP;
                            break;
                        else
                            tmp=cell(1,1);
                            tmp{1,1}=locP;
                            Joints{1,1}{j,1}=[Joints{1,1}{j,1};tmp];
                            break;
                        end
                    end
                end
                % header for each body
                fprintf(fid,'!---------------------------------------\n');
                fprintf(fid,['! Body: ',SegName,'\n']);
                fprintf(fid,'!---------------------------------------\n\n');
                
                % The point is a body joint pivot point
                fprintf(fid,['! ',PName,' graphics in Body ',SegName,'\n']);
                fprintf(fid,['Graphic=',SegName,'_',PName,' SPHERE x0=',num2str(locP(1)),'  y0=',num2str(locP(2)), ...
                             '  z0=',num2str(locP(3)),'  Radius=JointRadius;','  Facets=40  Material=red\n']);
                fprintf(fid,['Object=',SegName,'   Entity=',SegName,'_',PName,'   Motion=myMotion\n\n']);
            end
            save([R.ResultsPath,Filename],'Joints','-append'); %Save Joints information to .mat file <===== Ignacio
            fprintf(fid,'! ----------------------------------------------------\n');
            fprintf(fid,'! Markers in the Body(Rigidly fixed to Bodies) and Measured \n');
            fprintf(fid,'!----------------------------------------------------\n\n');

             % Header
            fprintf(fid,'!-------------------------------------------------------------\n');
            fprintf(fid,'! Graphics File:  Markers of the Bodies\n');
            fprintf(fid,'!-------------------------------------------------------------\n\n');
            
            %--------------------------------------------------
            % Graphics for the numbers
            %--------------------------------------------------
            fprintf(fid,['! Generic sphere for measured markers\n']);
            fprintf(fid,['Graphic=MeasuredMarker SPHERE x0=0.0  y0=0.0  z0=0.0  Radius=MarkerRadius;','  Facets=40  Material=green\n\n']);
            %--------------------------------------------------
            % Sphere for each marker rigidly attached
            %--------------------------------------------------
            ModelMarkersLocalCoords=cell(NSegments,1); %Cell for Model Markers local coordinates <===== Ignacio
            R.MarkerNames ={};
            for i=1:NSegments
                BodyName = R.Subject.Segments(i).Name;
                NMarkers = size(R.Subject.Segments(i).LocalMarkers,1);
                ModelMarkersLocalCoords{i,1}=cell(NMarkers,2); % <===== Ignacio
                for j=1:NMarkers
                    MarkerName = R.Subject.Segments(i).LocalMarkers(j).Point.Name;
                    R.MarkerNames = [R.MarkerNames,MarkerName];
                    locM = R.Subject.Segments(i).LocalMarkers(j).LocCoord;
                    fprintf(fid,['! ',MarkerName,' graphics in Body ',BodyName,'\n']);
                    fprintf(fid,['Graphic=',BodyName,'_',MarkerName,' SPHERE x0=',num2str(locM(1)),'  y0=',num2str(locM(2)), ...
                        '  z0=',num2str(locM(3)),'  Radius=MarkerRadius;','  Facets=40  Material=blue\n']);
                    fprintf(fid,['Object=',BodyName,'   Entity=',BodyName,'_',MarkerName,'   Motion=myMotion\n\n']);
                    ModelMarkersLocalCoords{i,1}{j,1}=locM; % <===== Ignacio
                    ModelMarkersLocalCoords{i,1}{j,2}=MarkerName; % <===== Ignacio
                end
            end
            ModelMarkers{1,1}=ModelMarkersLocalCoords;
            save([R.ResultsPath,Filename],'ModelMarkers','-append'); %Save ModelMarkers Local Coordinates to .mat file <===== Ignacio
            %------------------------------------------------------
            % Sphere for each measured Marker (motion capture)
            %------------------------------------------------------
            % Header
            fprintf(fid,'\n');
            fprintf(fid,'!-------------------------------------------------------------\n');
            fprintf(fid,'! Graphics File:  Markers measured with VICON\n');
            fprintf(fid,'!-------------------------------------------------------------\n\n');
            
            % 1) get number of inputs in z
            nInputs = size(R.z,1);
            if ~isempty(SensorMarkers)
                % Para markers del sensor ###
                fprintf(fid,['Graphic=ForceMarker SPHERE x0=0.0  y0=0.0  z0=0.0  Radius=MarkerRadius;  Facets=40  Material=gold\n']);
                fprintf(fid,['Graphic=SensorMarker SPHERE x0=0.0  y0=0.0  z0=0.0  Radius=MarkerRadius;  Facets=40  Material=cyan\n']);
                fprintf(fid,['Object=ForceOrigin  Entity=ForceMarker   Motion=myMotion\n']);
                fprintf(fid,['Object=MarkerXM  Entity=SensorMarker   Motion=myMotion\n']);
                fprintf(fid,['Object=MarkerXP  Entity=SensorMarker   Motion=myMotion\n']);
                fprintf(fid,['Object=MarkerYF  Entity=SensorMarker   Motion=myMotion\n']);
                fprintf(fid,['Object=MarkerYB  Entity=SensorMarker   Motion=myMotion\n']);
                fprintf(fid,['Object=MarkerZP  Entity=SensorMarker   Motion=myMotion\n']);
                fprintf(fid,['Object=MarkerYM  Entity=SensorMarker   Motion=myMotion\n']);
                fprintf(fid,['Object=MarkerFU  Entity=SensorMarker   Motion=myMotion\n']);
                fprintf(fid,['Object=MarkerBU  Entity=SensorMarker   Motion=myMotion\n']);
            end
            % 2) write the graphic object (sphere) for each experimental skin-marker
            % initialize vars
            InputIndex = 1;
            % write in file
            fprintf(fid,['! Spheres\n']);
            MeasuredMarkersLocalCoords=cell(nInputs/3,2); %Cell for Measured Markers local coordinates <===== Ignacio
            while InputIndex < nInputs
                InputName = R.z{InputIndex}; InputName = InputName(1:end-1);
                MarkerIndex = getVecIndex(InputName, R.Markers);
                if MarkerIndex ~= 0
                    MarkerName = R.Markers(MarkerIndex).Name;
                    MeasuredMarkersLocalCoords{MarkerIndex,1}=[0;0;0];  % <===== Ignacio
                    MeasuredMarkersLocalCoords{MarkerIndex,2}=MarkerName;    % <===== Ignacio
                    fprintf(fid,['Object=MeasuredMarker_',MarkerName, '  Entity=MeasuredMarker','   Motion=myMotion\n']);
                    % Add labels to marker spheres ONLY if graphic object exists
                
                    FileExist = exist([GraphicsPath,MarkerName,'.x'],'file');
                    if FileExist ~=2 % file does not exists
                        % Add here a error or warning if required
                    else                        
                        fprintf(fid,['Graphic=MeasuredMarker_',MarkerName,'    File=',RelPath,MarkerName,'.x', ...
                            '  Sx=',num2str(MarkerName_ScaleX), ' Sy=',num2str(MarkerName_ScaleY), ' Sz=',num2str(MarkerName_ScaleZ),...
                            ' Tx=',num2str(MarkerName_TransX),  ' Ty=',num2str(MarkerName_TransY),  ' Tz=',num2str(MarkerName_TransZ),...
                            ' Rx=',num2str(MarkerName_RotX),  ' Ry=',num2str(MarkerName_RotY),  ' Rz=',num2str(MarkerName_RotZ),'\n']);
                        fprintf(fid,['Object=MeasuredMarker_',MarkerName,'   Entity=MeasuredMarker_', MarkerName,'   Motion=myMotion\n']);
                    end
                    InputIndex = InputIndex + 2;
                end
                InputIndex = InputIndex + 1;
            end
            MeasuredMarkers{1,1}=MeasuredMarkersLocalCoords;
            save([R.ResultsPath,Filename],'MeasuredMarkers','-append'); %Save MeasuredMarkers Local Coordinates to .mat file <===== Ignacio
            %AUXILIAR BODIES
            %---------------------------
            %Cylinder between each experimental skin marker and each body skin marker
            AuxWire=cell(1,3);  %Cell to save AuxWire data to .mat file <===== Ignacio
            AuxWire{1,1}=AuxWireRadius;
            AuxWire{1,2}=40;
            AuxWire{1,3}='r';
            save([R.ResultsPath,Filename],'AuxWire','-append'); %Save AuxWire data to .mat file <===== Ignacio
            fprintf(fid,['!Cylinder between each experimental skin marker and each body skin marker ','\n']);
            fprintf(fid,['Graphic=Aux_Wireframe CYLINDER x0=0  y0=0  z0=0  x1=1  y1=0  z1=0  Radius=AuxWireRadius;','  Facets=40  Material=red\n']);
            [MarkeIndex, Nmarks] = getSmarkIndex(R.z, R.MarkerNames); %### vemos para mejorar
            for i=1:Nmarks
                Pname= R.z{MarkeIndex(i)}; Pname=Pname(1:end-1);
                fprintf(fid,['Object=Aux_',Pname,'   Entity=Aux_Wireframe   Motion=myMotion\n\n']);
            end
            %     -------------------------------
%             fprintf(fid,['include ',MarkerGraphFileExt,'\n']);
%             if nargin == 10
%                 fprintf(fid,['include ',MuscleGraphFileExt,'\n\n']);
%             else
%                 fprintf(fid,'\n');
%             end
            
            % close file
            fclose(fid);
        end
        function printCompFile(R,Filename, Data, VarNames, Type)
            % PRINTCOMPFILE creates a results file in COMPAMM format. It can be open
            % with COMPGRAPH
            %
            %   printCompFile(Filename, Path, Var, Data, Deltat, Type)
            %   Inputs:
            %     + Filename is the name of the the result file - string
            %     + Data is a double array(nFrames x nVars) with the value of each
            %       variable (generalized coordinate) at each sample time. For the sample
            %       time i, Datat(i,:) containts the value of each variable. The order of
            %       the variables is the same as in Var.
            %     + Type is the type of results file - string. It can be 'POSITION',
            %       'VELOCITY', 'ACCELERATION', etc.
            %   Outputs:
            %     NONE
            
            % Initialization
            Deltat= R.Deltat;
            Path  = R.ResultsPath;
 
            % sizes
            [NSamples, NVars] = size(Data);
%             NCols = size(Var,2);
            
            % open file
            fid = fopen([Path,Filename],'w');
            
            % Print header file
            fprintf(fid,[Type,' ',num2str(NVars+1),'\n']);
            fprintf(fid,'INDEX_BEGIN\n');
            fprintf(fid,'TIME 1\n');
            for i = 1 : NVars
%                 if NCols == 1
                    % Var is q, vector of generalized coordinates
                    fprintf(fid,[char(VarNames(i)),' 1\n']);
%                 else
%                     % VarNames is Muucle, cell with muscle data
%                     fprintf(fid,[VarNames{i,1},' 1\n']);
%                 end
            end
            fprintf(fid,'INDEX_END\n');
            
            % Print data for each step
            for i = 1 : NSamples
                % Print time
                fprintf(fid,'%10.4f\n', Deltat*(i-1));
                % Print each coordinate value
                fprintf(fid,'%10.6f\n',Data(i,:));
            end
            
            % close file
            fclose(fid);
        end
    end
    
end

