classdef EXPERIMENT <handle
    %EXPERIMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Model           % Model information                 Model.Path
                        %                                   Model.File / Name of model.m
        Experiment      % Motion and Subject Data Files     Experiment.Path / Experiment.Name
                        %                                   Experiment.DataFiles {i,1} SubDataFile / {i,2} Motion Files {j} CSV
        SolverType      % Type of solver used to resolve    char
        SubDataPath     % Subject Data path                 char
        MotionPath      % Motion path                       char
        ResultPath      % Result path                       char
        ModelEqPath     % Model Equations path              char
        AddCtrPath      % Additional Constraints  file Path char
        GuidedCoordPath % Guided Coordinates file Path      char
        GraphicPath     % Graphic path                      char
        Settings        % Struct that contains settings for motion reconstruction
                        % Settings.Type. Inverse Kinematics or Inverse Dynamics. Possible options: IK, ID for IK+ID
                        % Settings.Display. Define ammount of feedback to user: 0(little information) 1(complete information)
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
        function E = EXPERIMENT(varargin)
            
            global PathBar

            if nargin == 3 || nargin == 4
                E.Model =      varargin{1};
                E.Experiment = varargin{2};
                E.SolverType = 'weightedOTM';
                E.Settings = varargin{4};
                
                % Model.Guided contains a file name (*.m) Here the extension is removed
                E.Model.Guided = getFilenameAndExt(E.Model.Guided);

                % Experiment does not contain the experiment name. Here the experiment name
                % is taken from the experiment path.                
                BarPos = findstr(E.Experiment.Path,PathBar);
                E.Experiment.Name = [E.Experiment.Path(BarPos(end-1)+1:BarPos(end)-1)];
                
            end
            if  nargin == 4
                E.SolverType = varargin{3};
            end
            if isfield(E.Experiment,'Path')
                ExperimentPath = E.Experiment.Path;
                E.SubDataPath  = [ExperimentPath,'Subjects',PathBar];
                % Model Equations Path
                E.ModelEqPath = [ExperimentPath,'ModelEquations',PathBar];
                % Additional Constraints Path
                E.AddCtrPath = [ExperimentPath,'AdditionalEquations',PathBar];
                % Guided Coordinates Path
                E.GuidedCoordPath = [ExperimentPath,'AdditionalEquations',PathBar];
            else
                if isfield(E.Experiment,'SubjectPath')
                E.SubDataPath  = E.Experiment.SubjectPath;
                else
                    error('The struct "Experiment" does not contain field "SubjectPath"')
                end
                if isfield(E.Experiment,'MotionPath')
                E.MotionPath  = E.Experiment.MotionPath;
                else
                    error('The struct "Experiment" does not contain field "MotionPath"')
                end
                if isfield(E.Experiment,'ResultPath')
                E.ResultPath  = E.Experiment.ResultPath;
                else
                    error('The struct "Experiment" does not contain field "ResultPath"')
                end
                if isfield(E.Experiment,'AdditionalEquationsPath')
                E.AddCtrPath  = E.Experiment.AdditionalEquationsPath;
                else
                    error('The struct "Experiment" does not contain field "AdditionalEquationsPath"')
                end
                if isfield(E.Experiment,'ModelEquationsPath')
                E.ModelEqPath  = E.Experiment.ModelEquationsPath;
                else
                    error('The struct "Experiment" does not contain field "ModelEquationsPath"')
                end 
            end
            E.GraphicPath  = [E.Model.Path,'graphics',PathBar];
        end
        function doMotRec(E)
            
            global PathBar
            
            if ~(isdeployed)
                loadPath(E.ModelEqPath,E.AddCtrPath) %,E.GuidedCoordPath);
            end

            ExperVar.z =[];
            ExperVar.Wm= [];
            ExperVar.Markers = [];
            ExperVar.Ws = [];
            ExperVar.gs_t = [];
            ExperVar.PoszInq = [];
            ExperVar.FilePhiName = [];
            ExperVar.FilePhiqName = [];
            NSubjects = size(E.Experiment.DataFiles,1);
%             try
            if NSubjects == 0
                try
                    
                    % If default IK & ID solvers are used
                    if strcmpi(E.Settings.SolverIK, 'SODERKVIST')
                        % If not default IK solvers is used
                        ResultsPath = [E.Experiment.Path,'Results_Soderkvist',PathBar];
                    else
                        % If default IK & ID solvers (OTM & Recursive) are used
                        ResultsPath = [E.Experiment.Path,'Results',PathBar];
                    end
                    
                    if ~exist(ResultsPath,'dir')
                        mkdir(ResultsPath)
                    end
                    ExperLogFileId = fopen([ResultsPath,E.Experiment.Name,'_',E.Model.File(1:end-4),'.log'],'w');
                    Subject = E.getSubject(ExperLogFileId);
                    ExperVar.Par = [];
                    Motion = [];
                    PressurePath = [];
                    a=tic;
                    Reconstruction = RECONSTRUCTION(Subject,Motion,PressurePath,E.SolverType,E.ModelEqPath,E.AddCtrPath,E.Experiment.Path,...
                        E.Experiment.Name,ExperVar,ExperLogFileId,E.Settings,E.Model);
                    ExperVar = Reconstruction.doIK(E.Model.Guided,ResultsPath,E.GraphicPath);
                    if strcmpi(E.Settings.Type,'ID')
                        Reconstruction.doID(ResultsPath);
                    end
                    [H, MI, S] = second2HMS(toc(a));
                    str =   ['  Total Time: ',num2str(H),' hour(s) ',num2str(MI),' minute(s) and ',num2str(S),' second(s).'];
                    disp(str);
                    fprintf(ExperLogFileId,'%s\n', str);
                    fclose(ExperLogFileId);
                    
                catch ME
                    FirstLine = sprintf(['\n===================================================================================================\n', ...
                        ' ERROR during motion reconstruction:\n']);
                    
                    %ErrorMessage = ['   ',ME.message]; % for Matlab 2010b
                    ErrorMessage = sprintf(['   ',ME.message]); % for Matlab 2012b, 2014a 
                                        
                    NErrors=0;
                    if E.Settings.Display == 2
                        NErrors = length(ME.stack); % for debugging
                    end
                    
                    for k=1:NErrors
                        ErrorMessage = [ErrorMessage, sprintf(['\n      ',ME.stack(k).name,'  line:',num2str(ME.stack(k).line)])];
                    end
                    LastLine = sprintf('\n===================================================================================================\n');
                    FullMessage = [FirstLine,ErrorMessage,LastLine];
                    error(FullMessage); fprintf(ExperLogFileId, '%s\n', FullMessage);
                    
                    % close .log file
                    fclose(ExperLogFileId);                    
                end
            else
                for i=1:NSubjects
                    try
                        SubjectFile = E.Experiment.DataFiles{i,1};
                        if ~isempty(SubjectFile)
                            SubjecPath = [E.SubDataPath,SubjectFile(1:(end-4)),PathBar];
                            Motion.Path = [SubjecPath,'Motions',PathBar];
                            PressurePath = [SubjecPath,'PressureMaps',PathBar];
                            % If default IK & ID solvers are used
                            if strcmpi(E.Settings.SolverIK, 'SODERKVIST')
                                % If not default IK solvers is used
                                ResultsPath = [SubjecPath,'Results_Soderkvist',PathBar];
                            else
                                % If default IK & ID solvers (OTM & Recursive) are used
                                ResultsPath = [SubjecPath,'Results',PathBar];
                            end
                            if ~exist(ResultsPath,'dir')
                                mkdir(ResultsPath)
                            end
                            
                            AceptedExtension{1,1} = 'mat'; AceptedExtension{2,1} = 'exp'; AceptedExtension{3,1} = 'psp';
                            AceptedExtension{4,1} = 'rsp'; AceptedExtension{5,1} = 'par';
                            checkFileAndPath(SubjecPath,SubjectFile,AceptedExtension);
                            Subject = E.getSubject(SubjectFile,SubjecPath);
                            Subjectstr = [', subject: ',  SubjectFile(1:(end-4))];
                        else
                            Motion.Path = [E.Experiment.Path,'Motions',PathBar];
                            % If default IK & ID solvers are used
                            if strcmpi(E.Settings.SolverIK, 'SODERKVIST')
                                % If not default IK solvers is used
                                ResultsPath = [E.Experiment.Path,'Results_Soderkvist',PathBar];
                            else
                                % If default IK & ID solvers (OTM & Recursive) are used
                                ResultsPath = [E.Experiment.Path,'Results',PathBar];
                            end
                            if ~exist(ResultsPath,'dir')
                                mkdir(ResultsPath)
                            end
                            PressurePath = [];
                            Subject = E.getSubject();
                            Subjectstr = [];
                            
                        end
                        if strcmpi(E.Experiment.DataFiles{i,2}{1},'AllMotions')
                            Files = getFilesWithExtension(Motion.Path, 'csv', 0);
                            E.Experiment.DataFiles{i,2}= Files;
                        end
                        NMotions = size(E.Experiment.DataFiles{i,2},1);
                        ExperVar.Par = [];
                        
                    catch ME
                        FirstLine = sprintf(['\n===================================================================================================\n', ...
                            ' ERROR during motion reconstruction:\n']);
                        %ErrorMessage = ['   ',ME.message]; % for Matlab 2010b
                        ErrorMessage = sprintf(['   ',ME.message]); % for Matlab 2012b, 2014a
                        
                        NErrors=0;
                        if E.Settings.Display == 2
                            NErrors = length(ME.stack); % for debugging
                        end
                        
                        for k=1:NErrors
                            ErrorMessage = [ErrorMessage, sprintf(['\n      ',ME.stack(k).name,'  line:',num2str(ME.stack(k).line)])];
                        end
                        LastLine = sprintf('\n===================================================================================================\n');
                        ContLine = sprintf('CONTINUE with NEXT subject reconstruction...\n');
                        FullMessage = [FirstLine,ErrorMessage,LastLine,ContLine];
                        error(FullMessage); % stops all simulations
                        % disp(FullMessage); % continue with the next model simulation
                        % continue % uncomment this with disp
                    end
                    
                    for j=1:NMotions
                        try
                            
                            a=tic;
                            Motion.File = E.Experiment.DataFiles{i,2}{j};
                            ExperLogFileId = fopen([ResultsPath,E.Experiment.Name,'_',E.Model.File(1:end-4),'_',Motion.File(1:(end-4)),'.log'],'w');
                            
                            LineTop    = sprintf('----------------------------------------------------------------------------------------------------\n');
                            MinMessage = [' Motion Reconstruction for model: ',E.Model.File(1:(end-4)),Subjectstr,', motion: ',Motion.File(1:(end-4))];
                            LineBottom = sprintf('\n----------------------------------------------------------------------------------------------------');
                            FullMessage =[LineTop,MinMessage,LineBottom];
                            dispInfo(E.Settings.Display,FullMessage,MinMessage,ExperLogFileId)
                            %                             if(E.Settings.Display == 0)
                            %                                 disp(str1);
                            %                             elseif (E.Settings.Display == 1)
                            %                                 disp([str2,str1,str2]);
                            %                             end
                            %                             fprintf(ExperLogFileId, '%s\n', [str2,str1,str2]);
                            
                            Reconstruction = RECONSTRUCTION(Subject,Motion,PressurePath,E.SolverType,E.ModelEqPath,E.AddCtrPath,E.Experiment.Path,...
                                E.Experiment.Name,ExperVar,ExperLogFileId,E.Settings,E.Model);
                            ExperVar = Reconstruction.doIK(E.Model.Guided,ResultsPath,E.GraphicPath);
                            if strcmpi(E.Settings.Type,'ID')
                                Reconstruction.doID(ResultsPath);
                            end
                            [H, MI, S] = second2HMS(toc(a));
                            str =   ['  Total Time: ',num2str(H),' hour(s) ',num2str(MI),' minute(s) and ',num2str(S),' second(s).'];
                            disp(str);
                            fprintf(ExperLogFileId,'%s\n', str);
                            fclose(ExperLogFileId);
                            ExperLogFileId = [];
                            % disconnect the profiler and save profiling results
                            % % %                     ModelProfData = profile('info'); % save profiling results in structure
                            % % %                     % save the profiler in html format. Is static
                            % % %                     profsave(ModelProfData,[ResultsPath,E.Experiment.Name,E.Model.File,Motion.File,'_nosim']);
                            % % %                     save([ResultsPath,E.Experiment.Name,E.Model.File,Motion.File,'_nosim.pro'],'ModelProfData');
                            % % %                     profile off   % stop the profiler, profiling results ARE NOT cleared
                        catch ME
                            
                            FirstLine = sprintf(['\n===================================================================================================\n', ...
                                ' ERROR during motion reconstruction:\n']);
                            %ErrorMessage = ['   ',ME.message]; % for Matlab 2010b
                            ErrorMessage = sprintf(['   ',ME.message]); % for Matlab 2012b, 2014a
                            
                            NErrors=0;
                            if E.Settings.Display == 2
                                NErrors = length(ME.stack); % for debugging
                            end
                            
                            for k=1:NErrors
                                ErrorMessage = [ErrorMessage, sprintf(['\n      ',ME.stack(k).name,'  line:',num2str(ME.stack(k).line)])];
                            end
                            LastLine = sprintf('\n===================================================================================================\n');
                            ContLine = sprintf('CONTINUE with NEXT motion reconstruction\n');
                            %FullMessage = [FirstLine,ErrorMessage,LastLine,ContLine];
                            % sin linea "CONTINUE..."
                            FullMessage = [FirstLine,ErrorMessage,LastLine];
                            
                            fprintf(ExperLogFileId, '%s\n', FullMessage);
                            error(FullMessage);
                            %disp(FullMessage); fprintf(ExperLogFileId, '%s\n', FullMessage);
                            %continue % continue with the next model simulation
                            
                            % close .log file
                            fclose(ExperLogFileId);
                            ExperLogFileId = [];
                            
                        end
                    end
                end
            end
            
            if ~(isdeployed)
                removePath(E.ModelEqPath,E.AddCtrPath);
            end
        end
        function Subject = getSubject(E,varargin)
            % There is subject file but there is not motion file
            if nargin == (2+1) % we are counting variable E!
                SubjectFile     = varargin{1};
                SubjecPath      = varargin{2};
                % There is motion file subject file but there is not
            elseif nargin == (1+1) % we are counting variable E!
                ExperLogFileId  = varargin{1};
            end
            
            if nargin == (1+1)  % we are counting variable E!
                LineTop    = sprintf('----------------------------------------------------------------------------------------------------\n');
                MinMessage = [' Motion Recontruction for model: ',E.Model.File(1:(end-4))];
                LineBottom = sprintf('\n----------------------------------------------------------------------------------------------------');
                FullMessage =[LineTop,MinMessage,LineBottom];
                dispInfo(E.Settings.Display,FullMessage,MinMessage,ExperLogFileId)
                ModelMat = load([E.Model.Path, E.Model.File]);
                Subject = ModelMat.Human;
                Subject.setGraphicPars();
                %             Subject.calcMassMatrices();
            elseif nargin == (0+1)  % we are counting variable E!
                ModelMat = load([E.Model.Path, E.Model.File]);
                Subject = ModelMat.Human;
                Subject.setGraphicPars();
                %             Subject.calcMassMatrices();
            elseif nargin == (2+1) % we are counting variable E!
                ModelMat = load([E.Model.Path, E.Model.File]);
                Subject = ModelMat.Human;
                Subject.SubjectName = SubjectFile(1:(end-4));
                Subject.parseSubjectPars(SubjecPath,SubjectFile);
                %                 Subject.writeEXP(SubjecPath,SubjectFile); % Para los fichero .mat
                
                %             % ---------------------------------------
                %             % display info
                %             % ---------------------------------------
                %             str = sprintf(['------------------------------------------------------------------------------------------------\n', ...
                %                 ' Importing Anthropometric Parameters from file: ',SubjectFile, '\n', ...
                %                 '------------------------------------------------------------------------------------------------\n']);
                %             fprintf(E.ExperLogFileId, '%s\n', str);
                %             if E.InfoType == 1
                %                 disp(str);
                %             else
                %                 disp(' ');
                %                 disp([' Importing Anthropometric Parameters from file: ',SubjectFile]);
                %             end
                %             Subject.calcMassMatrices();
                %                 Subject.drawLCS();
                Subject.setGraphicPars();
                MarkersPhi = Subject.mkMarkerCtrs();
                Subject.mkPointPhiq(MarkersPhi);
            end
        end
    end
    
end

