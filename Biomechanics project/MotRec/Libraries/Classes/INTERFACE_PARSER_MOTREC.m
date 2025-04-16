classdef INTERFACE_PARSER_MOTREC < handle
    % INTERFACE_MOTREC_PARSE_XML reads input data for motion reconstruction
    % in XML file format
    
    properties
        Settings            % Struct with setting for motion reconstruction
        Model               % Struct with model definition data
        Experiment          % Struct with experiment definition data
        InterfaceDOM        % Document Object Model node representing the parsed XML file
        InputPath           % String with path
        InputFile           % String with file name and extension
        ModelPath           % String with model path
        ModelFile           % String with model file
        GraphicPath         % String with graphic path (contains graphic files for model)
        ExpPath             % String with experiment path
        GuidedVarsFile      % String with guided vars file
        Sep                 % String with separator for error messages. Initialization in Constructor
        DefaultSettings
        DefaultModel
    end
    
    methods
        function IPREC = INTERFACE_PARSER_MOTREC()
            % Constructor of class
            IPREC.Sep = sprintf('\n=========================================================================================\n');
            
            % DEFINE HERE DEFAULT VALUES FOR ALL VARIABLES ========
            
            % Default values for Type and InfoDisplay
            IPREC.DefaultSettings.Type = 'IK'; % Inv. Kinematics/Inv. Dynamics. options: IK, ID (ID includes IK+ID)
            IPREC.DefaultSettings.SolverIK = 'OTM';                                                
            IPREC.DefaultSettings.SolverID = 'RECURSIVE';                                                
            IPREC.DefaultSettings.Display = 0; % Define ammount of feedback to user: 0(minimum) 1(standard) 2(complete information)
            
            % Default values for Results
            IPREC.DefaultSettings.Results.Ramsis = 0; % Turns on/off Ramsis generation results.
            IPREC.DefaultSettings.Results.PAM = 0; % Turns on/off Ramsis generation results.
            IPREC.DefaultSettings.Results.InitPosture = 0; % Posture used for initialization in Compamm format (*_t0.pb & *_t0.sim)
            IPREC.DefaultSettings.Results.CompPlayback = 1; % Reconstructed motion in Compamm format (*.pb & *.sim)
            IPREC.DefaultSettings.Results.Position = 0; % Position of all model elements (markers, vectors, points) in Compamm format (*.pos)
            IPREC.DefaultSettings.Results.Velocity = 0; % Velocity of all model elements (markers, vectors, points) in Compamm format (*.vel)
            IPREC.DefaultSettings.Results.Acceleration = 0; % Acceleration of all model elements (markers, vectors, points) in Compamm format (*.acc)
            IPREC.DefaultSettings.Results.MarkerError = 0; % Distance between experimental-markers and model-markers in Compamm format (*.dis)
            IPREC.DefaultSettings.Results.RawMarkerTrajectory = 0; % Experimental marker trajectories in Compamm format (extracted from oringal data) (*.rmt)
            IPREC.DefaultSettings.Results.AcondMarkerTrajectory = 0; % Aconditioned (smoothed & gaps filled) marker trajectories in Compamm format (*.amt)
            IPREC.DefaultSettings.Results.Sensor = 0; % Variables measured by sensors in Compamm format (*.sen)
            IPREC.DefaultSettings.Results.JointEffort = 0; % Forces & torques at all the joints in the model. Includes reactions & motor efforts
            
            % Default values for Interpolation
            IPREC.DefaultSettings.Interpolation.Method = []; % Method for interpolating gaps in marker trajectories. Options:`'linear'
            IPREC.DefaultSettings.Interpolation.NInterFrames = 0; % Threshold (in frames) for interpolating missing markers. Options: any positive integer
            
            % Default values for Smoothing
            IPREC.DefaultSettings.Smoothing.Method = []; % Smoothing method for marker trajectories. Options: 'butter' (butterworth filter)
            IPREC.DefaultSettings.Smoothing.CutFreq = []; % Cutt-off frequency in Hz for butterworth filter

            % Default values for Model definition            
            IPREC.DefaultModel.Guided = [];
             % Default value for graphics path
             % It is defined in function parseModel(). Not possible to define it here
            
        end
        function checkJustOneElement(IPREC, Element, ElementName)
            if Element.getLength == 0
                error(['File "',IPREC.InputFile,'" -> element "',ElementName,'" not defined']);
            elseif Element.getLength > 1
                error(['File "',IPREC.InputFile,'" -> only ONE element "',ElementName,'" can be defined']);
            end
        end
        
        function nElem = checkMaxOneElement(IPREC, Element, ElementName)
            if Element.getLength == 0
                nElem = 0;
            elseif Element.getLength == 1
                nElem = 1;
            elseif Element.getLength > 1
                nElem = 2;
                error(['File "',IPREC.InputFile,'" -> only ONE element "',ElementName,'" can be defined']);
            end
        end
        
        function fillSingleElementValues(IPREC, Element, MemberVar, ElementValueList, InternalValueList)
            ElementValue = char(Element.item(0).getTextContent);
            ElementName  = char(Element.item(0).getTagName);
            nValues = length(ElementValueList);
            Done = 0; ValuesString = [];
            for i=1:nValues
                if strcmpi(ElementValueList{i}, ElementValue)
                    eval(['IPREC.',MemberVar,' = InternalValueList{i};']);
                    Done = 1;
                end
                ValuesString = [ValuesString, ' ',ElementValueList{i}];
            end
            if ~Done
                error(['File "',IPREC.InputFile,'" -> ',ElementName,' = ',ElementValue,' -> valid values are:',ValuesString]);
            end
        end
        
        function fillAttributeValue(IPREC, ElementName, AttribName, AttribValue, MemberVar, AttribValueList, InternalValueList)
            nValues = length(AttribValueList);
            Done = 0; ValuesString = [];
            for i=1:nValues
                if strcmpi(AttribValueList{i}, AttribValue)
                    eval(['IPREC.',MemberVar,' = InternalValueList{i};']);
                    Done = 1;
                end
                ValuesString = [ValuesString, ' ',AttribValueList{i}];
            end
            if ~Done
                error([ElementName,' -> Attribute "',AttribName,'"=',AttribValue,' -> Wrong value, valid values are:',ValuesString]);
            end
        end
        
        function [IsValid, MessageError] = isValidPath(IPREC, Path)
            global PathBar PathBar2
            PathInit = Path(1:3);
            if ~strcmpi(Path(end),PathBar) % check is Path (full or relative) ends with the right character
                IsValid = 0; % invalid path
                MessageError = ['path must finish with "',PathBar2,'".'];
            elseif strcmp(PathInit(1:2),PathBar2) || strcmp(PathInit(2:3),[':',PathBar]) || strcmp(PathInit(1:2),['.',PathBar])
                % IN WINDOWS:
                %   '\\' -> Full path (network path),
                %   ':\') % Drive path (e.g. C:\myfolder\) Full path
                %   '.\' -> Path relative to xml file
                IsValid = 1; % valid path
                MessageError = [];
            else
                MessageError = 'not valid path format';
                IsValid = 0; % invalid path
            end
        end
        function fillPath(IPREC, PathElement, MemberVars, ParentElementName)
            global PathBar
            Path = char(PathElement.item(0).getTextContent);            
            % Check and fix if path ends with \
            if ~strcmpi(Path(end),PathBar)
                Path = [Path,PathBar];
            end

            % check if path is valid
            [ValidPath, ErrorMessage] = IPREC.isValidPath(Path); % can be absolut or relative path
            if ~ValidPath
                error(['File "',IPREC.InputFile,'" -> ',ParentElementName,' -> Path -> ',ErrorMessage]);
            end
            
            % if path is relative build a full path
            if strcmp(Path(1:2),['.',PathBar])
                Path = [IPREC.InputPath, Path(3:end)];
            end
            
            % Check if path exists
            PathExist = exist(Path,'dir');
            if PathExist~=7
                error(['File "',IPREC.InputFile,'" -> ',ParentElementName,' -> Path -> path does not exist']);
            else
                % this allows for single elements to be string or cell
                if iscell(MemberVars), nMemberVars=length(MemberVars); else nMemberVars = 1; MemberVars = {MemberVars}; end
                %                 if nMemberVars == 1 % this allows for single elements to be string or cell
                %                     if ~iscell(MemberVars), MemberVars = {MemberVars}; end
                %                 end
                for i=1:nMemberVars
                    eval(['IPREC.',MemberVars{i},' = ''',Path,''';']);
                end
            end
        end
        function IsValid = isValidFileExt(IPREC, File, ValidExtensions)
            [FileName , Extension] = getFilenameAndExt(File);
            if sum(strcmpi(Extension, ValidExtensions)) > 0
                IsValid = 1;
            else
                IsValid = 0;
            end
        end
        
        function fillFile(IPREC, FileElement, Path, MemberVars, ParentElementName, ValidExtensions)
            
            File = char(FileElement.item(0).getTextContent);
            ElementName  = char(FileElement.item(0).getTagName);
            
            % is a valid extension?
            ValidFileExt = IPREC.isValidFileExt(File, ValidExtensions);
            if ~ValidFileExt
                nExtensions = length(ValidExtensions);
                ValuesString = [];
                for i=1:nExtensions
                    ValuesString = [ValuesString, ' ',ValidExtensions{i}];
                end
                error(['File "',IPREC.InputFile,'" -> ',ParentElementName,' -> ',...
                    ElementName,' -> file "',File,'"\n', '   valid extensions are:',ValuesString]);
            end
            
            % file exists?
            FileExist = exist([Path,File],'file');
            if FileExist ~=2
                PrintPath = getPrintPath(Path);
                error(['File "',IPREC.InputFile,'" -> ',ParentElementName,' -> ',...
                    ElementName,' -> file: "',File,'"\n', '   can not be found in path:\n    ',PrintPath]);
            else
                nMemberVars = size(MemberVars,1);
                if nMemberVars == 1 % this allows for single elements to be string or cell
                    if ~iscell(MemberVars), MemberVars = {MemberVars}; end
                end
                for i=1:nMemberVars
                    eval(['IPREC.',MemberVars{i},' = ''',File,''';']);
                end
            end
        end
        
        function readxml(IPREC,PathXML,FileXML)
            PathXML = checkAndCorrectFileAndPath(PathXML, FileXML, {'xml'});
            IPREC.InputPath = PathXML;
            IPREC.InputFile = FileXML;
            IPREC.InterfaceDOM = xmlread([PathXML,FileXML]);
        end
        
        function [RootElement] = getRootElement(IPREC)
            RootElement = IPREC.InterfaceDOM.getElementsByTagName('MotionReconstruction');
            if RootElement.getLength == 0
                error(['File "',IPREC.InputFile,'" -> wrong syntax for root element']);
            elseif RootElement.getLength > 1
                error(['File "',IPREC.InputFile,'" -> only one element "MotionReconstruction" is allowed']);
            end
            
            % parse attribute Type
            AttribValue = char(RootElement.item(0).getAttribute('Type'));
            if isempty(AttribValue) 
                error('MotionReconstruction -> Type -> attribute compulsory but not provided');
            else
                IPREC.fillAttributeValue('MotionReconstruction', 'Type', AttribValue, 'Settings.Type', {'IK','ID'}, {'IK','ID'});
            end
            
            % parse attribute SolverIK
            AttribValue = char(RootElement.item(0).getAttribute('SolverIK'));
            if isempty(AttribValue) 
                % if Attribute is optional, fill default value:
                %  IPREC.Settings.SolverIK = IPREC.DefaultSettings.SolverIK;
                % If attribute is compulsory, send error:
                error('MotionReconstruction -> SolverIK -> attribute compulsory but not provided');
            else
                IPREC.fillAttributeValue('MotionReconstruction', 'SolverIK', AttribValue, 'Settings.SolverIK', {'SODERKVIST','OTM'}, {'SODERKVIST','OTM'});
            end
            
            % parse attribute SolverID
            AttribValue = char(RootElement.item(0).getAttribute('SolverID'));
            if isempty(AttribValue) 
                % fill default value
                IPREC.Settings.SolverID = IPREC.DefaultSettings.SolverID;
            else
                IPREC.fillAttributeValue('MotionReconstruction', 'SolverID', AttribValue, 'Settings.SolverID', {'RECURSIVE'}, {'RECURSIVE'});
            end
        end
        
        function parseInfoDisplay(IPREC, RootElement)
            InfoElement = RootElement.item(0).getElementsByTagName('InfoDisplay');
            nElement = IPREC.checkMaxOneElement(InfoElement, 'InfoDisplay');
            if nElement == 0
                % fill default value
                IPREC.Settings.Display = IPREC.DefaultSettings.Display;
            elseif nElement == 1
                IPREC.fillSingleElementValues(InfoElement, 'Settings.Display', {'Minimum','Standard','Complete'}, {0,1,2});
            end
        end
        
        function parseResults_MotRec(IPREC, RootElement)
            % Results_MotRec element. Includes elements: Ramsis, PAM & Compamm
            ResultsMotRecElement = RootElement.item(0).getElementsByTagName('Results_MotRec');
            nElement = IPREC.checkMaxOneElement(ResultsMotRecElement, 'Results_MotRec');
            if nElement == 0
                % Not motion reconstrution results asked for, then if
                if strcmpi(IPREC.Settings.Type,'IK') || strcmpi(IPREC.Settings.Type,'ID')
                    error('No motion reconstruction results asked for. Therefore, reconstruction is not performed');
                end
                
            elseif nElement == 1 % fill settings from motion reconstruction analysis               
                % if analysis is SUBPAR then these results cannot be generated
                if strcmpi(IPREC.Settings.Type,'SUBPAR')
                    % warning results defined in Results_MotRec cannot be generated for analysis SUBPAR
                    error('Motion reconstruction results cannot be generated for analysis of Type SUBPAR');
                end
                
                % Ramsis settings
                RamsisElement = ResultsMotRecElement.item(0).getElementsByTagName('Ramsis');
                nElement = IPREC.checkMaxOneElement(RamsisElement, 'Ramsis');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.Ramsis = IPREC.DefaultSettings.Results.Ramsis;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(RamsisElement, 'Settings.Results.Ramsis', {'No','Yes'}, {0,1});
                end
                
                % PAM settings
                PamElement = ResultsMotRecElement.item(0).getElementsByTagName('PAM');
                nElement = IPREC.checkMaxOneElement(PamElement, 'PAM');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.PAM = IPREC.DefaultSettings.Results.PAM;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(PamElement, 'Settings.Results.PAM', {'No','Yes'}, {0,1});
                end
                
                % COMPAMM element. Includes elements: Playback, Position, Velocity, Sensor, etc.
                CompElement = ResultsMotRecElement.item(0).getElementsByTagName('Compamm');
                % InitializationPosture settings
                InitPostureElement = CompElement.item(0).getElementsByTagName('InitializationPosture');
                nElement = IPREC.checkMaxOneElement(InitPostureElement, 'InitializationPosture');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.InitPosture = IPREC.DefaultSettings.Results.InitPosture;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(InitPostureElement, 'Settings.Results.InitPosture', {'No','Yes'}, {0,1});
                end
                
                % Playback settings
                PlaybackElement = CompElement.item(0).getElementsByTagName('Playback');
                nElement = IPREC.checkMaxOneElement(PlaybackElement, 'Playback');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.CompPlayback = IPREC.DefaultSettings.Results.CompPlayback;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(PlaybackElement, 'Settings.Results.CompPlayback', {'No','Yes'}, {0,1});
                end
                
                % Position settings
                PositionElement = CompElement.item(0).getElementsByTagName('Position');
                nElement = IPREC.checkMaxOneElement(PositionElement, 'Position');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.Position = IPREC.DefaultSettings.Results.Position;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(PositionElement, 'Settings.Results.Position', {'No','Yes'}, {0,1});
                end
                
                % Velocity settings
                VelocityElement = CompElement.item(0).getElementsByTagName('Velocity');
                nElement = IPREC.checkMaxOneElement(VelocityElement, 'Velocity');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.Velocity = IPREC.DefaultSettings.Results.Velocity;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(VelocityElement, 'Settings.Results.Velocity', {'No','Yes'}, {0,1});
                end
                
                % Acceleration settings
                AccelerationElement = CompElement.item(0).getElementsByTagName('Acceleration');
                nElement = IPREC.checkMaxOneElement(AccelerationElement, 'Acceleration');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.Acceleration = IPREC.DefaultSettings.Results.Acceleration;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(AccelerationElement, 'Settings.Results.Acceleration', {'No','Yes'}, {0,1});
                end
                
                % MarkerError settings
                MarkerErrorElement = CompElement.item(0).getElementsByTagName('MarkerError');
                nElement = IPREC.checkMaxOneElement(MarkerErrorElement, 'MarkerError');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.MarkerError = IPREC.DefaultSettings.Results.MarkerError;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(MarkerErrorElement, 'Settings.Results.MarkerError', {'No','Yes'}, {0,1});
                end
                
                % RawMarkerTrajectory settings
                RawMarkerTrajectoryElement = CompElement.item(0).getElementsByTagName('RawMarkerTrajectory');
                nElement = IPREC.checkMaxOneElement(RawMarkerTrajectoryElement, 'RawMarkerTrajectory');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.RawMarkerTrajectory = IPREC.DefaultSettings.Results.RawMarkerTrajectory;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(RawMarkerTrajectoryElement, 'Settings.Results.RawMarkerTrajectory', {'No','Yes'}, {0,1});
                end
                
                % AcondMarkerTrajectory settings
                AcondMarkerTrajectoryElement = CompElement.item(0).getElementsByTagName('AcondMarkerTrajectory');
                nElement = IPREC.checkMaxOneElement(AcondMarkerTrajectoryElement, 'AcondMarkerTrajectory');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.AcondMarkerTrajectory = IPREC.DefaultSettings.Results.AcondMarkerTrajectory;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(AcondMarkerTrajectoryElement, 'Settings.Results.AcondMarkerTrajectory', {'No','Yes'}, {0,1});
                end
                
                % Sensor settings
                SensorElement = CompElement.item(0).getElementsByTagName('Sensor');
                nElement = IPREC.checkMaxOneElement(SensorElement, 'Sensor');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.Sensor = IPREC.DefaultSettings.Results.Sensor;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(SensorElement, 'Settings.Results.Sensor', {'No','Yes'}, {0,1});
                end
                
                % JointEffort settings
                JointEffortElement = CompElement.item(0).getElementsByTagName('JointEffort');
                nElement = IPREC.checkMaxOneElement(JointEffortElement, 'JointEffort');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.JointEffort = IPREC.DefaultSettings.Results.JointEffort;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(JointEffortElement, 'Settings.Results.JointEffort', {'No','Yes'}, {0,1});
                end
                
                % -------------------------------------------------------
                % Check Settings incompatibility
                % -------------------------------------------------------
                % If analysis type = IK, not possible to calculate Joint efforts!!
                if strcmpi(IPREC.Settings.Type,'IK') && (IPREC.Settings.Results.JointEffort == 1)
                    error(['Not possible to generate results for JointEffort.\n   ',...
                                  'Set Type to "ID" or results for JointEffort to "No".']);
                end
                
%                 % If analysis type = ID, and results for Joint efforts not demanded send a warning!!
%                 if strcmpi(IPREC.Settings.Type,'ID') && (IPREC.Settings.Results.JointEffort == 0)
%                     error(['Type set to "ID" but results for JointEffort set to "No".\n   ',...
%                                   'Set Type to "IK" or results for JointEffort to "Yes".']);
%                 end
                
                % If analysis type = SODERKVIST, initialization results have no meaning
                if strcmpi(IPREC.Settings.SolverIK,'SODERKVIST') && (IPREC.Settings.Results.InitPosture == 1)
                    WarnMessage = ['InitializationPosture set to "Yes" and Type set to "SODERKVIST". This analysis \n   ',...
                                   '        can not provide that kind of result. InitializationPosture set to "No" automatically.'];
                    dispWarning(WarnMessage,0);
                    IPREC.Settings.Results.InitPosture = 0;
                end                
                
                % If all results are set to NO, send a warning!!
                SettingsVec = [IPREC.Settings.Results.Ramsis; IPREC.Settings.Results.PAM; IPREC.Settings.Results.InitPosture; ...
                    IPREC.Settings.Results.CompPlayback; IPREC.Settings.Results.Position; IPREC.Settings.Results.Velocity; ...
                    IPREC.Settings.Results.Acceleration; IPREC.Settings.Results.MarkerError; IPREC.Settings.Results.RawMarkerTrajectory; ...
                    IPREC.Settings.Results.AcondMarkerTrajectory; IPREC.Settings.Results.Sensor; IPREC.Settings.Results.JointEffort];
                if sum(SettingsVec) == 0
                    error('No results asked for. Therefore, reconstruction is not performed');
                end
            end
            
        end
        function parseResults_SubjectPar(IPREC, RootElement)
            % Results_SubjectPar. Includes elements: ALM, RSP(=EXP) & PSP
            ResultsSubParElement = RootElement.item(0).getElementsByTagName('Results_SubjectPar');
            nElement = IPREC.checkMaxOneElement(ResultsSubParElement, 'Results_SubjectPar');

            if nElement == 0
                % Not Subject parameter results asked for, then if
                if strcmpi(IPREC.Settings.Type,'SUBPAR')
                    error('No subject parameters results asked for. Therefore, SUBPAR analysis is not performed');
                end
                
            elseif nElement == 1 % fill settings from SUBPAR analysis               
                % if analysis is Motion Reconstruction (IK or ID) then these results cannot be generated
                if strcmpi(IPREC.Settings.Type,'IK') || strcmpi(IPREC.Settings.Type,'ID')
                    % error results defined in Results_SubjectPar cannot be generated for analysis IK or ID
                    error('Subject paramter files cannot be generated for analysis of Type IK or ID');
                end
                
                % ALM settings
                ALM_Element = ResultsSubParElement.item(0).getElementsByTagName('ALM');
                nElement = IPREC.checkMaxOneElement(ALM_Element, 'ALM');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.ALM = IPREC.DefaultSettings.Results.ALM;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(ALM_Element, 'Settings.Results.ALM', {'No','Yes'}, {0,1});
                end
                
                % RSP settings
                RSP_Element = ResultsSubParElement.item(0).getElementsByTagName('RSP');
                nElement = IPREC.checkMaxOneElement(RSP_Element, 'RSP');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.RSP = IPREC.DefaultSettings.Results.RSP;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(RSP_Element, 'Settings.Results.RSP', {'No','Yes'}, {0,1});
                end
                
                % PSP settings
                PSP_Element = ResultsSubParElement.item(0).getElementsByTagName('PSP');
                nElement = IPREC.checkMaxOneElement(PSP_Element, 'PSP');
                if nElement == 0
                    % fill default value
                    IPREC.Settings.Results.PSP = IPREC.DefaultSettings.Results.PSP;
                elseif nElement == 1
                    IPREC.fillSingleElementValues(PSP_Element, 'Settings.Results.PSP', {'No','Yes'}, {0,1});
                end
            end
            
        end
        function parseInterpolation(IPREC, RootElement)
            
            % Interpolation element. Includes elements: Method, GapInterpolated
            InterpoElement = RootElement.item(0).getElementsByTagName('Interpolation');
            nElement = IPREC.checkMaxOneElement(InterpoElement, 'Interpolation');
            if nElement == 0
                % warning interpolation not applied
                WarningMessage = 'Interpolation method not defined. Trajectory gaps will not be interpolated';
                dispWarning(WarningMessage, 0)
                
                % fill default value
                IPREC.Settings.Interpolation.Method = IPREC.DefaultSettings.Interpolation.Method;
                IPREC.Settings.Interpolation.NInterFrames = IPREC.DefaultSettings.Interpolation.NInterFrames;
            elseif nElement == 1
                % Method element
                MethodElement = InterpoElement.item(0).getElementsByTagName('Method');
                IPREC.checkJustOneElement(MethodElement, 'Method');
                IPREC.fillSingleElementValues(MethodElement, 'Settings.Interpolation.Method', {'linear'}, {'linear'});
                
                % GapInterpolated element
                GapInterElement = InterpoElement.item(0).getElementsByTagName('GapInterpolated');
                IPREC.checkJustOneElement(GapInterElement, 'GapInterpolated');
                GapInterValue = str2double(GapInterElement.item(0).getTextContent);
                if ~isnan(GapInterValue) && GapInterValue >= 0
                    IPREC.Settings.Interpolation.NInterFrames = round(GapInterValue);
                    % Notice that GapInterpolated is named NInterFrames within the library!!!
                else
                    error(['File "',IPREC.InputFile,'" -> Interpolation -> GapInterpolated -> ',...
                        'valid values are: any positive integer']);
                end                
            end
        end
        function parseSmoothing(IPREC, RootElement)
            
            % Smoothing element. Includes elements: Method, CutFeq
            SmoothElement = RootElement.item(0).getElementsByTagName('Smoothing');
            nElement = IPREC.checkMaxOneElement(SmoothElement, 'Smoothing');
            if nElement == 0
                
                % warning smoothing not applied
                WarningMessage = 'Smoothing method not defined. Data will not be filtered';
                dispWarning(WarningMessage, 0)
                
                % fill default values
                IPREC.Settings.Smoothing.Method  = IPREC.DefaultSettings.Smoothing.Method;
                IPREC.Settings.Smoothing.CutFreq = IPREC.DefaultSettings.Smoothing.CutFreq;
                
            elseif nElement == 1
                
                % Method element
                MethodElement = SmoothElement.item(0).getElementsByTagName('Method');
                IPREC.checkJustOneElement(MethodElement, 'Method');
                IPREC.fillSingleElementValues(MethodElement, 'Settings.Smoothing.Method', {'butter'}, {'butter'});
                
                % CutFreq element (filter cutt-off frequency)
                CutFreqElement = SmoothElement.item(0).getElementsByTagName('CutFreq');
                IPREC.checkJustOneElement(CutFreqElement, 'CutFreq');
                CutFreqValue = str2double(CutFreqElement.item(0).getTextContent);
                if ~isnan(CutFreqValue) && CutFreqValue >= 0
                    IPREC.Settings.Smoothing.CutFreq = CutFreqValue;
                else
                    error(['File "',IPREC.InputFile,'" -> Smoothing -> CutFreq -> valid values are: ',...
                        'any positive number']);
                end
            end
        end
        function parseModel(IPREC, RootElement)
            global PathBar
            % Model element. Includes elements: Path, File
            ModelElement = RootElement.item(0).getElementsByTagName('Model');
            IPREC.checkJustOneElement(ModelElement, 'Model');
            
            % ModelPath element
            ModelPathElement = ModelElement.item(0).getElementsByTagName('ModelPath');
            IPREC.checkJustOneElement(ModelPathElement, 'ModelPath');
            IPREC.fillPath(ModelPathElement, {'Model.Path';'ModelPath'}, 'Model');
            
            % ModelFile element
            ModelFileElement = ModelElement.item(0).getElementsByTagName('ModelFile');
            IPREC.checkJustOneElement(ModelFileElement, 'ModelFile');
            IPREC.fillFile(ModelFileElement, IPREC.ModelPath, {'Model.File';'ModelFile'}, 'Model', {'xml','mat'});
            
            % GraphicPath element (optional)
            GraphicPathElement = ModelElement.item(0).getElementsByTagName('GraphicPath');
            nElement = IPREC.checkMaxOneElement(GraphicPathElement, 'GraphicPath');
            if nElement == 0
                % fill default value
                IPREC.Model.GraphicPath = [IPREC.Model.Path,'graphics',PathBar];
            elseif nElement == 1
                IPREC.fillPath(GraphicPathElement, {'Model.GraphicPath';'GraphicPath'}, 'Model');
            end
        end
        function parseExperiment(IPREC, RootElement)
            global PathBar
            % Experiment element. Includes elements: Path(only 1), GuidedVarsFile(only 1), Subject(several possible)
            ExpElement = RootElement.item(0).getElementsByTagName('Experiment');
            IPREC.checkJustOneElement(ExpElement, 'Experiment');
            
            % ExpPath element
            PathElement = ExpElement.item(0).getElementsByTagName('ExpPath');
            IPREC.checkJustOneElement(PathElement, 'ExpPath');
            IPREC.fillPath(PathElement, {'Experiment.Path';'ExpPath'}, 'Experiment');
            
            % GuidedVarsFile element
            GuidedVarsElement = ExpElement.item(0).getElementsByTagName('GuidedVarsFile');
            nElement = IPREC.checkMaxOneElement(GuidedVarsElement, 'GuidedVarsFile');
            if nElement == 0
                % fill default value
                IPREC.Model.Guided = IPREC.DefaultModel.Guided;               
            elseif nElement == 1
                IPREC.fillFile(GuidedVarsElement, [IPREC.ExpPath,'AdditionalEquations',PathBar], {'Model.Guided';'GuidedVarsFile'}, 'Experiment', {'m'});
            end
            
            % Subject element (there can many any number)
            AllSubjectElements = ExpElement.item(0).getElementsByTagName('Subject');
            nSubjects = AllSubjectElements.getLength;
            
            if (nSubjects == 0)  % only 1 subject whose parameters are defined in the model def. file (.XML)
                
                % check that no SubPar elements exist. They are not allowed for this option
                ParamElement = ExpElement.item(0).getElementsByTagName('SubPar');
                if ParamElement.getLength >=1
                    error(['File "',IPREC.InputFile,'" -> Experiment -> \n',...
                        '   element "SubPar" can ONLY be defined inside a "Subject" element.']);
                end
                
                % read motion files
                MotionElement = ExpElement.item(0).getElementsByTagName('Motion');
                nMotions = MotionElement.getLength;
                if nMotions >= 1
                    SubjectMotions = cell(nMotions,1);
                    for j=0:nMotions-1
                        
                        MotionFile_j = char(MotionElement.item(j).getTextContent);
                        
                        % Does it have the right extension?
                        if ~IPREC.isValidFileExt(MotionFile_j, {'csv';'c3d';'mat'})
                            error(['File "',IPREC.InputFile,'" -> Experiment ->  Motion: "', ...
                                MotionFile_j,'"\n','   wrong extension, valid extensions are: C3D CSV MAT']);
                        end
                        
                        % Does file exist in the right path?
                        MotionPath = [IPREC.ExpPath,'Motions',PathBar];
                        FileExist = exist([MotionPath,MotionFile_j],'file');
                        if FileExist ~=2
                            error(['File "',IPREC.InputFile,'" -> Experiment ->  Motion "',...
                                MotionFile_j,'"\n','   file can not be found in path:\n    ',getPrintPath(MotionPath)]);
                        end
                        % Stote file in motion files cell
                        SubjectMotions{j+1} = MotionFile_j;
                    end
                    IPREC.Experiment.DataFiles{1,1} = []; % SubjectPar file
                    IPREC.Experiment.DataFiles{1,2} = SubjectMotions; % All motion in one cell
                else
                    IPREC.Experiment.DataFiles = [];
                    % warning smoothing not applied
                    WarnMessage = ['Experiment -> no Motion files provided. Motion must be defined in GuidedVarsFile.'];
                    dispWarning(WarnMessage, 0);
                    % if results for RawMarkerTrajectory as set to Yes send warning
                    if IPREC.Settings.Results.RawMarkerTrajectory == 1
                        WarnMessage = ['Results_MotRec -> RawMarkerTrajectory set to YES but Motion file not provided.\n',...
                        '           It is not possible to generate these results. RawMarkerTrajectory set to NO automatically'];
                        dispWarning(WarnMessage, 0);
                        IPREC.Settings.Results.RawMarkerTrajectory = 0;
                    end                    
                end                
            else % if (nSubjects >= 1)
                
                for i=0:nSubjects-1
                    
                    % ================================================================
                    % A === Subject parameter definition. Options:
                    % ================================================================
                    % 1) Define RSP/EXP/PSP file (processed parameter file) e.g. <SubPar>01_SD.rsp</SubPar>
                    % 2) Define ALM file (palpation file) e.g. <SubPar>01_SD.alm</SubPar>
                    % 3) Define subject code (if ALM or RSP/EXP/PSP not available <SubPar>01_SD</SubPar>
                    % All subject files (ALM/EXP/RSP/PSP) must be located in subject's folder e.g. ./Subjects/01_SD/ -->
                    SubParElement = AllSubjectElements.item(i).getElementsByTagName('SubPar');
                    nSubjectParFiles = SubParElement.getLength;
                    
                    if nSubjectParFiles == 0 % if element Subject is defined, element SubPar is compulsory
                        error(['File "',IPREC.InputFile,'" -> Experiment -> Subject ',...
                            num2str(i+1),' -> no Subject data defined in "SubPar"']);
                        
                    elseif nSubjectParFiles == 1
                        SubDataFile = char(SubParElement.item(0).getTextContent);
                        % check if SubDataFile has valid extension
                        % Valid extension are:
                        %  EXP (Ramsis parameter file created by Ramsis software)
                        %  RSP (Ramsis parameter file created by CEIT software)
                        %  PSP (PAM parameter file)
                        %  XML (CEIT parameter file. Marker, Point or Vector local coord. are given in CEIT model definition format)
                        %  ALM (Palpation file). Palpation toolbox has to be called to generate EXP/RSP & PSP
                        %  --- (No extension means a SUBJECT CODE is given.
                        %       Then the Palpation toolbox has to be called to generate EXP/RSP & PSP
                        ValidFileExt = IPREC.isValidFileExt(SubDataFile, {'exp';'rsp';'psp';'alm'});
                        [SubjectCode, FileExt] = getFilenameAndExt(SubDataFile);
                        
                        if ~ValidFileExt
                            % check if SubjectCode is provided
                            [FileName , Ext] = getFilenameAndExt(SubDataFile);
                            
                            if isempty(SubjectCode) && ~isempty(FileExt) % subject code is provided
                                SubjectCode = FileExt;
                                disp(sprintf(['  Subject ',SubDataFile,': subject parmeter file not provided, ONLY subject code.\n'...
                                    '    Subject parameter file will be generated using the palpation data if available']));
                            else
                                % neither SubjectCode nor Subject par file with valid extension provided
                                error(['File ',IPREC.InputFile,' -> Experiment -> Subject ',...
                                    num2str(i+1),' ',SubDataFile,' -> valid file ext. (exp rsp psp alm)']);
                            end
                        else
                            % check if files exist in subject folder
                            SubjectPath = [IPREC.ExpPath,'Subjects\',SubjectCode,PathBar];
                            SubjectFileExist = exist([SubjectPath,SubDataFile],'file');
                            if SubjectFileExist ~=2
                                error(['File "',IPREC.InputFile,'" -> Experiment -> Subject ',...
                                    SubjectCode,' -> \n','     file ',SubDataFile,' can not be found in path: ',...
                                    getPrintPath(SubjectPath)]);
                            end
                        end
                        % SubDataFile is stored at the end.
                        
                    elseif nSubjectParFiles > 1
                        error(['File "',IPREC.InputFile,'" -> Experiment -> Subject ',num2str(i+1),' -> ',...
                            'only 1 element "SubPar" per subject']);
                    else
                        error(['File "',IPREC.InputFile,'" -> Experiment -> Subject ',num2str(i+1),' -> ',...
                            'unexpected error for element "SubPar"']);
                    end % end of Subject parameters
                    
                    % ================================================================
                    % B === Motion files for subject (CSV)
                    % ================================================================
                    SubjectMotionElement = AllSubjectElements.item(i).getElementsByTagName('Motion');
                    nMotionsForSubject = SubjectMotionElement.getLength;
                    SubjectMotions = cell(nMotionsForSubject,1);
                    % first motion file
                    MotionFile = char(SubjectMotionElement.item(0).getTextContent);
                    
                    if  nMotionsForSubject == 0
                        error(['File "',IPREC.InputFile,'" -> Experiment -> Subject ',...
                            SubjectCode,' -> no motions defined']);
                        
                    elseif nMotionsForSubject == 1 && (strcmp(MotionFile,'ALL_IN_MOTIONS_FOLDER') || ...
                            strcmp(MotionFile,'ALL_C3D_IN_MOTIONS_FOLDER') || strcmp(MotionFile,'ALL_CSV_IN_MOTIONS_FOLDER'))
                        
                        MotionPath = [IPREC.ExpPath,'Subjects',PathBar,SubjectCode,PathBar,'Motions',PathBar];
                        
                        if strcmp(MotionFile,'ALL_C3D_IN_MOTIONS_FOLDER')
                            FilesC3D = getFilesWithExtension(MotionPath, 'c3d', 0);
                            SubjectMotions = FilesC3D;
                            
                        elseif strcmp(MotionFile,'ALL_CSV_IN_MOTIONS_FOLDER')
                            FilesCSV = getFilesWithExtension(MotionPath, 'csv', 0);
                            SubjectMotions = FilesCSV;
                            
                        elseif strcmp(MotionFile,'ALL_IN_MOTIONS_FOLDER')
                            % All motions in folder have been selected
                            % Option 1) Fill SubjectMotions with 'AllMotions' and
                            % pass the job to another function. This option is already implemented
                            %SubjectMotions = {'AllMotions'};
                            % Option 2) Search for files right here in this function and
                            % fill SubjectMotions with the motions. This is compatible with existing functions
                            FilesCSV = getFilesWithExtension(MotionPath, 'csv', 0);
                            FilesC3D = getFilesWithExtension(MotionPath, 'c3d', 0);
                            SubjectMotions = [FilesCSV; FilesC3D];
                        end
                        if isempty(SubjectMotions)
                            error(['File "',IPREC.InputFile,'" -> Experiment -> Subject ',...
                                SubjectCode,' -> Motion "',MotionFile,'"\n   -> ',...
                                'no valid motion files found']);
                        end
                        
                    else % if nMotionsForSubject > 1 OR
                        %    nMotionsForSubject == 1 && "Motion" does not contains keyword *_IN_MOTIONS_FOLDER
                        for j=0:nMotionsForSubject-1
                            
                            MotionFile_j = char(SubjectMotionElement.item(j).getTextContent);
                            MotionPath_j = [IPREC.ExpPath,'Subjects',PathBar,SubjectCode,PathBar,'Motions',PathBar];
                            
                            if strcmp(MotionFile_j,'ALL_IN_MOTIONS_FOLDER') || ...
                                    strcmp(MotionFile_j,'ALL_C3D_IN_MOTIONS_FOLDER') || ...
                                    strcmp(MotionFile_j,'ALL_CSV_IN_MOTIONS_FOLDER')
                                error(['File "',IPREC.InputFile,'" -> Experiment -> Subject ',...
                                    SubjectCode,' -> Motion "',MotionFile_j,'"\n   -> ',...
                                    'When option ',MotionFile_j,' is used there can only be 1 "Motion" element']);
                            end
                            
                            % Does it have the right extension?
                            if ~IPREC.isValidFileExt(MotionFile_j, {'csv';'c3d'})
                                error(['File "',IPREC.InputFile,'" -> Experiment -> Subject ',...
                                    SubjectCode,'\n   -> Motion "',MotionFile_j,'" -> wrong extension, valid ext. are: C3D CSV']);
                            end
                            
                            % Does file exist in the right path?
                            FileExist = exist([MotionPath_j,MotionFile_j],'file');
                            if FileExist ~=2
                                error(['File "',IPREC.InputFile,'" -> Experiment -> Subject ',...
                                    SubjectCode,'\n   -> Motion file "',MotionFile_j,'" -> file can not be found in path \n    ',...
                                    getPrintPath(MotionPath_j)]);
                            end
                            % Stote file in motion files cell
                            SubjectMotions{j+1} = MotionFile_j;
                        end
                    end
                    IPREC.Experiment.DataFiles{i+1,1} = SubDataFile; % SubjectPar file
                    IPREC.Experiment.DataFiles{i+1,2} = SubjectMotions; % All motion in one cell
                    
                    % Options IMPLEMENTED.
                    %  1) Several subjects (par file or subjet code) with several motions(csv,c3d) per subject
                    %  2) No subject par file and model data defined in XML model file
                    %     Then, No Subject element and Motions given at Experiment element level
                    %  3) No subject par file & no motions file. The motion must be defined in file *_guidedVars
                    %  4) Consider options: ALL_IN_MOTIONS_FOLDER, ALL_C3D_IN_MOTIONS_FOLDER,  ALL_CSV_IN_MOTIONS_FOLDER
                    
                end % for i nSubject
            end % if
            
        end % function parseExperiment()
    end % method
    
end % class

