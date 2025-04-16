classdef INTERFACE_PARSER_MOTSIM < handle
    % INTERFACE_PARSER_MOTSIM reads input data for motion simulation
    % in XML file format (SID-simulation interface definition file)
    
    properties
        Settings            % Struct with setting for motion reconstruction
        Model               % Struct with model definition data
        Experiment          % Struct with experiment definition data
        InterfaceDOM        % Document Object Model node representing the parsed XML file
        InputPath           % String with path
        InputFile           % String with file name and extension
        ModelPath           % String with model path
        ModelFile           % String with model file
        DefaultSettings
    end
    
    methods
        function IPSIM = INTERFACE_PARSER_MOTSIM()
            
            % DEFINE HERE DEFAULT VALUES FOR ALL VARIABLES ========
            
            % Default values for Type and InfoDisplay
            IPSIM.DefaultSettings.Type = 'IK'; % Inverse Kinematics/Inverse Dynamics/Subject parameter estimation.
                                                % Possible options: IK, ID, SUBPAR, STATIC. ID includes IK+ID
            IPSIM.DefaultSettings.Display = 0; % Define ammount of feedback to user: 0(minimum) 1(standard) 2(complete information)
            
            % Default values for Results
            IPSIM.DefaultSettings.Results.CompPlayback = 1; % Reconstructed motion in Compamm format (*.pb & *.sim)
            IPSIM.DefaultSettings.Results.Position = 0; % Position of all model elements (markers, vectors, points) in Compamm format (*.pos)
            IPSIM.DefaultSettings.Results.Velocity = 0; % Velocity of all model elements (markers, vectors, points) in Compamm format (*.vel)
            IPSIM.DefaultSettings.Results.Acceleration = 0; % Acceleration of all model elements (markers, vectors, points) in Compamm format (*.acc)
            IPSIM.DefaultSettings.Results.EndEffTraj = 0; % Experimental marker trajectories in Compamm format (extracted from oringal data) (*.traj)
            IPSIM.DefaultSettings.Results.JointEffort = 0; % Forces & torques at all the joints in the model. Includes reactions & motor efforts
            
            % Default values for Tolerance
            IPSIM.DefaultSettings.Tolerances.TolFun = [];
            IPSIM.DefaultSettings.Tolerances.TolX = [];
            
            % Default values for Weights
            IPSIM.DefaultSettings.Weights.Tau = 1e-06;
            IPSIM.DefaultSettings.Weights.q = 1;
            IPSIM.DefaultSettings.Weights.qTrans = 0;
            IPSIM.DefaultSettings.Weights.EE = 1000;
            IPSIM.DefaultSettings.Weights.FRoot = 0;
            IPSIM.DefaultSettings.Weights.MRoot = 0;
            IPSIM.DefaultSettings.Weights.Power = 0;
            IPSIM.DefaultSettings.Weights.Energy = 0;
        end
        function nElem = checkMaxOneElement(IPSIM, Element)            
            ElementName = char(Element.item(0).getTagName);
            if Element.getLength == 0
                nElem = 0;
            elseif Element.getLength == 1
                nElem = 1;
            elseif Element.getLength > 1
                nElem = 2;
                error(['Only ONE element "',ElementName,'" can be defined']);
            end
        end
        function checkJustOneElement(IPSIM, Element)
            ElementName = char(Element.item(0).getTagName);
            if Element.getLength == 0
                error(['Element "',ElementName,'" not defined']);
            elseif Element.getLength > 1
                error(['Only ONE element "',ElementName,'" can be defined']);
            end
        end
        function [IsValid, MessageError] = isValidPath(IPSIM, Path)            
            global PathBar PathBar2
            
            PathInit = Path(1:3);
            if ~strcmpi(Path(end),PathBar) % check is Path (full or relative) ends with the right character
                IsValid = 0; % invalid path
                MessageError = ['path must finish with ',PathBar2];
            elseif strcmp(PathInit(1:2),PathBar2) || strcmp(PathInit(2:3),[':',PathBar]) || strcmp(PathInit(1:2),['.',PathBar])
                % IN WINDOWS
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
        function fillPath(IPSIM, Path, MemberVars, Track2Error)
            global PathBar
            [ValidPath, ErrorMessage] = IPSIM.isValidPath(Path); % can be absolut or relative path
            if ~ValidPath
                error(sprintf([Track2Error,' -> Path -> ',ErrorMessage]));
            end
            
            % if path is relative error, PENDING: build a full path
            if strcmp(Path(1:2),['.',PathBar])
                error(sprintf([Track2Error,' -> Path -> relative path not allowed']));
                %Path = [IPSIM.InputPath, Path(3:end)];
            end
            
            % Check if path exists
            PathExist = exist(Path,'dir');
            if PathExist~=7
                error(sprintf([Track2Error,' -> Path -> path does not exist']));
            else
                % this allows for single elements to be string or cell
                if iscell(MemberVars), nMemberVars=length(MemberVars); else nMemberVars = 1; MemberVars = {MemberVars}; end
                for i=1:nMemberVars
                    eval(['IPSIM.',MemberVars{i},' = ''',Path,''';']);
                end
            end
        end
        function IsValid = isValidFileExt(IPSIM, File, ValidExtensions)
            [FileName , Extension] = getFilenameAndExt(File);
            if sum(strcmpi(Extension, ValidExtensions)) > 0
                IsValid = 1;
            else
                IsValid = 0;
            end
        end
        function fillFile(IPSIM, File, Path, ValidExtensions, MemberVars, Track2Error)

            % is a valid extension?
            ValidFileExt = IPSIM.isValidFileExt(File, ValidExtensions);
            if ~ValidFileExt
                nExtensions = length(ValidExtensions);
                ValuesString = [];
                for i=1:nExtensions
                    ValuesString = [ValuesString, ' ',ValidExtensions{i}];
                end
                error(sprintf([Track2Error,' -> File ',File,'\n','   valid extensions are:',ValuesString]));
            end
            
            % file exists?
            FileExist = exist([Path,File],'file');
            if FileExist ~=2
                PrintPath = getPrintPath(Path);
                error(sprintf([Track2Error,' -> File: ',File,'\n', '   can not be found in path:  ',PrintPath]));
            else
                nMemberVars = size(MemberVars,1);
                if nMemberVars == 1 % this allows for single elements to be string or cell
                    if ~iscell(MemberVars), MemberVars = {MemberVars}; end
                end
                for i=1:nMemberVars
                    eval(['IPSIM.',MemberVars{i},' = ''',File,''';']);
                end
            end
        end
        function readxml(IPSIM,PathXML,FileXML)
            global PathBar PathBar2
            if ~strcmpi(PathXML(end),PathBar)
                PathXML = [PathXML,PathBar];
                dispWarning('Path for file .SID does not finish with ',PathBar2,' Error corrected automatically',0)
            end
            PathXML = checkAndCorrectFileAndPath(PathXML, FileXML, {'sid'}); % define between {} the accepted extensions separated by ;
            IPSIM.InputPath = PathXML;
            IPSIM.InputFile = FileXML;
            IPSIM.InterfaceDOM = xmlread([PathXML,FileXML]);
        end
        function Done = fillAttribute(IPSIM, AttributeContent, MemberVar, ElementValueList, InternalValueList)
            nValues = length(ElementValueList);
            Done = 0; ValuesString = [];
            for i=1:nValues
                if strcmpi(ElementValueList{i}, AttributeContent)
                    eval(['IPSIM.',MemberVar,' = InternalValueList{i};']);
                    Done = 1;
                end
                ValuesString = [ValuesString, ' ',ElementValueList{i}];
            end
        end
        function RootElement = getRootElement(IPSIM)
            RootElement = IPSIM.InterfaceDOM.getElementsByTagName('SimulationDefinition');
            if RootElement.getLength == 0
                error('Wrong syntax for root element');
            elseif RootElement.getLength > 1
                error('Only one element "SimulationDefinition" is allowed');
            end
            % read attribute Type
            Type = char(RootElement.item(0).getAttribute('Type'));
            if isempty(Type) % use default value or send error
                IPSIM.Settings.Type = IPSIM.DefaultSettings.Type;
            else
                Done = IPSIM.fillAttribute(Type, 'Settings.Type', {'KS','DS'}, {'KS','DS'});                
                if ~Done
                    error(['File ',IPSIM.InputFile,' -> SimulationDefinition -> Type -> \n',...
                           '       ',Type,' is not a valid value. Valid values are: KS or DS']);
                end                
            end
            % read attribute InfoDisplay
            InfoDisplay = char(RootElement.item(0).getAttribute('InfoDisplay'));
            if isempty(InfoDisplay) % use default value or send error
                IPSIM.Settings.InfoDisplay = IPSIM.DefaultSettings.InfoDisplay;
            else
                Done = IPSIM.fillAttribute(InfoDisplay, 'Settings.InfoDisplay', {'Minimum','Standard'}, {0,1});                
                if ~Done
                    error(['File ',IPSIM.InputFile,' -> SimulationDefinition -> InfoDisplay -> \n',...
                           '       ',Type,' is not a valid value. Valid values are: Minimum or Standard']);
                end                
            end
        end
        function parseAdditionalResults(IPSIM, RootElement)
            % AdditionalResults element. Includes elements: Compamm
            AdditResultsElement = RootElement.item(0).getElementsByTagName('AdditionalResults');
            nElement = IPSIM.checkMaxOneElement(AdditResultsElement);
            
            if nElement == 1 
                
                % COMPAMM element. Includes elements: Playback, Position, Velocity, Sensor, etc.
                CompElement = AdditResultsElement.item(0).getElementsByTagName('Compamm');
                nCompElement = IPSIM.checkMaxOneElement(CompElement);
                
                % read attribute Playback
                Playback = char(CompElement.item(0).getAttribute('Playback'));
                if isempty(Playback) % use default value or send error
                    IPSIM.Settings.Results.CompPlayback = IPSIM.DefaultSettings.Results.CompPlayback;
                else
                    Done = IPSIM.fillAttribute(Playback, 'Settings.Results.CompPlayback', {'Yes','No'}, {1,0});
                    if ~Done
                        error(['AdditionalResults -> Compamm -> Playback -> \n',...
                            '       ',Playback,' is not a valid value. Valid values are: Yes or No']);
                    end
                end
                
                % read attribute Position
                Position = char(CompElement.item(0).getAttribute('Position'));
                if isempty(Position) % use default value or send error
                    IPSIM.Settings.Results.Position = IPSIM.DefaultSettings.Results.Position;
                else
                    Done = IPSIM.fillAttribute(Position, 'Settings.Results.Position', {'Yes','No'}, {1,0});
                    if ~Done
                        error(['AdditionalResults -> Compamm -> Position -> \n',...
                            '       ',Position,' is not a valid value. Valid values are: Yes or No']);
                    end
                end

                % read attribute Velocity
                Velocity = char(CompElement.item(0).getAttribute('Velocity'));
                if isempty(Velocity) % use default value or send error
                    IPSIM.Settings.Results.Velocity = IPSIM.DefaultSettings.Results.Velocity;
                else
                    Done = IPSIM.fillAttribute(Velocity, 'Settings.Results.Velocity', {'Yes','No'}, {1,0});
                    if ~Done
                        error(['AdditionalResults -> Compamm -> Velocity -> \n',...
                            '       ',Velocity,' is not a valid value. Valid values are: Yes or No']);
                    end
                end

                % read attribute Acceleration
                Acceleration = char(CompElement.item(0).getAttribute('Acceleration'));
                if isempty(Acceleration) % use default value or send error
                    IPSIM.Settings.Results.Acceleration = IPSIM.DefaultSettings.Results.Acceleration;
                else
                    Done = IPSIM.fillAttribute(Acceleration, 'Settings.Results.Acceleration', {'Yes','No'}, {1,0});
                    if ~Done
                        error(['AdditionalResults -> Compamm -> Acceleration -> \n',...
                            '       ',Acceleration,' is not a valid value. Valid values are: Yes or No']);
                    end
                end

                % read attribute EndEffectorTrajectory
                EndEffTraj = char(CompElement.item(0).getAttribute('EndEffectorTrajectory'));
                if isempty(EndEffTraj) % use default value or send error
                    IPSIM.Settings.Results.EndEffTraj = IPSIM.DefaultSettings.Results.EndEffTraj;
                else
                    Done = IPSIM.fillAttribute(EndEffTraj, 'Settings.Results.EndEffTraj', {'Yes','No'}, {1,0});
                    if ~Done
                        error(['AdditionalResults -> Compamm -> EndEffectorTrajectory -> \n',...
                            '       ',EndEffTraj,' is not a valid value. Valid values are: Yes or No']);
                    end
                end

                % read attribute JointEffort
                JointEffort = char(CompElement.item(0).getAttribute('JointEffort'));
                if isempty(JointEffort) % use default value or send error
                    IPSIM.Settings.Results.JointEffort = IPSIM.DefaultSettings.Results.JointEffort;
                else
                    Done = IPSIM.fillAttribute(JointEffort, 'Settings.Results.JointEffort', {'Yes','No'}, {1,0});
                    if ~Done
                        error(['AdditionalResults -> Compamm -> JointEffort -> \n',...
                            '       ',JointEffort,' is not a valid value. Valid values are: Yes or No']);
                    end
                end
                
                % Check Settings incompatibility
                if strcmpi(IPSIM.Settings.Type,'KS') && (IPSIM.Settings.Results.JointEffort == 1)
                    error('It is not possible to generate results for JointEffort. Set Type to "DS"')
                end                
            end            
        end % function
        function parseSimEnvironment(IPSIM, RootElement)
            
            % --------------------------------------------------------------------------
            % SimEnvironment element. Includes elements: Seat, Pedal & OtherObjects
            % --------------------------------------------------------------------------
            SimEnvElement = RootElement.item(0).getElementsByTagName('SimEnvironment');
            IPSIM.checkJustOneElement(SimEnvElement);

            % read attribute Name
            EnvName = char(SimEnvElement.item(0).getAttribute('Name'));
            if isempty(EnvName) 
                error('SimEnvironment -> Name -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Name = EnvName;
            end
            
            % read attribute GraphicsPath
            GraphicsPath = char(SimEnvElement.item(0).getAttribute('GraphicsPath'));
            if isempty(GraphicsPath)
                error('SimEnvironment -> GraphicsPath -> attribute compulsory but not provided');
            else
                IPSIM.fillPath(GraphicsPath, 'Experiment.SimEnvironment.GraphicsPath', 'SimEnvironment -> GraphicPath');
            end
            
            
            % --------------------------------------------------------------------------
            % Seat element
            % --------------------------------------------------------------------------
            SeatElement = SimEnvElement.item(0).getElementsByTagName('Seat');
            IPSIM.checkJustOneElement(SeatElement);

            % read attribute Unit
            Unit = char(SeatElement.item(0).getAttribute('Unit'));
            if isempty(Unit) 
                error('SimEnvironment -> Seat -> Unit -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Seat.Unit = Unit;
            end            
            
            % read attribute HPoint
            HPoint = str2num(char(SeatElement.item(0).getAttribute('H-point')));
            if isempty(HPoint) 
                error('SimEnvironment -> Seat -> HPoint -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Seat.HPoint = HPoint;
            end

            % read attribute SEA_BL
            SEA_BL = str2num(char(SeatElement.item(0).getAttribute('SEA_BL')));
            if isempty(SEA_BL) 
                error('SimEnvironment -> Seat -> SEA_BL -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Seat.SEA_BL = SEA_BL;
            end            
            % read attribute SEA_FL
            SEA_FL = str2num(char(SeatElement.item(0).getAttribute('SEA_FL')));
            if isempty(SEA_FL) 
                error('SimEnvironment -> Seat -> SEA_FL -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Seat.SEA_FL = SEA_FL;
            end
            % read attribute SEA_BR
            SEA_BR = str2num(char(SeatElement.item(0).getAttribute('SEA_BR')));
            if isempty(SEA_BR) 
                error('SimEnvironment -> Seat -> SEA_BR -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Seat.SEA_BR = SEA_BR;
            end            
            % read attribute SEA_BU
            SEA_BU = str2num(char(SeatElement.item(0).getAttribute('SEA_BU')));
            if isempty(SEA_BU) 
                error('SimEnvironment -> Seat -> SEA_BU -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Seat.SEA_BU = SEA_BU;
            end
            
            % read attribute GraphicFile
            GraphicFile = char(SeatElement.item(0).getAttribute('GraphicFile'));
            if isempty(GraphicFile) 
                error('SimEnvironment -> Seat -> GraphicFile -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Seat.GraphicFile = GraphicFile;
            end            
            
            
            % --------------------------------------------------------------------------
            % Pedal element
            % --------------------------------------------------------------------------
            PedalElement = SimEnvElement.item(0).getElementsByTagName('Pedal');
            IPSIM.checkJustOneElement(PedalElement);
            
            % read attribute Unit
            Unit = char(PedalElement.item(0).getAttribute('Unit'));
            if isempty(Unit) 
                error('SimEnvironment -> Pedal -> Unit -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Pedal.Unit = Unit;
            end            

            % read attribute UnpressedPos
            UnpressedPos = str2num(char(PedalElement.item(0).getAttribute('UnpressedPos')));
            if isempty(UnpressedPos) 
                error('SimEnvironment -> Pedal -> UnpressedPos -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Pedal.UnpressedPos = UnpressedPos;
            end            
            
            % read attribute TravelLength
            TravelLength = str2double(PedalElement.item(0).getAttribute('TravelLength'));
            if isempty(TravelLength) 
                error('SimEnvironment -> Pedal -> TravelLength -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Pedal.TravelLength = TravelLength;
            end            
            
            % read attribute TravelAngle
            TravelAngle = str2double(PedalElement.item(0).getAttribute('TravelAngle'));
            if isempty(TravelAngle) 
                error('SimEnvironment -> Pedal -> TravelAngle -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Pedal.TravelAngle = TravelAngle;
            end            
            
            % read attribute LeverArm
            LeverArm = str2double(PedalElement.item(0).getAttribute('LeverArm'));
            if isempty(LeverArm) 
                error('SimEnvironment -> Pedal -> LeverArm -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Pedal.LeverArm = LeverArm;
            end            
            
            % read attribute RotAxisVertexL
            RotAxisVertexL = str2num(char(PedalElement.item(0).getAttribute('RotAxisVertexL')));
            if isempty(RotAxisVertexL) 
                error('SimEnvironment -> Pedal -> RotAxisVertexL -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Pedal.RotAxisVertexL = RotAxisVertexL;
            end            
            
            % read attribute RotAxisVertexR
            RotAxisVertexR = str2num(char(PedalElement.item(0).getAttribute('RotAxisVertexR')));
            if isempty(RotAxisVertexR) 
                error('SimEnvironment -> Pedal -> RotAxisVertexR -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Pedal.RotAxisVertexR = RotAxisVertexR;
            end            
            
            % read attribute GraphicFile
            GraphicFile = char(PedalElement.item(0).getAttribute('GraphicFile'));
            if isempty(GraphicFile) 
                error('SimEnvironment -> Pedal -> GraphicFile -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Pedal.GraphicFile = GraphicFile;
            end            
            
            % --------------------------------------------------------------------------
            % Stiffness element
            % --------------------------------------------------------------------------
            StiffElement = PedalElement.item(0).getElementsByTagName('Stiffness');
            IPSIM.checkJustOneElement(StiffElement);
            
            % read attribute Path
            Path = char(StiffElement.item(0).getAttribute('Path'));
            if isempty(Path)
                error('SimEnvironment -> Pedal -> Stiffness -> Path -> attribute compulsory but not provided');
            else
                IPSIM.fillPath(Path, 'Experiment.SimEnvironment.Pedal.Stiffness.Path', 'SimEnvironment -> Pedal -> Stiffness');
            end
            
            % read attribute File
            File = char(StiffElement.item(0).getAttribute('File'));
            if isempty(File) 
                error('SimEnvironment -> Pedal -> Stiffness -> File -> attribute compulsory but not provided');
            else
                IPSIM.fillFile(File, Path, {'xml'}, 'Experiment.SimEnvironment.Pedal.Stiffness.File', ...
                               'SimEnvironment -> Pedal -> Stiffness');
            end            
            
            % --------------------------------------------------------------------------
            % OtherObjects element
            % --------------------------------------------------------------------------
            OtherObjElement = SimEnvElement.item(0).getElementsByTagName('OtherObjects');
            % read Graphic Files
            GraphFileElement = OtherObjElement.item(0).getElementsByTagName('GraphicFile');            
            nGraphicFiles = GraphFileElement.getLength;
            
            if nGraphicFiles >= 1
                OtherObjects = cell(nGraphicFiles,1);
                for j=0:nGraphicFiles-1                    
                    GraphicFile_j = char(GraphFileElement.item(j).getTextContent);                    
                    % Does it have the right extension?
                      % TODO
                    % Does file exist in the right path?
                      % TODO
                    % Stote file in motion files cell
                    OtherObjects{j+1} = GraphicFile_j;
                end
                IPSIM.Experiment.SimEnvironment.OtherObjects = OtherObjects;
            else
                IPSIM.Experiment.SimEnvironment.OtherObjects = {};
            end
        end
        function parseSimMotions(IPSIM, RootElement)
            
            % SimMotions element. Includes elements: ActiveModel, ComplementaryPosture, 
            %                                        SimEquations, SimPeriod & SimMotion ======
            SimMotionsElement = RootElement.item(0).getElementsByTagName('SimMotions');
            IPSIM.checkJustOneElement(SimMotionsElement);
                        
            % --------------------------------------------------------------------------
            % read element ActiveModel
            % --------------------------------------------------------------------------
            ActiveModelElement = SimMotionsElement.item(0).getElementsByTagName('ActiveModel');
            IPSIM.checkJustOneElement(ActiveModelElement);
            
            % read attribute Path
            Path = char(ActiveModelElement.item(0).getAttribute('Path'));
            if isempty(Path)
                error('SimMotions -> ActiveModel -> Path -> attribute compulsory but not provided');
            else
                IPSIM.fillPath(Path, 'Model.Path', 'SimMotions -> ActiveModel');
            end
            
            % read attribute File
            File = char(ActiveModelElement.item(0).getAttribute('File'));
            if isempty(File) 
                error('SimMotions -> ActiveModel -> File -> attribute compulsory but not provided');
            else
                IPSIM.fillFile(File, Path, {'xml';'mat'}, 'Model.File', 'SimMotions -> ActiveModel');
            end            
            
            % --------------------------------------------------------------------------
            % read element ComplementaryPosture
            % --------------------------------------------------------------------------
            CompPostureElement = SimMotionsElement.item(0).getElementsByTagName('ComplementaryPosture');
            IPSIM.checkJustOneElement(CompPostureElement);
            
            % read attribute Path
            Path = char(CompPostureElement.item(0).getAttribute('Path'));
            if isempty(Path)
                error('SimMotions -> ComplementaryPosture -> Path -> attribute compulsory but not provided');
            else
                IPSIM.fillPath(Path, 'Experiment.CompPosture.Path', 'SimMotions -> ComplementaryPosture');
            end
            
            % read attribute File
            File = char(CompPostureElement.item(0).getAttribute('File'));
            if isempty(File) 
                error('SimMotions -> ComplementaryPosture -> File -> attribute compulsory but not provided');
            else
                IPSIM.fillFile(File, Path, {'kfr';'xml'}, 'Experiment.CompPosture.File', 'SimMotions -> ComplementaryPosture');
            end            

            % --------------------------------------------------------------------------
            % read element SimEquations
            % --------------------------------------------------------------------------
            SimEquationsElement = SimMotionsElement.item(0).getElementsByTagName('SimEquations');
            IPSIM.checkJustOneElement(SimEquationsElement);
            
            % read attribute Path
            Path = char(SimEquationsElement.item(0).getAttribute('Path'));
            if isempty(Path)
                error('SimMotions -> SimEquations -> Path -> attribute compulsory but not provided');
            else
                IPSIM.fillPath(Path, 'Experiment.EquationsPath', 'SimMotions -> SimEquations');
            end

            % --------------------------------------------------------------------------
            % read element SimPeriod
            % --------------------------------------------------------------------------
            SimPeriodElement = SimMotionsElement.item(0).getElementsByTagName('SimPeriod');
            IPSIM.checkJustOneElement(SimPeriodElement);
            
            % read element SimStart -----------------------------------
            % read element SimFinish -----------------------------------
            
            % read element KeyFrame -----------------------------------
            KeyFrameElement = SimPeriodElement.item(0).getElementsByTagName('KeyFrame');
            IPSIM.checkJustOneElement(KeyFrameElement);
            % read attribute Name
            Name = char(KeyFrameElement.item(0).getAttribute('Name'));
            if isempty(Name)
                error('SimMotions -> SimPeriod -> KeyFrame -> Name -> attribute compulsory but not provided');
            else
                IPSIM.Experiment.SimEnvironment.Name = EnvName;
                
            end

            
            
            
            
            
            
            
            
            
        end
     end % method
    
end % class
