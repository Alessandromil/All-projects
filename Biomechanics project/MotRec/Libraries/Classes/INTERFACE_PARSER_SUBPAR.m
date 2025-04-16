classdef INTERFACE_PARSER_SUBPAR < handle
    % INTERFACE_MOTREC_PARSE_XML reads input data for motion reconstruction
    % in XML file format
    
    properties
        ExpPath             % String with experiment path
        Subjects            % Struct with subject data
                            %   Subjects(i).SubPar -> subject code/file
                            %   Subjects(i).Dist_FM1_FM2 -> subject specific par
        Settings            % Struct with setting for motion reconstruction
                            %   Settings.Type -> {'SUBPAR'}
                            %   Settings.Display -> {'Minimum','Standard'}
                            %   Settings.Results.ALM
                            %   Settings.Results.XML
                            %   Settings.Results.RSP
                            %   Settings.Results.PSP
        InterfaceDOM        % Document Object Model node representing the parsed XML file
        InterfaceFilePath   % String with path
        InterfaceFile       % String with file name and extension of the SubPar Interface Definition file
        DefaultSettings     % Struct with default settings
    end
    
    methods
        function ISP = INTERFACE_PARSER_SUBPAR()
            
            ISP.DefaultSettings.Display = 1;  
            
            % Default values for Results. Yes is 1, No is 0
            ISP.DefaultSettings.Results.ALM = 0; % Turns on/off generation of ALM file
            ISP.DefaultSettings.Results.XML = 0; % Turns on/off generation of XML file
            ISP.DefaultSettings.Results.RSP = 0; % Turns on/off generation of RSP file
            ISP.DefaultSettings.Results.PSP = 0; % Turns on/off generation of PSP file

        end
        function checkJustOneElement2(ISP, Element, ElementName)
            if Element.getLength == 0
                error(['File ',ISP.InterfaceFile,' -> element "',ElementName,'" not defined']);
            elseif Element.getLength > 1
                error(['File ',ISP.InterfaceFile,' -> only ONE element "',ElementName,'" can be defined']);
            end
        end
        function checkJustOneElement(ISP, Element, TagPath, TagName)
            if Element.getLength == 0
                error(['File ',ISP.InterfaceFile,' -> ',TagPath,' one element "',TagName,'" must be defined']);
            elseif Element.getLength > 1
                error(['File ',ISP.InterfaceFile,' -> ',TagPath,' only ONE element "',TagName,'" can be defined']);
            end
        end
        function nElem = checkMaxOneElement(ISP, Element, TagPath, TagName)
            if Element.getLength == 0
                nElem = 0;
            elseif Element.getLength == 1
                nElem = 1;
            elseif Element.getLength > 1
                nElem = 2;
                error(['File ',ISP.InterfaceFile,' -> ',TagPath,' only ONE element "',TagName,'" can be defined']);
            end
        end
        function fillSingleElementValues(ISP, Element, MemberVar, ElementValueList, InternalValueList)
            ElementValue = char(Element.item(0).getTextContent);
            ElementName  = char(Element.item(0).getTagName);
            nValues = length(ElementValueList);
            Done = 0; ValuesString = [];
            for i=1:nValues
                if strcmpi(ElementValueList{i}, ElementValue)
                    eval(['ISP.',MemberVar,' = InternalValueList{i};']);
                    Done = 1;
                end
                ValuesString = [ValuesString, ' ',ElementValueList{i}];
            end
            if ~Done
                error(['File ',ISP.InterfaceFile,' -> ',ElementName,' = ',ElementValue,' -> valid values are:',ValuesString]);
            end
        end
        function [IsValid, MessageError] = isValidPath(ISP, Path)
            global PathBar PathBar2
            PathInit = Path(1:3);
            if ~strcmpi(Path(end),PathBar) % check is Path (full or relative) ends with the right character
                IsValid = 0; % invalid path
                MessageError = ['path must finish with "',PathBar2,'".'];
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
        function fillPath(ISP, PathElement, MemberVars, ParentElementName)
            global PathBar
            Path = char(PathElement.item(0).getTextContent);
            [ValidPath, ErrorMessage] = ISP.isValidPath(Path); % can be absolut or relative path
            if ~ValidPath
                error(['File ',ISP.InterfaceFile,' -> ',ParentElementName,' -> Path -> ',ErrorMessage]);
            end
            
            % if path is relative build a full path
            if strcmp(Path(1:2), ['.',PathBar])
                Path = [ISP.InterfaceFilePath, Path(3:end)];
            end
            
            % Check if path exists
            PathExist = exist(Path,'dir');
            if PathExist~=7
                error(['File ',ISP.InterfaceFile,' -> ',ParentElementName,' -> Path -> path does not exist']);
            else
                % this allows for single elements to be string or cell
                if iscell(MemberVars), nMemberVars=length(MemberVars); else nMemberVars = 1; MemberVars = {MemberVars}; end
                %                 if nMemberVars == 1 % this allows for single elements to be string or cell
                %                     if ~iscell(MemberVars), MemberVars = {MemberVars}; end
                %                 end
                for i=1:nMemberVars
                    eval(['ISP.',MemberVars{i},' = ''',Path,''';']);
                end
            end
        end
        function IsValid = isValidFileExt(ISP, File, ValidExtensions)
            [FileName , Extension] = getFilenameAndExt(File);
            if sum(strcmpi(Extension, ValidExtensions)) > 0
                IsValid = 1;
            else
                IsValid = 0;
            end
        end
        function fillFile(ISP, FileElement, Path, MemberVars, ParentElementName, ValidExtensions)
            
            File = char(FileElement.item(0).getTextContent);
            ElementName  = char(FileElement.item(0).getTagName);
            
            % is a valid extension?
            ValidFileExt = ISP.isValidFileExt(File, ValidExtensions);
            if ~ValidFileExt
                nExtensions = length(ValidExtensions);
                ValuesString = [];
                for i=1:nExtensions
                    ValuesString = [ValuesString, ' ',ValidExtensions{i}];
                end
                error(sprintf(['File ',ISP.InterfaceFile,' -> ',ParentElementName,' -> ',...
                    ElementName,' -> file ',File,'\n', '   valid extensions are:',ValuesString]));
            end
            
            % file exists?
            FileExist = exist([Path,File],'file');
            if FileExist ~=2
                PrintPath = getPrintPath(Path);
                error(sprintf(['File ',ISP.InterfaceFile,' -> ',ParentElementName,' -> ',...
                    ElementName,' -> file: ',File,'\n', '   can not be found in path:\n  ',PrintPath]));
            else
                nMemberVars = length(MemberVars);
                if nMemberVars == 1 % this allows for single elements to be string or cell
                    if ~iscell(MemberVars), MemberVars = {MemberVars}; end
                end
                for i=1:nMemberVars
                    eval(['ISP.',MemberVars{i},' = ''',File,''';']);
                end
            end
        end
        function readxml(ISP, PathXML, FileXML)
            PathXML = checkAndCorrectFileAndPath(PathXML, FileXML, {'xml'});
            ISP.InterfaceFilePath = PathXML;
            ISP.InterfaceFile = FileXML;
            ISP.InterfaceDOM = xmlread([PathXML, FileXML]);
        end
        function RootElement = getRootElement(ISP)
            RootElement = ISP.InterfaceDOM.getElementsByTagName('MotionReconstruction');
            if RootElement.getLength == 0
                error(['File ',ISP.InterfaceFile,' -> wrong syntax for root element']);
            elseif RootElement.getLength > 1
                error(['File ',ISP.InterfaceFile,' -> only one element "MotionReconstruction" is allowed']);
            end
        end
        function parseReconstructionType(ISP, RootElement)
            RecTypeElement = RootElement.item(0).getElementsByTagName('Type');
            ISP.checkJustOneElement(RecTypeElement, 'MotionReconstruction ->', 'Type'); 
            ISP.fillSingleElementValues(RecTypeElement, 'Settings.Type', {'SUBPAR'}, {'SUBPAR'});
        end
        function parseInfoDisplay(ISP, RootElement)
            InfoElement = RootElement.item(0).getElementsByTagName('InfoDisplay');
            nElement = ISP.checkMaxOneElement(InfoElement, 'MotionReconstruction ->', 'InfoDisplay');
            if nElement == 0
                % fill default value
                ISP.Settings.Display = ISP.DefaultSettings.Display;
            elseif nElement == 1
                ISP.fillSingleElementValues(InfoElement, 'Settings.Display', {'Minimum','Standard'}, {0,1});
            end
        end
        function parseResults_SubjectPar(ISP, RootElement)
            % Results_SubPar. Includes elements: ALM, XML, RSP(=EXP) & PSP
            ResultsSubParElement = RootElement.item(0).getElementsByTagName('Results_SubPar');
            nElement = ISP.checkMaxOneElement(ResultsSubParElement, 'MotionReconstruction ->', 'Results_SubPar');

            if nElement == 0
                % Not Subject parameter results asked for, then if
                error('No subject parameters results asked for. Therefore, SUBPAR analysis is not performed');
                
            elseif nElement == 1 % fill settings from SUBPAR analysis               
                % if analysis is SUBPAR (IK or ID) then these results cannot be generated
                if strcmpi(ISP.Settings.Type,'IK') || strcmpi(ISP.Settings.Type,'ID') || strcmpi(ISP.Settings.Type,'STATIC')
                    % error results defined in Results_SubjectPar cannot be generated for analysis IK, ID or STATIC
                    error('Subject paramter files cannot be generated for analysis of Type IK, ID or STATIC');
                end
                
                % ALM settings
                ALM_Element = ResultsSubParElement.item(0).getElementsByTagName('FileALM');
                nElement = ISP.checkMaxOneElement(ALM_Element, 'MotionReconstruction -> Results_SubPar ->', 'FileALM');
                if nElement == 0
                    % fill default value
                    ISP.Settings.Results.ALM = ISP.DefaultSettings.Results.ALM;
                elseif nElement == 1
                    ISP.fillSingleElementValues(ALM_Element, 'Settings.Results.ALM', {'No','Yes'}, {0,1});
                end
                
                % XML settings
                XML_Element = ResultsSubParElement.item(0).getElementsByTagName('FileXML');
                nElement = ISP.checkMaxOneElement(XML_Element, 'MotionReconstruction -> Results_SubPar ->', 'FileXML');
                if nElement == 0
                    % fill default value
                    ISP.Settings.Results.XML = ISP.DefaultSettings.Results.XML;
                elseif nElement == 1
                    ISP.fillSingleElementValues(XML_Element, 'Settings.Results.XML', {'No','Yes'}, {0,1});
                end
                
                % RSP settings
                RSP_Element = ResultsSubParElement.item(0).getElementsByTagName('FileRSP');
                nElement = ISP.checkMaxOneElement(RSP_Element, 'MotionReconstruction -> Results_SubPar ->', 'FileRSP');
                if nElement == 0
                    % fill default value
                    ISP.Settings.Results.RSP = ISP.DefaultSettings.Results.RSP;
                elseif nElement == 1
                    ISP.fillSingleElementValues(RSP_Element, 'Settings.Results.RSP', {'No','Yes'}, {0,1});
                end
                
                % PSP settings
                PSP_Element = ResultsSubParElement.item(0).getElementsByTagName('FilePSP');
                nElement = ISP.checkMaxOneElement(PSP_Element, 'MotionReconstruction -> Results_SubPar ->', 'FilePSP');
                if nElement == 0
                    % fill default value
                    ISP.Settings.Results.PSP = ISP.DefaultSettings.Results.PSP;
                elseif nElement == 1
                    ISP.fillSingleElementValues(PSP_Element, 'Settings.Results.PSP', {'No','Yes'}, {0,1});
                end
            end
            
        end
        function parseExperimentSubPar(ISP, RootElement)
            
            % Experiment element. Includes elements: Path(only 1), GuidedVarsFile(only 1), Subject(several possible)
            ExpElement = RootElement.item(0).getElementsByTagName('Experiment');
            ISP.checkJustOneElement(ExpElement, 'MotionReconstruction ->', 'Experiment');
            
            % ExpPath element
            PathElement = ExpElement.item(0).getElementsByTagName('ExpPath');
            ISP.checkJustOneElement(PathElement, 'MotionReconstruction -> Experiment ->', 'ExpPath');
            ISP.fillPath(PathElement, {'ExpPath'}, 'Experiment');
            
            % Subject element (there can many any number)
            AllSubjectElements = ExpElement.item(0).getElementsByTagName('Subject');
            nSubjects = AllSubjectElements.getLength;
            
            if (nSubjects == 0)  % no subject was defined for estimating his parameters
                error(['File ',ISP.InterfaceFile,' -> Experiment -> \n',...
                    '   at least one element "Subject" must be defined.']);
                
            else % if (nSubjects >= 1)
                
                %ISP.SubjectList = cell(nSubjects,1);
                % FOR each element "Subject"
                for i=0:nSubjects-1
                    
                    % A) Parse element SubPar (contains subject code)-----------------------------
                    SubParElement = AllSubjectElements.item(i).getElementsByTagName('SubPar');
                    ISP.checkJustOneElement(SubParElement, ['MotionReconstruction -> Experiment -> Subject ',...
                           num2str(i+1),' ->'], 'SubPar');
                    SubPar = char(SubParElement.item(0).getTextContent);
                    ISP.Subjects(i+1).SubPar = SubPar; % it can be a subject code, alm file or xml file
                    
                    % B) Parse element AdditionalPar ------------------------------------------------
                    AddiParElement = AllSubjectElements.item(i).getElementsByTagName('AdditionalPar');
                    nAddiParElement = AddiParElement.getLength;
                        
                    if nAddiParElement == 1 % if nAddiParElement==0 no error, it is OK
                        % here look for all possible addtional parameters.
                        
                        % 1) Additional parameter: Dist_RFM1_RFM2
                        FootParElement = AddiParElement.item(0).getElementsByTagName('Dist_RFM1_RFM2');
                        nElement = ISP.checkMaxOneElement(FootParElement, ['MotionReconstruction -> Experiment -> Subject ',...
                               num2str(i+1),' -> AdditionalPar ->'], 'Dist_RFM1_RFM2');
                        if nElement == 1
                            ISP.Subjects(i+1).Dist_RFM1_RFM2 = str2double(FootParElement.item(0).getTextContent);
                        end
                           
                        % 2) Additional parameter: Dist_LFM1_LFM2
                        FootParElement = AddiParElement.item(0).getElementsByTagName('Dist_LFM1_LFM2');
                        nElement = ISP.checkMaxOneElement(FootParElement, ['MotionReconstruction -> Experiment -> Subject ',...
                               num2str(i+1),' -> AdditionalPar ->'], 'Dist_LFM1_LFM2');
                        if nElement == 1
                            ISP.Subjects(i+1).Dist_LFM1_LFM2 = str2double(FootParElement.item(0).getTextContent);
                        end

                        % x) Additional parameter: XXX to be completed with more parameters
                        
                    elseif nAddiParElement > 1
                        error(['File ',ISP.InterfaceFile,' -> Experiment -> Subject ',num2str(i+1),' -> ',...
                            'only 1 element "AdditionalPar" per subject']);
                    end
                    
                end % for i nSubject
            end % if
            
        end % function parseExperiment()
    end % method
    
end % class

