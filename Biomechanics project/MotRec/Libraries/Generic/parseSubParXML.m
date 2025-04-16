function [ExpPath, Subjects, Settings] = parseSubParXML(Path, FileXML)

% Son necesarios dos interfaces: Uno para MotSim y otro MotRec
try
    InterfaceSubPar = INTERFACE_PARSER_SUBPAR();
    
    InterfaceSubPar.readxml(Path, FileXML)
    
    RootElement = InterfaceSubPar.getRootElement();
    InterfaceSubPar.parseReconstructionType(RootElement);
    InterfaceSubPar.parseInfoDisplay(RootElement);
    InterfaceSubPar.parseResults_SubjectPar(RootElement);
    InterfaceSubPar.parseExperimentSubPar(RootElement);
    
    % fill output vars
    ExpPath  = InterfaceSubPar.ExpPath;
    Subjects = InterfaceSubPar.Subjects;
    Settings = InterfaceSubPar.Settings;
    
catch ME
    
    FirstLine = sprintf([...
        '\n===================================================================================================\n', ...
        ' ERROR while reading interface file ',FileXML,' for subject parameter estimation:\n']);
    %ErrorMessage = ['   ',ME.message]; % for Matlab 2010b
    ErrorMessage = sprintf(['   ',ME.message]); % for Matlab 2012b, 2014a
    
    for i=1:length(ME.stack)
        ErrorMessage = [ErrorMessage, sprintf(['\n      ',ME.stack(i).name,'  line:',num2str(ME.stack(i).line)])];
    end
    LastLine = sprintf('\n===================================================================================================');
    error([FirstLine,ErrorMessage,LastLine]);
    % fill outputs
    Model = []; Experiment= []; Settings = [];
end


