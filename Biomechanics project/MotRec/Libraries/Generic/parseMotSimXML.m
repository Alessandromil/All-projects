function [Model,Experiment,Settings] = parseMotSimXML(Path, FileXML)



% Son necesarios dos interfaces: Uno para MotSim y otro MotRec
try
    InterfaceMotSim = INTERFACE_PARSER_MOTSIM();
    
    InterfaceMotSim.readxml(Path, FileXML)
    
    RootElement = InterfaceMotSim.getRootElement();
    InterfaceMotSim.parseAdditionalResults(RootElement);
    
    %InterfaceMotSim.parseSolver(RootElement);
    InterfaceMotSim.parseSimEnvironment(RootElement)    
    InterfaceMotSim.parseSimMotions(RootElement);
    
    % fill output vars
    Settings   = InterfaceMotSim.Settings;
    Model      = InterfaceMotSim.Model;
    Experiment = InterfaceMotSim.Experiment;
    
catch ME
    
    FirstLine = sprintf([...
        '\n===================================================================================================\n', ...
        ' ERROR while reading interface file ',FileXML,' for motion simulation:\n']);
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

