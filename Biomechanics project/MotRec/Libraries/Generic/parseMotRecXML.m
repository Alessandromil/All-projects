function [Model,Experiment,Settings] = parseMotRecXML(Path, FileXML)


% Son necesarios dos interfaces: Uno para MotSim y otro MotRec
try    
    % display messages
    LineTop    = sprintf('----------------------------------------------------------------------------------------------------\n');
    HeadMessage = sprintf(' Parser for motion reconstruction');
    LineBottom = sprintf('\n----------------------------------------------------------------------------------------------------\n');
    BodyMessage = sprintf([' Parsing file ',FileXML,'\n','  in folder ',getPrintPath(Path)]);
    FullMessage = [LineTop,HeadMessage,LineBottom,BodyMessage];
    dispInfo(1,FullMessage)
    
    % create class
    InterfaceMotRec = INTERFACE_PARSER_MOTREC();
    
    % read & parse file
    InterfaceMotRec.readxml(Path, FileXML)    
    RootElement = InterfaceMotRec.getRootElement();
    %InterfaceMotRec.parseReconstructionType(RootElement);
    InterfaceMotRec.parseInfoDisplay(RootElement);
    InterfaceMotRec.parseResults_MotRec(RootElement);
    InterfaceMotRec.parseInterpolation(RootElement);
    InterfaceMotRec.parseSmoothing(RootElement);
    InterfaceMotRec.parseModel(RootElement);
    InterfaceMotRec.parseExperiment(RootElement);
    
    % fill output vars
    Settings   = InterfaceMotRec.Settings;
    Model      = InterfaceMotRec.Model;
    Experiment = InterfaceMotRec.Experiment;    
       
catch ME
    
    if strfind(ME.message,'Java')
        % jave error produced by the xml reader
        error(sprintf(['===================================================\n'...
                       ' Java error. See firt red error line\n'...
                       '===================================================\n']))
    else
        FirstLine = sprintf([...
            '\n===================================================================================================\n', ...
            ' ERROR while reading Main file "',FileXML,'":\n']);
        
        %ErrorMessage = ['   ',ME.message]; % for Matlab 2010b
        ErrorMessage = sprintf(['   ',ME.message]); % for Matlab 2012b, 2014a
        
        NErrors=0;
        if isfield(InterfaceMotRec.Settings,'Display') && InterfaceMotRec.Settings.Display == 2
            NErrors = length(ME.stack); % for debugging
        end
        
        for i=1:NErrors
            ErrorMessage = [ErrorMessage, sprintf(['\n      ',ME.stack(i).name,'  line:',num2str(ME.stack(i).line)])];
        end
        LastLine = sprintf('\n===================================================================================================');
        error([FirstLine,ErrorMessage,LastLine]);
        % fill outputs
        Model = []; Experiment= []; Settings = [];
    end
end


