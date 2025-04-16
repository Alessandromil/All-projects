function InterfaceSubPar(Path, FileXML)
global PathBar

% Interface function that reads xml file defining which subject 
% parameters must be generated from palpation measurements


% read XML file
[ExpPath, Subjects, SubParSettings] = parseSubParXML(Path, FileXML);

% Para poder generar PSP hay que suministrar dos Paths.
% Estos paths habra que darlos en xml
% SubParSettings.PAM_OutputPath = '';
% SubParSettings.PAM_DeskPath = '';


% Write ALM from palpation data. Two options to call the function
% OPTION 1)Path and Subjects are selected by the user in a dialog window 
% [ExpPath, Subjects] = writeALM();
% OPTION 2)Path and Subjects are provided by user as inputs
if SubParSettings.Results.ALM == 1
    writeALM(ExpPath, Subjects);
    
elseif SubParSettings.Results.ALM == 0
    
    nSubjects = length(Subjects);
    ErrMessage = sprintf('ALM file(s) for the following subject(s) not found:\n');
    ErrorYES = 0;
    
    for i=1:nSubjects
        SubjectCode = Subjects(i).SubPar;
        SubParFile = [SubjectCode,'.alm'];
        SubParPath = [ExpPath,'Subjects',PathBar,SubjectCode,PathBar];
        if exist([SubParPath,SubParFile],'file') ~=2
            ErrMessage = [ErrMessage, sprintf(['  File ',SubParFile,' not found in ',getPrintPath(SubParPath),'\n'])];
            ErrorYES = 1;
        end
    end
    
    if ErrorYES
        ErrMessage = [ErrMessage, 'Set Results_SubPar element ALM to YES or add an ALM file to subject path'];
        error(ErrMessage);
    end
end


% Write Subject Parameter files:
%   file XML (intermediate subject parameter file)
%   file RSP (Ramsis Subject Parameters file) 
%   file PSP (PAM Subject Parameters file) 
writeSubPar(ExpPath, Subjects, SubParSettings)

