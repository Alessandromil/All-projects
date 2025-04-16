function writeSubPar(ExperimentPath, Subjects, SubParSettings)

global PathBar

disp(' ');
disp('-----------------------------------------------------------------------------------------')
disp(' Writing Subject parameter files');
disp('-----------------------------------------------------------------------------------------')


% SubjectCodes can be a cell or a struct. If struct it must be translated to cell.
if isstruct(Subjects) % data sent from xml interface file.
    SubjectCodes = {Subjects(:).SubPar};
    SubjectCodes = SubjectCodes';
end

nSubParFiles = length(SubjectCodes);
for i=1:nSubParFiles
    
    try
        SubjectCode = SubjectCodes{i};
        SubjectPath = [ExperimentPath,'Subjects',PathBar,SubjectCode,PathBar];
        disp(['  Subject ',SubjectCode]);
        
        % check file
        AceptedExtension{1,1} = 'alm';
        File = [SubjectCode,'.alm'];
        checkFileAndPath(SubjectPath,File,AceptedExtension);
        
        % read file .alm
        RealHuman = HUMAN_MODEL();
        RealHuman.parseAlmFile(SubjectPath,File);
        RealHuman.setSubAdditPar(Subjects(i)); % fill subject additional parameters (not palpated!!)
        % do calculations
        RealHuman.addModelJoints();
        RealHuman.calcTCS();
        RealHuman.calcTransferedPoint();
        RealHuman.calcLCS();
        RealHuman.calcInertiaPars();
        RealHuman.calcALLCS_T_RamsisLCS();
        
        % write Subject Parameters Intermediate file (.XML)
        disp('    writing Subject Parameters Intermediate file (.XML) ...');
        RealHuman.writeSubjectParsIntermed(SubjectCode, SubjectPath);
        
        % write Ramsis Subject Parameters (.RSP)
        if SubParSettings.Results.RSP == 1
            disp('    writing Ramsis Subject Parameters file (.RSP) ...');
            RealHuman.writeRSP(SubjectCode, SubjectPath);
        end
        
        % write PAM Subject Parameters (.PSP)
        if SubParSettings.Results.PSP == 1
            PAM_OutputPath = ''; % to be define in interface xml file
            PAM_DeskPath = ''; % to be define in interface xml file
            disp('    writing PAM Subject Parameters file (.PSP) ...');
            RealHuman.writePSP(SubjectCode, SubjectPath, ExperimentPath, PAM_OutputPath, PAM_DeskPath);
        end
    catch ME
        FirstLine = sprintf([...
            '    ==================================================================================\n', ...
            '     ERROR while writing subject parameter files:\n']);
        %ErrorMessage = ['   ',ME.message]; % for Matlab 2010b
        ErrorMessage = sprintf(['   ',ME.message]); % for Matlab 2012b, 2014a
        for j=1:length(ME.stack)
            ErrorMessage = [ErrorMessage, sprintf(['\n      ',ME.stack(j).name,'  line:',num2str(ME.stack(j).line)])];
        end
        LastLine = sprintf('\n    ==================================================================================');
        disp([FirstLine,ErrorMessage,LastLine]);
        continue % with next subject
        
    end
end

