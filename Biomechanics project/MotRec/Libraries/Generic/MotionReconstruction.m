function MotionReconstruction(Models,ExperimentData,Settings)

% Settings pending of being programmed
Settings.Results.CentreOfMass = 0;     % Position of centre of mass in Compamm format (*.com)

NModels = size(Models,1);
Solver = 'weightedOTM';
for i=1:NModels

    Model_i = Models(i);
    ModelPath_i = Model_i.Path;
    ModelFile_i = Model_i.File;
    
    % -----------------------------------------------------------------
    % Model creation (if needed)
    % -----------------------------------------------------------------
    
    % Check if model file is .MAT (already created) or .XML (to be created)
    [ModelFileName_i, ModelFileExt_i] = getFilenameAndExt(ModelFile_i);
    if strcmpi(ModelFileExt_i,'xml') % model must be created next
        
        try
            c=tic;            
            if (Settings.Display == 1 || Settings.Display == 2)
                disp('----------------------------------------------------------------------------------------------------');
                disp(' Model Creation');
                disp('----------------------------------------------------------------------------------------------------');
            end
            disp([' Creating model ',ModelFileName_i]);
            disp(['  in folder ',ModelPath_i]);
            Human = HUMAN_MODEL();
            Human.parseXmlModelFile(ModelPath_i,ModelFile_i);
            Human.fillq;
            Human.addMarkerInq();
            % if reconstruction method is SODERKVIST do not create constraints
            if ~strcmpi(Settings.Type,'SODERKVIST')
                Human.mkModelCtrs();
            end
            
            save([ModelPath_i,ModelFileName_i],'Human');
            [H, MI, S] = second2HMS(toc(c));
            printElapsedTime(H, MI, S, 3, '  Model creation time: ');
            % updata model file extension
            Model_i.File = [ModelFileName_i,'.mat'];
            
        catch ME
            
            if strfind(ME.message,'Java')
                % jave error produced by the xml reader
                error(sprintf(['===================================================\n'...
                               ' Java error. See firt red error line\n'...
                               '===================================================\n']))
            else                
                FirstLine = sprintf([...
                    '\n===================================================================================================\n', ...
                    ' ERROR while reading model definition file "',ModelFile_i,'":\n']);

                %ErrorMessage = ['   ',ME.message]; % for Matlab 2010b
                ErrorMessage = sprintf(['   ',ME.message]); % for Matlab 2012b, 2014a 
                
                NErrors=0;
                if Settings.Display == 2 
                    NErrors = length(ME.stack); % for debugging
                end
                
                for j=1:NErrors 
                    ErrorMessage = [ErrorMessage, sprintf(['\n      ',ME.stack(j).name,'  line:',num2str(ME.stack(j).line)])];
                end
                LastLine = sprintf('\n===================================================================================================');
                error([FirstLine,ErrorMessage,LastLine]);
            end
        end
    end
    
    % -----------------------------------------------------------------
    % Motion Reconstruction
    % -----------------------------------------------------------------
    Experiment = EXPERIMENT(Model_i, ExperimentData, Solver, Settings);
    Experiment.doMotRec();
    ExperimentData = {};

end
