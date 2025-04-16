function InverseKinematic(Models,Experiment,Settings)

global PathBar

Settings.ResultsOpt = 0;
Settings.DHErgoOutput = 0;

NModels = size(Models,1);
Solver = 'weightedOTM';
if isempty(Models)
    error('You must chose a model')
end
for i=1:NModels
    Model = Models(i);
    
    % 1) Check file model definition (.mat or .xml)
    AceptedExtension{1,1} = 'mat'; AceptedExtension{2,1} = 'xml';
    checkFileAndPath(Model.Path,Model.File,AceptedExtension);

    % 2) Check file model guided (.m)
    AceptedExtension{1,1} = 'm';
    checkFileAndPath([Experiment.Path,'AdditionalEquations',PathBar],[Model.Guided],AceptedExtension);
    [Model.Guided , ModelExtension] = getFilenameAndExt(Model.Guided);

    % If the program does not enter into any 'if' means that 
    % strutc 'Experiment' contains a field "Name" defined by the user
    % 3) Check Experiment definition
    if isfield(Experiment,'Path')
        if ~strcmp(Experiment.Path(end),PathBar)
            error('The experiment path must be finish with "',PathBar,'".')
        end
        BarPos = findstr(Experiment.Path,PathBar);
        Experiment.Name = [Experiment.Path(BarPos(end-1)+1:BarPos(end)-1)];
    elseif ~isfield(Experiment,'Name')
        error('The struct "Experiment" contains neither field Path nor Name');
    end 
    
    % -----------------------------------------------------------------
    % Motion Reconstruction
    % -----------------------------------------------------------------
    Experiment = EXPERIMENT(Model, Experiment, Solver,Settings);
    Experiment.doMotRec();
end
