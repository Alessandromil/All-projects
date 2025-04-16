function InterfaceMotSim(Path, FileXML)

% Define parameters
% Two options: define them in XML file or in M-File
% Here only XML file option because this is for DHErgo Demonstrator interface
[Models,Experiment,Settings] = parseMotSimXML(Path, FileXML);

% Do motion simulationreconstruction
DynamicSimulation(Models,Experiment,Settings);