function [Markers, Frequency, IndexAllVisible, MarkerNames, MarkerCoords] = readC3D(InputPath, InputFile)

% READFILEC3D read a motion capture file in .C3D format using C3Dserver (Shriners Hospital, 2006)
%   [Frequency,MarkerNames, MarkerCoords] = readFileC3D()
%   Outputs:
%     + MarkerNames is a cell (nMarkers x 1) with the names of the skin-markers. Each
%       element of the cell is a string.
%     + MarkerCoords is a double array(nFrames x 3*nMarkers) with the 3D coordinates of
%       skin-markers. The coordinates x,y,z of skin-marker i are in position
%       MarkerCoords(3*i-2,3*i-1,3*i)
%     + Frequency is a double array (1x1) with the frequency in Hertz at which the motion
%       was recorded.
%     + Markers is a structure with nMarkers fields. The name of each field is the name of
%       the marker. For example, for marker 'M4' it coordinates are Markers.M4, this is an
%       array (nFramesx3) with M4 x,y,z coordinates for each frame

% display info
disp(['Reading file ',InputFile,' ...']);

checkFileAndPath(InputPath,InputFile,{'C3D'})

% Activate C3Dserver as a COM object
SERVER = c3dserver();

% Load file C3D in C3Dserver
openc3d(SERVER, 0, [InputPath, InputFile]);

% Read marker coordinates
Markers = get3dtargets(SERVER);

% Read Video Frame Rate (VFR)
% Frequency = getVFR(SERVER); En el interior de la función
VFR_Index = SERVER.GetParameterIndex('POINT', 'RATE');
Frequency = SERVER.GetParameterValue(VFR_Index, 0);

% close file and C3Dserver
closec3d(SERVER);

% Units of marker coordinates
Units = Markers.units;
% Remove field 'units', only markers are left
Markers = rmfield(Markers, 'units');

% Extraer nombres de los markers
MarkerNames = fieldnames(Markers);
nMarkers    = length(MarkerNames);

% Convert marker coordinates units to metres if necessary.
if strcmpi(Units,'mm')
    for i=1:nMarkers
        Markers.(MarkerNames{i}) = (1/1000) * Markers.(MarkerNames{i});
    end   
    
elseif ~strcmpi(Units,'m')
    error('Marker coordinate units are neither metres [m] nor milimetres [mm]. Appropriated conversion should be implemented');
end

% Crear matriz MarkerCoords array(nFrames x 3*nMarkers)
MarkerCoords = [];
for i=1:nMarkers
    MarkerCoords = [MarkerCoords, Markers.(MarkerNames{i})];
end

% Cut out of trajectories. In the 1st frame must be at least one marker (not all may be lost)
% We also calculate the first frame where they appear all markers
Index = 0;
IndexAllVisible = 0;
nFrames  = size(MarkerCoords,1);
i = 1;
while IndexAllVisible == 0  &&  i<=nFrames
    if ((sum(isnan(MarkerCoords(i,:))) < 3*nMarkers - 3) && (Index==0)) % At least 1 visible Marker  
        Index = i;
    elseif (sum(isnan(MarkerCoords(i,:))) == 0)                         % All markers visibles
        IndexAllVisible = i;      
    else
        i = i + 1;
    end
end
if IndexAllVisible == 0
    warning('Not found a frame where all the markers are visible');
end
MarkerCoords = MarkerCoords(Index:end,:);
MarkerCoords = double(MarkerCoords);
for i=1:nMarkers
    Markers.(MarkerNames{i}) = double(Markers.(MarkerNames{i})(Index:end,:));    
end

