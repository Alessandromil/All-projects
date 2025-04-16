function SmarkCoords = result_CSVnoise(FilePath, FileName, SmarkNames, SmarkCoords, Amax, Wmax, Phasemax, Bmax)


disp(['Writing computer-simulated marker trajectories in file ',FileName,' ...']);

% sizes
[nSamples, nCoords] = size(SmarkCoords);
nSmarks = nCoords/3;

% open file
fid = fopen([FilePath, FileName],'w');

% m 2 mm
SmarkCoords = 1000 * SmarkCoords;


% ------------------------------------------------------
% add random noise
% ------------------------------------------------------
% There are two possible noise models
% * Sergio noise model
%     B + A * sin(w*t + phase)
%     B defines constant error due to errors in model parameter estimation
%     A defines skin movement artifact.
% * Cheze noise model
%     A * sin(w*t + phase)
%
% The actual value of A,w & phase for each coord. of each marker
% is selected randomly in the iterval [0,Amax]; [0,Wmax]; [0,Phasemax] respectively

A  = [];
B  = [];
W  = [];
Ph = [];
if ~isempty(Amax)
    for i=1:nSmarks
        if isempty(B) % Noise model A*sin(w*t + phase)
            A  = [A;  (1/sqrt(3))*Amax(i)*rand(3,1)]; % error interval for marker [-Amax, Amax]
            W  = [W;  Wmax(i)*rand(3,1)];
            Ph = [Ph; Phasemax(i)*rand(3,1)];
        else % Noise model  B + A*sin(w*t + phase)
            A  = [A;  (1/sqrt(3))*Amax(i)*rand(3,1)]; % error interval for marker [-Amax, Amax]
            B  = [B;  (2/sqrt(3))*Bmax(i)*(rand(3,1)-0.5)];  % error interval for marker [-Bmax, Bmax]
            W  = [W;  Wmax(i)*rand(3,1)];
            Ph = [Ph; Phasemax(i)*rand(3,1)];
        end
    end
    
    % Add noise to TRUE marker trajectories
    for i = 1:nSamples
        if isempty(B) % Noise model A*sin(w*t + phase)
            SmarkCoords(i,:) = SmarkCoords(i,:) + (A .* sin(W*0.02*(i-1) + Ph))';
        else % Noise model  B + A*sin(w*t + phase)
            SmarkCoords(i,:) = SmarkCoords(i,:) + B' + (A .* sin(W*0.02*(i-1) + Ph))';
        end
    end
end

% ------------------------------------------------------
% Write .csv file with smark coord.
% ------------------------------------------------------
%if strcmp(RawdataFileExt,'csv') | strcmp(RawdataFileExt,'CSV')
    % header ----------------------------------
    fprintf(fid,'TRAJECTOIRES\r\n');
    fprintf(fid,'50,Hz\r\n');
    SmarkString = [',',SmarkNames{1}];
    for i = 2 : nSmarks
        SmarkString = strcat(SmarkString,[',,,',SmarkNames{i}]);
    end
    fprintf(fid,[SmarkString,'\r\n']);
    FieldString = 'Field #';
    for i = 2 : nSmarks
        FieldString = strcat(FieldString,',X,Y,Z');
    end
    fprintf(fid,[FieldString,'\r\n']);
    % Skin-markers for each sample time -------
    for i = 1: nSamples
        NumberSample = num2str(i);
        SmarkString  = sprintf(',%f',SmarkCoords(i,:));
        fprintf(fid,[NumberSample,SmarkString,'\r\n']);
    end
    fprintf(fid,'\r\n');
 
%else
%    error(['results files with extension "',RawdataFileExt,'" not implemented yet']);
%end

% close file
fclose(fid);
