% This function performs three steps in order to generate computer-simulated motions
%  1) Motion reconstruction of an original motion (*.csv)
%	  This file can contain gaps in marker trajectories of markers missing during the
%     whole motion. From the reconstructed motion, the model-marker trajectories give
%     new trajectories without gaps. With these trajectories a new *.csv file can be created
% 
%  2) Motion reconstruction of the new motion file (*.csv) without missing markers.
%     Este movimiento se considerará como el true motion.
%
%  3) Generation of computer-simulated motions. 
%     Noise has to be added to the motion reconstructed in point 2)


REC_ORIG  = 0; % 1) Reconstruct original motion from csv & generate CSV with model-marker trajectories (no gaps in data)
GEN_noGAP = 0; % 2) Genera motion file using model-markers reconstructed from original motion (no missing markers!)
REC_TRUE  = 0; % 3) Reconstruct true motion from csv without gaps (markers in csv from model-markers)
                  % Reconstructed motion is TRUE motion.
GEN_NOISE = 0; % 4) Create computer-simulated trajectories by adding noise to TRUE trajectory
REC_NOISE = 1; % 5) Reconstruct computer-simulated trajectories 

% ============================================================================================
% ============================================================================================
% 1) Motion reconstruction of an original motion (*.csv)
% ============================================================================================
% ============================================================================================

% --------------------------------------------------------------------
% A) Definition of original motion file, par file & model 
% --------------------------------------------------------------------
% A1) Nombre del fichero CSV original
FileNameEXP = 's21u2m01_ve_swap'; % En este fichero se han correguido un error de swapping en un marker
% A2) Original motion data 
RawdataFiles{1,1} = [pwd,'\Experiments\ElodiePhD_ie5motions\CSV\'];
RawdataFiles{1,2} = {[FileNameEXP,'.csv']};
% A3) Inverse kinematic model and path (comun para todos los movimientos)
ModelIKPath     = [pwd,'\BioModels\RamsisLoc_ieElderly\'];
%ModelIKFileName = 'RamsisLoc_ie_IK';
ModelIKFileName = 'RamsisLoc_ie_LLimbs_IK';

%ModelIKPath     = [pwd,'\BioModels\RamsisLoc_ie2Elderly\'];
%ModelIKFileName = 'RamsisLoc_ie2_IK';
% A4) Subject-specific parameters
CalibFiles{1,1} = [pwd,'\Experiments\ElodiePhD_ie5motions\EXP\'];
CalibFiles{1,2} = 's21.ramsisDimensions.mat';
CalibFiles{1,3} = {[pwd,'\Experiments\ElodiePhD_ie5motions\EXP\'], 's21u2m01_ve.rmr', [], []}; % Initialization
% A5) Path for the original motion results
ResultsPath_ORIG = [pwd,'\Experiments\ElodiePhD_ie5motions\RMR_OTM\'];

% --------------------------------------------------------------------
% B) Reconstruct original motion file
% --------------------------------------------------------------------
if REC_ORIG == 1
    SimType = '12';
    %SimType = '2';
    global TOL; TOL = 1e-5;   % Defines the tolerance used for the simulation
    INFO_TYPE = 1; % INFO_TYPE = 0 only basic info; INFO_TYPE = 1 all info
    %                Log RMR  PB POS PBt0 .pro
    ResultOptions = [ 1;  1;  1;  1;   1;   0];
    % Reconstruct motion
    runModel(SimType, INFO_TYPE, ModelIKFileName, ModelIKPath, CalibFiles, RawdataFiles, ResultsPath_ORIG, 1, ResultOptions);
end


% ============================================================================================
% ============================================================================================
% 2) Generate motion file (*.csv) using model-markers reconstructed
%    from original motion (no missing markers!)
% ============================================================================================
% ============================================================================================
% A1) Nombre del fichero CSV sin missing markers
FileNameEXPnoGaps = [FileNameEXP,'_NoGaps'];
if GEN_noGAP == 1
    % A2) get all markers in original CSV file
    [fs, SmarkNamesCSV, SmarkCoordsCSV]  = readMotionData(RawdataFiles{1,1}, RawdataFiles{1,2}{1});
    % A3) Get variable names in *.POS
    [DataPos, VarNamesPOS] = compread([ResultsPath_ORIG, FileNameEXP, '.pos']);
    % A4) Find which markers of CSV have been included in the model (they are in *.POS)
    [SmarkCoordsCSVinModel, SmarkNamesModel, SmarkCoordsModel] = fillSmarkData(SmarkNamesCSV, SmarkCoordsCSV, VarNamesPOS, DataPos);
    nSmarksModel = length(SmarkNamesModel);
    % A5) Save new CSV file with all markers visible obtained from model-marker trajectories in .pos
    wfileCSV(RawdataFiles{1,1}, [FileNameEXPnoGaps,'.csv'], SmarkNamesModel, 1000*SmarkCoordsModel, 50)
end


% ============================================================================================
% ============================================================================================
% 3) Reconstruct true motion from csv without gaps (markers in csv from model-markers)
% ============================================================================================
% ============================================================================================
% --------------------------------------------------------------------
% A) Definition of motion file without gaps, par file & model 
% --------------------------------------------------------------------
% motion data with no gaps
RawdataFiles{1,1} = [pwd,'\Experiments\ElodiePhD_ie5motions\CSV\'];
RawdataFiles{1,2} = {[FileNameEXPnoGaps,'.csv']};
% THE OTHER SETTINGS ARE THE SAME AS IN 1-A

% --------------------------------------------------------------------
% B) Reconstruct motion file without gaps. Reconstruted mot. is TRUE motion
% --------------------------------------------------------------------
if REC_TRUE == 1
    SimType = '2';
    global TOL; TOL = 1e-5;   % Defines the tolerance used for the simulation
    INFO_TYPE = 1; % INFO_TYPE = 0 only basic info; INFO_TYPE = 1 all info
    %                Log RMR  PB POS PBt0 .pro
    ResultOptions = [ 1;  1;  1;  1;   1;   0];
    % Reconstruct motion
    runModel(SimType, INFO_TYPE, ModelIKFileName, ModelIKPath, CalibFiles, RawdataFiles, ResultsPath_ORIG, 1, ResultOptions);
end


% ==============================================================================================================================
% ==============================================================================================================================
% 4) Create computer-simulated trajectories (*.CSV) by adding noise to TRUE trajectory
% ==============================================================================================================================
% ==============================================================================================================================
% D1) Nombres ficheros con trayectorias simuladas -----------------
NoisedataFiles{1,1} = [pwd,'\Experiments\ElodiePhD_ie5motions\CSV_noise\'];
%NoisedataFiles{1,2} = {[FileNameEXPnoGaps,'_N1.csv']};
NoisedataFiles{1,2} = {[FileNameEXPnoGaps,'_N1.csv']; [FileNameEXPnoGaps,'_N2.csv']; [FileNameEXPnoGaps,'_N3.csv']; ...
                       [FileNameEXPnoGaps,'_N4.csv']; [FileNameEXPnoGaps,'_N5.csv']; [FileNameEXPnoGaps,'_N6.csv']; ...
                       [FileNameEXPnoGaps,'_N7.csv']; [FileNameEXPnoGaps,'_N8.csv']; [FileNameEXPnoGaps,'_N9.csv']; ...
                       [FileNameEXPnoGaps,'_N10.csv']}; 

if GEN_NOISE == 1
    % -----------------------------------------------------------------------
    % 2.A) Extract TRUE variables from *.POS and *.RMR (joint angles)
    % -----------------------------------------------------------------------
    % get all markers in CSV file
    [fs, SmarkNamesCSV, SmarkCoordsCSV]  = readMotionData(RawdataFiles{1,1}, RawdataFiles{1,2}{1});
    % Get variable names in *.POS (TRUE data)
    [DataPos, VarNamesPOS] = compread([ResultsPath_ORIG, FileNameEXPnoGaps, '.pos']);
    % Find which markers of CSV have been included in the model (they are in *.POS)
    [SmarkCoordsCSVinModel, SmarkNamesModel, SmarkCoordsModel] = fillSmarkData(SmarkNamesCSV, SmarkCoordsCSV, VarNamesPOS, DataPos);
    nSmarksModel = length(SmarkNamesModel);
    % Get joint angles from RMR file
    RMRstruct = loadrmr([ResultsPath_ORIG, FileNameEXPnoGaps,'.rmr']);
    % Save all data in a cell.
    TrueData(1,:) = {SmarkNamesModel, SmarkCoordsCSVinModel, SmarkCoordsModel, DataPos, VarNamesPOS, RMRstruct};
    
    % -----------------------------------------------------------------------
    % 2.B) Add noise to TRUE variables for generating simulated trajectories
    % -----------------------------------------------------------------------    
    % Define noise for the Motion capture files (*.csv). Two posible noise models
    % * Sergio noise model
    %     B + A * sin(w*t + phase)
    %       B defines constant error due to model parameters estimation error
    %       A defines skin movement artifact.
    % * Cheze noise model
    %     A * sin(w*t + phase)
    %
    % Note that Amax is defined for each coordinates and Amax for marker is different
    % Amax for marker = sqrt(3) * Amax for coord.
    % Amax for coord -> Amax for marker  |  Amax for marker -> Amax for coord.
    %    5                     8.66      |        5                 2.88
    %    7.5                  12.99      |       10                 5.77
    %   10                    17.32      |       15                 8.66
    %   12.5                  21.65      |       20                11.55
    %   15                    25.98      |       25                14.43
    %   17.5                  30.31      |       30                17.32
    %   20                    34.64      |       35                20.21

    Amax = [ ...
        %Body     BEC         OSR        OSR        OSR        OSR        OSR      USR         USR         USR      FUR      FUR
        %Marker  EIAD  GR_TROCH_D  ROT_H_D_1  ROT_H_D_2  GEN_INT_D  GEN_EXT_D  ROT_J_D  CHEV_INT_D  CHEV_EXT_D  TALON_D  TARS5_D
                   12;          7;        7;        12;         7;         7;      12;         12;          7;       7;        7; ...
        %Body           OSL        OSL        OSL        OSL        OSL      USL         USL         USL      FUL      FUL      FUL
        %Marker  GR_TROCH_G  ROT_H_G_1  ROT_H_G_2  GEN_INT_G  GEN_EXT_G  ROT_J_G  CHEV_INT_G  CHEV_EXT_G  TALON_G  TARS1_G  TARS5_G
                         12;         7;         7;        12;         7;       7;         12;         12;       7;       7;       7; ...
        %Body         OAR        OAR         OAR          UAR          UAR         UAR         UAR     HAR
        %Marker ROT_B_D_1  ROT_B_D_2  COUD_EXT_D  ROT_AVB_D_1  ROT_AVB_D_2  POIG_INT_D  POIG_EXT_D  META_D
                      12;         7;          7;          12;           7;          7;         12;      12; ...
        %Body         OAL        OAL         OAL         OAL          UAL          UAL         UAL         UAL     HAL
        %Marker ROT_B_G_1  ROT_B_G_2  COUD_EXT_G  COUD_INT_G  ROT_AVB_G_1  ROT_AVB_G_2  POIG_EXT_G  POIG_INT_G  META_G
                      10;        12;         15;         15;          10;          10;         10;         15;     10; ...
        %Body     OBW   OBW  OHW     SBR     SBL      KO     KO     KO     KO
        %Marker XIPHO  MANU   C7  ACRO_D  ACRO_G  VERTEX  TEM_G  TEM_D  FRONT
                  10;    10;  10;     7;      7;     10;    10;    10;    10];

    nMarkers = length(Amax);
    Wmax     = 3*(pi/2)*ones(nMarkers,1);
    Phasemax = 2*pi*ones(nMarkers,1);

    % Add noise
    nFiles = length(NoisedataFiles{1,2});
    for i=1:nFiles
        CSVnoisePath   = NoisedataFiles{1,1};
        CSVnoiseFile_i = NoisedataFiles{1,2}{i};
        result_CSVnoise(CSVnoisePath, CSVnoiseFile_i, SmarkNamesModel(:,1), SmarkCoordsModel, Amax, Wmax, Phasemax, []);
        % se usa la funcion result_CSVnoise de la version de BioToolbox carga en el PATH
    end    
end % End of GEN_NOISE



% ==============================================================================================================================
% ==============================================================================================================================
% 5) Reconstruct and analyze computer-simulated motions
% ==============================================================================================================================
% ==============================================================================================================================
% Path for the results of computer-simulated CSV files using OTM
ResultsPath_OTM = [pwd,'\Experiments\ElodiePhD_ie5motions\RMR_OTM_CompSimMot2\'];

if REC_NOISE == 1
    % --------------------------------------------------------------------------------------------
    % 3.A) Motion reconstruction of computer-simulated motions
    % --------------------------------------------------------------------------------------------
    %SimType = '12';
    SimType = '2';
    global TOL; TOL = 1e-5;   % Defines the tolerance used for the simulation
    INFO_TYPE = 1; % INFO_TYPE = 0 only basic information; INFO_TYPE = 1 all the information displayed
    %                Log RMR  PB POS PBt0 .pro
    ResultOptions = [ 1;  1;  1;  1;   0;   0];

    % reconstruct motions
    runModel(SimType, INFO_TYPE, ModelIKFileName, ModelIKPath, CalibFiles, NoisedataFiles, ResultsPath_OTM, 1, ResultOptions);    
end



