function SmarkCoords_filt = filtTrajNaN(SmarkCoords, SmoothPar)

% FILTTRAJNAN  smooths trajectories using cubic splines. The trajectories
% can contain NaN values. If the trajectory contains NaNs, first the 
% trajectory is divided in pieces with no NaNs, then each piece is smoothed
% independently, finally the trajectory is recomposed respecting the 
% original NaNs.
%  e.g. SmarkCoords_filt = filtTrajNaN(SmarkCoords, p)
%
% Inputs:
%   SmarkCoords: array(nFrames x 3*nSmarks) with x,y,z marker coordinates
%                for each marker.
%   SmoothPar cell{2x1}:
%      SmoothPar{1} string: name of the smoothing method
%      SmoothPar{2} array (nPar x 1): contain parameter values for the method
%      SmoothPar{3} sampling freq. of motion capture system
%
% Output:
%    SmarkCoords_filt array(nFrames x 3*nSmarks) with smoothed trajectories
%

% sizes
[nFrames, nCoords] = size(SmarkCoords);
nSmarks = nCoords/3;

% Initialize
SmarkCoords_filt = zeros(nFrames, nCoords);

for i=1:nSmarks
    % get 3 coords of marker i
    TrajSmarki = SmarkCoords(:,[3*i-2:3*i]);
    
    % initialize
    TrajSmarki_filt = TrajSmarki; % if nothing "happens" original traj. is given
    
    % count total number of missing values
    nNaNs_TrajSmarki = sum( isnan(TrajSmarki(:,1)) );

    % If the trajectory has not gaps (NaNs) smooth directly
    if nNaNs_TrajSmarki == 0; % x,y & z have the same NaNs. We check only x coord.
        TrajSmarki_filt = filtTraj(TrajSmarki, SmoothPar);
        
    else        
        % If the trajectory has gaps (NaNs)
        if (nFrames - nNaNs_TrajSmarki) > 3 ; % If the trajectory has at least 3 non NaN values
            
            % Initialize trajectory
            TrajSmarki_filt = NaN * ones(nFrames,3);
            
            % divide trajectory in pieces without NaNs  
            TrajNx = getVecNoNaN(TrajSmarki(:,1));
            TrajNy = getVecNoNaN(TrajSmarki(:,2));
            TrajNz = getVecNoNaN(TrajSmarki(:,3));
            
            % sizes
            nPieces = length(TrajNx); % number of pieces without NaNs
%             Inside = 0;
            for j=1:nPieces
                % If piece j has more than 2 values then smooth
                if length(TrajNx(j).Val) >= 2; %### Si utilizo spline recordar que debe ser mayor que 2
                    Inside=1;
                    Id = TrajNx(j).Id; % indices of points for current piece
                    Piece_j = [TrajNx(j).Val', TrajNy(j).Val', TrajNz(j).Val'];
                    % smooth Piece_j of marker's trajectory i                    
                    TrajSmarki_filt(Id,:) = filtTraj(Piece_j, SmoothPar);

                end
            end
                
        end
    end
    
    % store smoothed trajectory of marker i
    SmarkCoords_filt(:,[3*i-2:3*i]) = TrajSmarki_filt;

end % of for i

