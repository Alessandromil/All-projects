classdef SOLVER < handle
    %SOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Reconstruction   % All variables to solve reconstruction                        RECONSTRUCTION
        q_t              % Value of eachgeneralized coordinate at each sample time.     double [NFrames x NVars]
        Epsilon = 1e-6;  % Epsilon is the precission of the solution.
        MarkerNames      % Names of Markers ### sólo para ke funcione
        ExperLogFileId
        Settings        % Struct that contains settings for motion reconstruction
                        % Settings.Type. Inverse Kinematics or Inverse Dynamics. Possible options: IK, ID for IK+ID
                        % Settings.Display. Define ammount of feedback to user: 0(little information) 1(complete information)
                        %   ALL RESULTS VARIABLES CAN HAVE 2 POSSIBLE INTEGER VALUES: 1(Yes)/0(No).                 
                        % Settings.Results.Ramsis. Turns on/off Ramsis generation results. 
                        % Settings.Results.PAM. Turns on/off Ramsis generation results.
                        % Settings.Results.InitPosture. Posture used for initialization in Compamm format (*_t0.pb & *_t0.sim)
                        % Settings.Results.CompPlayback. Reconstructed motion in Compamm format (*.pb & *.sim)
                        % Settings.Results.Position. Position of all model elements (markers, vectors, points) in Compamm format (*.pos)
                        % Settings.Results.Velocity. Velocity of all model elements (markers, vectors, points) in Compamm format (*.vel)
                        % Settings.Results.Acceleration. Acceleration of all model elements (markers, vectors, points) in Compamm format (*.acc)
                        % Settings.Results.MarkerError. Distance between experimental-markers and model-markers in Compamm format (*.dis)
                        % Settings.Results.RawMarkerTrajectory. Experimental marker trajectories in Compamm format (extracted from oringal data) (*.traj)
                        % Settings.Results.AcondMarkerTrajectory. Aconditioned (smoothed & gaps filled) marker trajectories in Compamm format (*.amt)
                        % Settings.Results.Sensor. Variables measured by sensors in Compamm format (*.sen)
                        % Settings.Results.JointEffort. Forces & torques at all the joints in the model. Includes reactions & motor efforts
                        % Settings.Interpolation.Method. Method for interpolating gaps in marker trajectories. Options:`'linear'
                        % Settings.Interpolation.NInterFrames. Threshold (in frames) for interpolating missing markers. Options: any positive integer
                        % Settings.Smoothing.Method. Smoothing method for marker trajectories. Options: 'butter' (butterworth filter)
                        % Settings.Smoothing.CutFreq. Cutt-off frequency in Hz for butterworth filter
        
    end
    
    methods
        function S = SOLVER(Reconstruction,ExperLogFileId,Settings)
            S.Reconstruction = Reconstruction;
            S.ExperLogFileId = ExperLogFileId;
            S.Settings = Settings;
        end
        function w = autoWeight(S,i, w0_previous, g_t_coord_j, n, MaxWeightValue)

            % AUTOWEIGHT autoWeight automatic weighting factor for a given driven coordinates
            % The value of the automatic weighting factor is based on previous and future values
            % of the driven coordinate.
            %
            % The weighting scheme implemented in this function is the following:
            %   If a marker is missing w = 0
            %   If a marker is visible in frame i but it was missing in frame i-1 then w = 1/n
            %   If a marker is visible in frame i but it was missing in frame i-2 then w = 2/n
            %   If a marker is visible in frame i but it was missing in frame i-n then w = n/n
            % When approaching a missing marker the scheme is similar:
            %   If a marker is visible but will be missing in frame    i+1  then w = 1/n
            %   If a marker is visible but will be missing in frame    i+2  then w = 2/n
            %   If a marker is visible but will be missing in frame i+(n-1) then w = (n-1)/n
            %
            %   w = autoWeight(i, w0_previous, g_t_coord_j, n, MaxWeightValue);
            %
            %   Inputs:
            %     + i          : is the number of the current frame. It goes from 1 to nFrames. Integer value
            %     + w0_previous: is the weight associated to driven coord. in the previous frame. Scalar value
            %     + g_t_coord_j: is the value of the driven coord. for each sample time (frame). Double array (nFrames x 1)
            %     + n          : is the number of steps required to reach weight value equal to 'MaxWeightValue' Examples:
            %                      n=2 MaxWeightValue*(0.5;1.0)                              MaxWeightValue*(1/2, 2/2)
            %                      n=3 MaxWeightValue*(0.333;0.666;1.0)                      MaxWeightValue*(1/3, 2/3, 3/3)
            %                      n=4 MaxWeightValue*(0.25;0.5;0.75;1.0)                    MaxWeightValue*(1/4, 2/4, 3/4, 4/4)
            %                      n=5 MaxWeightValue*(0.2;0.4;0.6;0.8;1.0)                  MaxWeightValue*(1/5, 2/5, 3/5, 4/5, 5/5)
            %                      n=6 MaxWeightValue*(0.16667;0.33333;0.5;0.66667;0.83333)  MaxWeightValue*(1/6, 2/6, 3/6, 4/6, 5/6, 6/6)
            %     + MaxWeightValue is the maximum weight associated to a driven coordinate
            %   Outputs:
            %     + w weighting factor for a coordinate in frame i taking into account pass and future values of the coordinate
            % EXAMPLE FOR n = 5 ----------------------------------------------------------------------
            %             if w0_previous == 0   % previous value
            %                 w = 0.2; % always possible
            %             elseif w0_previous == 0.2
            %                 if isnan(g_t_coord_j(i+1))
            %                     w = 0.2;
            %                 else
            %                     w = 0.4;
            %                 end
            %             elseif w0_previous == 0.4
            %                 if isnan(g_t_coord_j(i+1))
            %                     w = 0.2;
            %                 else
            %                     if isnan(g_t_coord_j(i+2))
            %                         w = 0.4;
            %                     else
            %                         w = 0.6;
            %                     end
            %                 end
            %             elseif w0_previous == 0.6
            %                 if isnan(g_t_coord_j(i+2))
            %                     w = 0.4;
            %                 else
            %                     if isnan(g_t_coord_j(i+3))
            %                         w = 0.6;
            %                     else
            %                         w = 0.8;
            %                     end
            %                 end
            %             elseif w0_previous == 0.8
            %                 if isnan(g_t_coord_j(i+3))
            %                     w = 0.6;
            %                 else
            %                     if isnan(g_t_coord_j(i+4))
            %                         w = 0.8;
            %                     else
            %                         w = 1.0;
            %                     end
            %                 end
            %             elseif w0_previous == 1.0
            %                 if isnan(g_t_coord_j(i+4))
            %                     w = 0.8;
            %                 else
            %                     w = 1.0;
            %                 end
            %             end
            % -------------------------------------------------------------------------------------
            
            % INITIALIZE -------------------------------
            w = [];
            
            nSamples = length(g_t_coord_j);
            if isnan(g_t_coord_j(i)) % coordinate is missing.
                w = 0;
            else
                % we now now that in the current frame the driven coordinate is not NaN
                if i == 1 % first frame!!!
                    % Initialize all weights taking into account value of coordinates in next frames
                    if nSamples > n % the most common case
                        for j = 1:n-1
                            if isnan(g_t_coord_j(i+j))
                                w = MaxWeightValue*(j/n);
                                break; % end FOR loop
                            end
                        end
                        if isempty(w)
                            w = MaxWeightValue;
                        end
                    else
                        % if nSamples <= n
                        % For a small number of frames it has no sense to increase and decrease
                        % the weighting factors of the missing markers.
                        % If a marker is missing weight=0 and if it is visible weight=1
                        w = MaxWeightValue;
                    end
                    
                else % all of the frame except for the firs one
                    if nSamples > n
                        % n=2 (0.5;1.0) (1/2, 2/2)
                        % n=3 (0.333;0.666;1.0) (1/3, 2/3, 3/3)
                        % n=4 (0.25;0.5;0.75;1.0) (1/4, 2/4, 3/4, 4/4)
                        % n=5 (0.2;0.4;0.6;0.8;1.0) (1/5, 2/5, 3/5, 4/5, 5/5)
                        % n=6 (0.16667;0.33333;0.5;0.66667;0.83333) (1/6, 2/6, 3/6, 4/6, 5/6, 6/6)
                        
                        % The current value of the driven coordinate is not NaN!! checked before
                        % Assign weights depending of the previous weight (i-1) and next frames.
                        % The next frames checked depend on the value of the previous weight.
                        % The higher the value of the previous weight the farther the function check
                        % the next frames because more frame are required to decrease the value of the
                        % weight to zero is a NaN is found.
                        
                        % Previous weighting factor = 0 --------------------------------
                        if w0_previous == 0   % previous value
                            w = MaxWeightValue*(1/n); % always possible
                            % Previous weighting factor = 1/n ------------------------------
                        elseif w0_previous == MaxWeightValue*(1 /n)
                            if (nSamples-i) >= 1
                                if isnan(g_t_coord_j(i+1))  % check frame (i+1)
                                    w = MaxWeightValue*(1/n);
                                else
                                    w = MaxWeightValue*(2/n);
                                end
                            else
                                w = MaxWeightValue*(2/n);
                            end
                        end
                        % Previous weighting factor = (j-1)/n --------------------------
                        for j=3:n
                            if w0_previous == MaxWeightValue*((j-1)/n)
                                if (nSamples-i) >= (j-2)
                                    if isnan(g_t_coord_j(i+(j-2)))  % check frame (i+(j-2))
                                        w = MaxWeightValue*((j-2)/n);
                                    else
                                        if (nSamples-i) >= (j-1)
                                            if isnan(g_t_coord_j(i+(j-1)))  % check frame (i+(j-1))
                                                w = MaxWeightValue*((j-1)/n);
                                            else
                                                w = MaxWeightValue*(j/n);
                                            end
                                        else
                                            w = MaxWeightValue*(j/n);
                                        end
                                    end
                                else
                                    w = MaxWeightValue*(j/n);
                                end
                            end
                        end
                        % Previous weighting factor = 1 ------------------------------
                        if w0_previous == MaxWeightValue*(n/n)
                            if (nSamples-i) >= (n-1)
                                if isnan(g_t_coord_j(i+(n-1)))  % check frame (i+(n-1))
                                    w = MaxWeightValue*((n-1)/n);
                                else
                                    w = MaxWeightValue*(n/n);
                                end
                            else
                                w = MaxWeightValue*(n/n);
                            end
                        end
                        
                    else
                        % For a small number of frames it has no sense to increase and decrease
                        % the weighting factors of the missing markers.
                        % If a marker is missing weight=0 and if it is visible weight='MaxWeightValue'
                        w = MaxWeightValue;
                    end
                end % end of else
            end % end of else
        end
        function beta = lineSearchW(S,q0, Deltaq, g0, wm0, gs0, ws0)
            % LINESEARCH returns the optimal length of the step
            % Bang-bang line search
            
            alpha = (1:-0.05:0)';
            n = length(alpha);
            nInps = length(g0);
            
            % initialize
            beta = 0;
            
            % error for different values of beta
            for i = 1:n
                q1 = q0 + alpha(i) * Deltaq;
%                 Phi_q1 = feval(S.Reconstruction.FilePhiName, q1,S.Reconstruction.Subject); %, ws0, gs0);
                Phi_q1 = feval(S.Reconstruction.FilePhiName, q1,S.Reconstruction.Par); %, ws0, gs0);
                z1 = q1(S.Reconstruction.PoszInq);
                % deal with NaN is g0
                Psi = (z1 - g0);
                for j = 1:nInps
                    if wm0(j) == 0 % eliminate driver constraint
                        Psi(j) = 0;
                    else
                        Psi(j) = wm0(j) * Psi(j);
                    end
                end
                objfun  = 0.5 * Psi' * Psi;
                normPhi = norm(Phi_q1);
                
                if i == 1
                    Min_objfun  = objfun;
                    Min_normPhi = normPhi;
                else
                    if objfun >= Min_objfun
                        beta = alpha(i-1);
                        break
                    elseif normPhi >= Min_normPhi
                        beta = alpha(i-1);
                        break
                    else
                        Min_objfun = objfun;
                        Min_normPhi = normPhi;
                    end
                end
            end
            % --------------------------------------------------------------------------
            % If there is not a minimum ...
            % --------------------------------------------------------------------------
            if beta == 0
                beta = 0.101;
            end
        end
        function beta = lineSearchW2(S,q0, Deltaq, g0, wm0, gs0, ws0)

            % LINESEARCH returns the optimal length of the step
            % Merit function line search
            
            alpha = (1:-0.05:0)';
            n = length(alpha);
            nInps = length(g0);
            
            % initialize
            beta  = 0;
            gamma = 1;
            
            % error for different values of beta
            for i = 1:n
                q1 = q0 + alpha(i) * Deltaq;
%                 Phi_q1 = feval(S.Reconstruction.FilePhiName, q1, S.Reconstruction.Subject);%%, ws0, gs0);
                Phi_q1 = feval(S.Reconstruction.FilePhiName, q1, S.Reconstruction.Par);%%, ws0, gs0);
                z1 = q1(S.Reconstruction.PoszInq);
                % deal with NaN is g0
                Psi = (z1 - g0);
                for j = 1:nInps
                    if wm0(j) == 0 % eliminate driver constraint
                        Psi(j) = 0;
                    else
                        Psi(j) = wm0(j) * Psi(j);
                    end
                end
                objfun  = 0.5 * Psi' * Psi;
                normPhi = norm(Phi_q1);
                
                if i == 1
                    if objfun == 0
                        mu = Inf;
                    else
                        mu = normPhi/objfun;
                    end
                    Min_meritFun = objfun + gamma *(1/mu) * normPhi;
                else
                    meritFun = objfun + gamma *(1/mu) * normPhi;
                    if meritFun >= Min_meritFun
                        beta = alpha(i-1);
                        break
                    else
                        Min_meritFun = meritFun;
                    end
                end
            end
            % --------------------------------------------------------------------------
            % If there is not a minimum ...
            % --------------------------------------------------------------------------
            if beta == 0
                beta = 0.101;
            end
        end
        function [wm0,wm0_all] = setWeightingFactors(S,g0,i,wm0_all)
            missing = isnan(g0); % vector with 1's in NaN positions and 0's in non-NaN positions
            NWeights_Wm = size(S.Reconstruction.Wm,1);
            wm0_tmp  = zeros(NWeights_Wm,1);
%             wm0_all  = zeros(NFrames, NWeights_Wm);
            for j = 1:NWeights_Wm
                wm_j = S.Reconstruction.Wm{j};
                if isempty(wm_j)
                    % auto weigthing factors for missing marker problem.
                    MaxWeightValue = 1.0;
                    Nf = 8;
                    if i==1 % The first frame
                        wm0_tmp(j) = S.autoWeight(i, wm0_tmp(j), S.Reconstruction.g_t(:,j), Nf, MaxWeightValue);
                    else
                        wm0_tmp(j) = S.autoWeight(i, wm0_all(i-1,j), S.Reconstruction.g_t(:,j), Nf, MaxWeightValue);
                    end
%                     if ~isempty(WEIGHT_TYPE) && ( strcmpi(WEIGHT_TYPE,'lu') |  strcmpi(WEIGHT_TYPE,'exp') )
%                         % Weighting factors proportional to measurement error.
%                         if i > 5 | (StartWeighting & wm0(j) > 0)  % wm0 == 0 when marker is missing / if you change i>1 loi
%                             wm0(j) = wm0_tmp(j) * wm0_propErr(j);
%                         else
%                             wm0(j) = wm0_tmp(j);
%                         end
%                         
%                     elseif ~isempty(WEIGHT_TYPE) && strcmpi(WEIGHT_TYPE,'const')
%                         wm0(j) = wm0_tmp(j);
%                         
%                     elseif isempty(WEIGHT_TYPE) % if weighting is not defined
                         wm0(j) = wm0_tmp(j);
%                         
%                     else
%                         error('option WEIGHT_TYPE = ',WEIGHT_TYPE,' not valid');
%                         
%                     end
%                     
%                     
                elseif isscalar(wm_j)
                    if missing(j) % value of driven coordinate is NaN in the current frame
                        wm0(j) = 0;
                    else
                        wm0(j) = wm_j;
                    end
%                     
%                 elseif isvector(wm_j) && ~iscell(wm_j) % This is a double array (nFrames x 1)
%                     if missing(j) % value of driven coordinate is NaN in the current frame
%                         wm0(j) = 0;
%                     else
%                         wm0(j) = wm_j(i); % weight for the current sample time
%                     end
%                     
%                 elseif isvector(wm_j) && iscell(wm_j) && length(wm_j) == 3  % This is a cell (3 x 1)
%                     if missing(j) % value of driven coordinate is NaN in the current frame
%                         wm0(j) = 0;
%                     else
%                         % ----------------------------
%                         % Set weighting factor value
%                         % ----------------------------
%                         % this constraint is active when marker wm_j{1} is missing
%                         SmarkFirstCoordIndex = wm_j{1};
%                         Wmissing             = wm_j{2};
%                         Wvisible             = wm_j{3};
%                         
%                         % If first skin-marker coord in z is Nan all are NaN
%                         if SmarkFirstCoordIndex == 0 | isnan(g0(SmarkFirstCoordIndex)) % first coordinate of the skin-marker in z
%                             % NOTE1: SmarkFirstCoordIndex == 0 is the marker is not in the file of model parameters
%                             %       and it appears in the weighting factor definition.
%                             % NOTE2: isnan(g0(SmarkFirstCoordIndex)) == 1 if the marker is the motion file is missing for the
%                             %        current frame
%                             
%                             % If the skin-marker associated to this driver constraint is MISSING
%                             wm0(j) = Wmissing;
%                             % ----------------------------
%                             % Set value of the coordinate
%                             % ----------------------------
%                             % The coordinate is guided to its value in the previous step.
%                             g0(j) = q0(SOL.PoszInq(j));
%                         else
%                             % If the skin-marker associated to this driver constraint is VISIBLE
%                             wm0(j) = Wvisible;
%                         end
%                         
%                     end
%                     
%                 elseif isvector(wm_j) && iscell(wm_j) && length(wm_j) == 4  % This is a cell (4 x 1)
%                     % ----------------------------
%                     % Set weighting factor value
%                     % ----------------------------
%                     % this constraint is active when marker wm_j{1} is missing
%                     VecSmarkIndex = wm_j{1};
%                     VarCoordIndex = wm_j{2}; % var z(i) will be assigned the value of VarName at previous frame
%                     Wmissing      = wm_j{3};
%                     Wvisible      = wm_j{4};
%                     
%                     % initialize
%                     AnySmarkIsMissing = 0;
%                     % Find if some Smark of the list is missing in order to active driving constraint
%                     nMissSmarks = length(VecSmarkIndex);
%                     for k=1:nMissSmarks
%                         SmarkFirstCoordIndex = VecSmarkIndex(k);
%                         % Smark k is missing? | If first skin-marker coord. in z is Nan all are NaN
%                         if SmarkFirstCoordIndex == 0 | isnan(g0(SmarkFirstCoordIndex)) % first coordinate of the skin-marker in z
%                             % NOTE1: SmarkFirstCoordIndex == 0 is the marker is not in the file of model parameters
%                             %       and it appears in the weighting factor definition.
%                             % NOTE2: isnan(g0(SmarkFirstCoordIndex)) == 1 if the marker is the motion file is missing for the
%                             %        current frame
%                             AnySmarkIsMissing = 1;
%                         end
%                     end
%                     
%                     if AnySmarkIsMissing
%                         % If the skin-marker associated to this driver constraint is MISSING
%                         wm0(j) = Wmissing;
%                         
%                         % set "guiding-value" of z(i) to the value of VarName at previous frame
%                         if VarCoordIndex < 0
%                             g0(j) = - q0( abs(VarCoordIndex) );
%                         else
%                             g0(j) = q0(VarCoordIndex);
%                         end
%                     else
%                         % If the skin-marker associated to this driver constraint is VISIBLE
%                         wm0(j) = Wvisible;
%                     end
%                     
                else
                    error(['Unknown data format for the definion of weighting factors in wm(',num2str(j),')']);
                end
            end
            % store the value of the weight for each frame
            wm0_all(i,:) = wm0';
    end
        function [q1, ErrorMin, NIterations, IterData] = stepWeightedOTM(S,q0, g0, wm0, gs0, ws0,LsearchON, ReordMeth)
            % sizes
            NInps     = length(S.Reconstruction.PoszInq);
            NCoords   = length(q0);
            
            % initialize
            NormPhi     = 10;
            ErrorMin    = 10;
            NIterations = 0;
            Deltaq      = zeros(NCoords,1);
            IterData    = [];
            
            % initialize Hw - weighted Hessian matrix
            Hw = sparse(NCoords, NCoords); % NCoords-by-NCoords all zero sparse matrix
            
            % initialize matrix of zeros
            ZeroMat = sparse(zeros(NCoords, NCoords));
            
            % Data for the first iteration
%             Phi_q0 = feval(S.Reconstruction.FilePhiName, q0, S.Reconstruction.Subject);
            Phi_q0 = feval(S.Reconstruction.FilePhiName, q0, S.Reconstruction.Par);
            
            % ----------------------------------------------------------------------------------------
            % iteration
            % ----------------------------------------------------------------------------------------
            % NLP problem is:
            %       minimize     f(q)   = (1/2) * PsiT(q) * W * Psi(q) = (1/2)* deltaqT * H * Deltaq + gT(q) * Deltaq + f(q)
            %       subject to   Phi(q) = 0
            % where
            %       Phi = vector of kinematic constraints
            %       Psi = vector of driver constraints. Psi = S*q-g
            %
            % The resultant linear system of equations is:
            %   [H  A][Deltaq] = [g] ->  C*y=d
            %   [A  0][Lambda]   [b]
            %
            % where
            %    A =  PhiqT(q0) * Phiq(q0)
            %    b = -PhiqT(q0) * Phi(q0)
            %    H = S'*S, Hessian matrix of f(q)
            %    g = S'*(Sq-d), Gradient vector of f(q)
            % with weighting matrix W then
            %    Hw = S'* W * S,  weighted Hessian matrix of f(q)
            %    gw = S'* W *(Sq-d), weighted Gradient vector of f(q)
            
            
            while NormPhi > S.Epsilon
                % Initialize cpu timer
                Tini = cputime;
                
                % fill Hw (weighted Hessian matrix) -------------------------
                for i = 1:NInps
                    Hw(S.Reconstruction.PoszInq(i),S.Reconstruction.PoszInq(i)) = wm0(i);
                end
                % fill C -------------------------
%                 Phiq_q0 = feval(S.Reconstruction.FilePhiqName, q0,S.Reconstruction.Subject);
                Phiq_q0 = feval(S.Reconstruction.FilePhiqName, q0,S.Reconstruction.Par);
                A = Phiq_q0' * Phiq_q0;
                C = [Hw, A; A, ZeroMat]; % H  is constant
                
                % fill hw (weighted gradient vector) ------------------------
                hw = zeros(NCoords, 1);
                z0 = q0(S.Reconstruction.PoszInq); % cambia en cada ITER
                f0 = g0 - z0;     % cambia en cada ITER
                for i = 1:NInps
                    if wm0(i) == 0 % eliminate driver constraint
                        f0(i) = 0;
                    else
                        f0(i) = wm0(i) * f0(i);
                    end
                end
                hw(S.Reconstruction.PoszInq)= f0;
                
                % fill b --------------------------
                b = -Phiq_q0' * Phi_q0;
                
                % fill d --------------------------
                d = [hw; b];
                
                
                % ------------------------------------------------------------------------
                % Solve the linear system of equations C * y = d
                % ------------------------------------------------------------------------
                
%                 if ~isempty(SOLVER)
%                     Method  = SOLVER(1);
                     Tol = eps;
%                     % Extension according to method
%                     if strcmp(Method,'A')
                        % Sparse QR method (single step) ------------------------------
%                         WarnON = 1;%
                        [y, rankC] = S.spSolverQR(C,d,Tol,ReordMeth);
                        % Extract Deltaq and lambda form y
                        Deltaq = y(1:NCoords); % x
                        %Lambda = y(NCoords+1:end);
                        
%                     elseif strcmp(Method,'B')
%                         % Sparse QR method in two steps (2s) --------------------------
%                         [Deltaq, rankC] = S.spSolverQR2s(Hw,A,d,Tol);
%                         
%                     elseif strcmp(Method,'C')
%                         % Null-space method based on QR factorizaton ------------------
%                         [Deltaq, rankC] = spSolverNS(Hw,A,-hw,b,Tol,ReordMeth);
%                         %ONLY ReordMeth = 0 gives satisfactory results!!!
%                         
%                     elseif strcmp(Method,'D')
%                         % Sparse LU method (two steps, delete LD ctrts in Phiq) -------
%                         Tol1 =  1e-8; % controls the detection of LD constraints in Phiq
%                         Tol2 = 1e-14; % controls LU factorisation of C
%                         [y, rankC] = spSolverLU2s(Hw,-hw,Phiq_q0,Phi_q0,Tol1,Tol2);
%                         %Extract Deltaq and lambda form y
%                         Deltaq = y(1:NCoords); % x
%                         %Lambda = y(NCoords+1:end);
%                         
%                     end
%                 else
%                     % if SOLVER is not defined the default Method is:
%                     % Sparse QR method (single step) ------------------------------
%                     WarnON = 1;
%                     Tol = eps;
%                     [y, rankC] = S.spSolverQR(C,d,Tol, ReordMeth);
%                     % Extract Deltaq and lambda form y
%                     Deltaq = y(1:NCoords); % x
%                     %Lambda = y(NCoords+1:end);
%                     
%                 end
                
                
                % ---------------------------------------------------------------
                % calc q0 & Phi(q0) for the next iteration
                % ---------------------------------------------------------------
                if LsearchON == 1
                    % Opcion 1: Line search bang-bang
                    beta = S.lineSearchW(q0, Deltaq, g0, wm0, gs0, ws0);
                    q0 = q0 + beta * Deltaq;
                elseif LsearchON == 2
                    % Opcion 2: Line search Meritfun
                    beta = S.lineSearchW2(q0, Deltaq, g0, wm0, gs0, ws0);
%                     beta = S.lineSearchW2_old(q0, Deltaq, g0, wm0, gs0, ws0);
                    q0 = q0 + beta * Deltaq;
                else
                    % Opcion 0: FULL step
                    q0 = q0 + Deltaq;
                end
%                 Phi_q0 = feval(S.Reconstruction.FilePhiName, q0, S.Reconstruction.Subject);
                Phi_q0 = feval(S.Reconstruction.FilePhiName, q0, S.Reconstruction.Par);
                
                
                % ---------------------------------------------------------------
                % Objective function
                % ---------------------------------------------------------------
                Psi = (q0(S.Reconstruction.PoszInq) - g0);
                for j = 1:NInps
                    if wm0(j) == 0 % eliminate driver constraint
                        Psi(j) = 0;
                    else
                        Psi(j) = wm0(j) * Psi(j);
                    end
                end
                objfun = 0.5 * Psi' * Psi;
                
                
                % ---------------------------------------------------------------
                % calculate error
                % ---------------------------------------------------------------
                NormPhi    = norm(Phi_q0);
                %normInfPhi = norm(Phi_q0,inf);
                normDeltaq = norm(Deltaq);
                [MaxPhiValue, MaxPhiIndex]    = max(abs(Phi_q0));
                [MaxPhiValues, MaxPhiIndexes] = sort(abs(Phi_q0));
                [MaxObjfun, MaxObjfunIndex]   = max(abs(Psi));
                [ObjfunValues, ObjfunIndexes] = sort(abs(Psi));
                MaxPhiIndexes = flipud(MaxPhiIndexes);
                MaxPhiValues  = flipud(MaxPhiValues);
                ObjfunIndexes = flipud(ObjfunIndexes);
                ObjfunValues  = flipud(ObjfunValues);
                
                
                % ---------------------------------------------------------------
                % Update NIterations
                % ---------------------------------------------------------------
                NIterations = NIterations + 1;
                
                
                % ---------------------------------------------------------------
                % Check errors
                % ---------------------------------------------------------------
                if (NormPhi < ErrorMin)  ||  (NIterations == 1)
                    ErrorMin = NormPhi;
                    q0Min = q0;
                end
                
                % ---------------------------------------------------------------
                % Measure iteration time
                % ---------------------------------------------------------------
                Tend = cputime-Tini;
                
                
                % ---------------------------------------------------------------
                % Display information
                % ---------------------------------------------------------------
                if LsearchON
                    % store data for efficiency analysis of the method
                    IterData = [IterData; [NormPhi, objfun, normDeltaq, beta, rankC, Tend]];
                    
                    str = ['      Iter Nº = ',num2str(NIterations,'%2d'),',  NormPhi = ',num2str(NormPhi,'%6.3e'), ...
                        ',  objfun = ',num2str(objfun,'%6.3e'),',  normDeltaq = ',num2str(normDeltaq,'%6.3e'), ...
                        ',  beta = ',num2str(beta),',  rank(C) = ',num2str(rankC,'%2d')];
                else
                    % store data for efficiency analysis of the method
                    IterData = [IterData; [NormPhi, objfun, normDeltaq, rankC, Tend]];
                    
                    str = ['      Iter Nº = ',num2str(NIterations,'%2d'),',  NormPhi = ',num2str(NormPhi,'%6.6e'), ...
                        ',  objfun = ',num2str(objfun,'%6.6e'),',  normDeltaq = ',num2str(normDeltaq,'%6.3e'), ...
                        ',  rank(C) = ',num2str(rankC,'%2d')];
                end
                if(S.Settings.Display == 2)
                    disp(str);
                end
                fprintf(S.ExperLogFileId, '%s\n', str);

%                       for k=1:6
%                           disp(['       max error in phi(',num2str(MaxPhiIndexes(k)),') = ' num2str(MaxPhiValues(k))]);
%                       end
                
                if (NIterations >= 25) && (NormPhi > S.Epsilon)
                    % Maximum number of steps reached and NormPhi is still > Epsilon
                    % Select solution with minimum error
                    % The second condition in the 'if' is to avoid entering the 'if'
                    % when the NormPhi is less than Epsilon in the last iteration
                    str = sprintf(['      Number of iterations > 25\n', ...
                        '    Solution with the minimum error selected. Error = ',num2str(ErrorMin)]);
                    if(S.Settings.Display == 1 || S.Settings.Display == 2)
                        disp(str);
                    end
                    fprintf(S.ExperLogFileId, '%s\n', str);
                    q0 = q0Min;
                    break
                end
                
                
                %     % --------------------------------------------------------
                %     % Check if the simulation has been stoped by the user
                %     % --------------------------------------------------------
                %     if ~isempty(STOPSIM)
                %         % The simManager window is not refreshed HERE because then
                %         % it is not possible to stop the simulation.
                %         % I ignore the reason why it works but it works
                %         if STOPSIM == 1
                %             q1 = q0';
                %             % Return control to previous function
                %             return
                %         end
                %     end
                
            end  % end of iteration
            
            % solution found
            q1 = q0';

    
        end
        function [y, rankC] = spSolverQR(S,varargin)

            % SPSOLVERQR solves a sparse linear system C * y = d using
            % a sparse solver based on the QR decomposition. See function
            % QR for more details.
            %
            %   Inputs:
            %     + C is the sparse system matrix
            %     + d is the vector of independent terms
            %     + Tol is the tolerance
            %     + WarnON is an integer (1 or 0) to connect or disconnect warnings
            %     + ReordMeth defines the reordering method used before the factorization
            %       of the linear system of equations. There are four possibilities:
            %          0 - No reordering
            %          1 - colperm (columns are ordered according to increasing count of nonzero entries)
            %          2 - symrcm (reverse Cuthill-McKee ordering)
            %          3 - symamd (approximate minimum degree ordering)
            %   Outputs:
            %     + y is the solution vector
            %     + rankC is the rank of the system matrix C
            
            % set input variables
            if nargin == 4
                C         = varargin{1};
                d         = varargin{2};
                Tol       = varargin{3};
                ReordMeth = 0; % Default - No reordering
            elseif nargin == 5
                C         = varargin{1};
                d         = varargin{2};
                Tol       = varargin{3};
                ReordMeth = varargin{4};
            end
            
            % size
            [nRows, nCols] = size(C);
            
            % check sizes
            if nRows ~= nCols
                error('Algorithm only valid for square matrices');
            end
            
            % 0) Reordering method and QR decomposition of C ---------------------------------------------------------
            %    QR function returns only d0 = Q'*d and R upper triangular matrix
            if ReordMeth == 0              % No reordering
                [d0, R] = qr(C,d);
            elseif  ReordMeth == 1
                p = colperm(C);  % colperm (columns are ordered according to increasing count of nonzero entries)
                [d0, R] = qr(C(p,p),d(p));
            elseif  ReordMeth == 2
                p = symrcm(C);   % symrcm  (reverse Cuthill-McKee ordering)
                [d0, R] = qr(C(p,p),d(p));
            elseif  ReordMeth == 3
                p = symamd(C);   % symamd  (approximate minimum degree ordering)
                [d0, R] = qr(C(p,p),d(p));
            end
            
            % Initialize variables
            pivotColIndex    = [];
            noPivotColIndex  = [];
            nPivots          = 0;
            
            % 1) Find "Pivot column index" vector and "noPivot column index" vector ----------------------------------
            for j = 1 : nCols
                i = nPivots + 1;
                if(abs(R(i,j)) < Tol)
                    noPivotColIndex = [noPivotColIndex; j];
                else
                    nPivots = nPivots + 1;
                    pivotColIndex = [pivotColIndex; j];
                end
            end
            
            % 2) Build permutation matrix (column reordering) --------------------------------------------------------
            %    First: Pivots Cols. After noPivot Cols are set to zero
            E = speye(nCols);
            E1 = [E(:,pivotColIndex) 0*E(:,noPivotColIndex)];
            
            % 3) Reorder R to get an upper diagonal matrix, set noPivots Cols to zero --------------------------------
            R1 = R * E1;
            
            % 4) Set variables associated to noPivots Cols to 1. This is done by: ------------------------------------
            %       1- Setting pivot to 1
            %       2- Setting independent term in d1 to 0
            d1 = d0;
            for i = (nPivots + 1) : nRows
                R1(i,i) = 1; % Pivot of "noPivot Cols" are set to 1
                d1(i) = 0;   % associated independent term set to 0
            end

            % 5) System resolution -----------------------------------------------------------------------------------
            y1 = R1 \ d1;
            y = E1 * y1;
            % Note that E1 sets some variables in y1 to zero. Those who have been set to zero.
            
            % if ReordMeth == 0 no reordering was performed and
            % there is no need to undo the reordering
            if  ReordMeth == 1 || ReordMeth == 2 || ReordMeth == 3
                % inverse of Reordering method aplied before
                r(p) = 1:length(p);
                y = y(r);
            end
            
            % 6) Find rank of the coeffient matrix C ------------------------------------------------------------------
            rankC = nPivots;
            
            % 7) Check if the problem is under-guided -----------------------------------------------------------------
            if ReordMeth == 0%#### 0
                nVars = nRows/2;
                nVarsUnder = 0;
                ValueVariable = [];
                for i=1:length(noPivotColIndex)
                    if noPivotColIndex(i) <= nVars
                        nVarsUnder = nVarsUnder + 1;
                        ValueVariable = [ValueVariable;noPivotColIndex(i)];
                    end
                end
                if nVarsUnder > 0
                    str1 = sprintf(['      WARNING: There could be ',num2str(nVarsUnder),' variables under-guided !!!\n']);
                             str2 = '               The underguided variables are probably: ';
                    for i=1:nVarsUnder
                        str2 = [str2,S.Reconstruction.Subject.q{ValueVariable(i)},' '];
                    end
                    str3 = [str1,str2];
                    if(S.Settings.Display == 1 || S.Settings.Display == 2)
                        disp(str3);
                    end
                    fprintf(S.ExperLogFileId, '%s\n', str3);
                end
            end
        end
        function [Deltaq, rankC] = spSolverQR2s(S,H,A,d,Tol)

            % SPSOLVERQR2S solves a sparse linear system C * y = d using
            % a sparse solver based on the QR decomposition in TWO STEPS (2s)
            % See function QR for more details on QR decomposition.
            %
            %   Deltaq = spSolverQR2s(H,A,d,Tol)
            %
            %   Inputs:
            %     + H is the Hessian matrix of the objective function
            %     + A is Phiq'*Phiq
            %     + d is the vector of independent terms
            %     + Tol is the tolerance
            %   Outputs:
            %     + Deltaq is the solution vector
            %     + rankC is the rank of the system matrix C
            
            
            % The linear system of equations is:
            %   [H  A][ Deltaq] = [-g] ->  C*y=d
            %   [A  0][-Lambda]   [ b]
            % Consider the following matrices
            %   E=[H]  F=[A] d = [-g]
            %     [A]    [0]     [ b]
            % and that Lambda = -Lambda. The original variable can be obtained
            % by changing the sign of the new Lambda
            %
            % Then
            %    E*Deltaq + F*Lambda = d [1]
            %


            % 1) Form matrices E, F, d -> E*deltaq + F*Lambda = d [1]   ----------------------------------
            % sizes
            nCoords = size(A,1);
            % Definitions
            E = [H; A];
            F = [A; sparse(zeros(nCoords,nCoords))];
            % d is an input of the function


            % 2) QR decomposition of E:  E = QR ----------------------------------------------------------
            % Note that orthogonal transformations are applied to [F d] producing
            %  Q'*[F d] = [Q'*F  Q'*d] without computing Q.
            %
            % From [1] Q*R*Deltaq + F*Lambda = d
            %            R*Deltaq + Q'*F*Lambda = Q'd
            %            R*Deltaq + F0*Lambda = d0
            [Fd,R]= qr(E,[F, d]);
            d0 = Fd(:,end);
            F0 = Fd(:,1:end-1);


            % 3) Analizing the structure of [1] we can get two linear system of equations  ---------------
            %            [ R1 ] * deltaq + [ F1 ] * lambda = [ d1 ]
            %            [  0 ]            [ F2 ]            [ d2 ]
            %
            %   R1*Deltaq + F1*lambda = d1 [2]
            %               F2*lambda = d2 [3]
            
            % Definition of matrices in linear system [2]
            R1 = R(1:nCoords,:);
            F1 = F0(1:nCoords,:);
            d1 = d0(1:nCoords);

            % Definition of matrices in linear system [3]
            F2 = F0(nCoords+1:end,:);
            d2 = d0(nCoords+1:end);
            
            % Check if the problem is well formulated, i.e. all DoFs are effectively driven
            diagR1    = full(diag(R1));
            NoPivot   = find(abs(diagR1)<Tol);
            nNoPivots = length(NoPivot);
%             if  nNoPivots > 0
%                 error(['The model is not properly guided. There are ',num2str(nNoPivots),' DoFs not guided']);
%             end


            % 4) QR decomposition of matrix F2 for linear system [3]  ------------------------------------
            %    F2 = Q2*R2
            %
            %    Q2*R2*Lambda = d2        [3]
            %    R2*Lambda = Q2'*d2 = d3  [3]
            [d3, R2] = qr(F2, d2);


            % 5) Analizing the structure of [3] ----------------------------------------------------------
            %                                   [ R4| D4 ] * Lambda = [d4]
            %     R2*Lambda = Q2'*d2 = d3 ->    [  0 |  0 ]            [ 0]
            %
            % If R2 is ordered as an upper triangular matrix then we can distinguish R3 and D3
            % Note that after reordering the vector of Lambdas will change.
            %
            %    R2*Lambda = Q2'*d2 = d3 -> R2 * E2 * (E2' * Lambda) = d3
            %
            % with new variables:
            %    Lambda2 = E2' * Lambda
            %    R3      = R2  * E2
            %
            % Then
            %    R3*Lambda2 = d3
            
            % 5.a) Calculating permutation (column reordering) matrix E2 ***************
            % size
            [nRows_R2, nCols_R2] = size(R2);

            % Initialize variables
            pivotColIndex    = [];
            noPivotColIndex  = [];
            nPivots          =  0;
            
            % Find Pivot column index and noPivot column index
            for j = 1 : nCols_R2
                i = nPivots + 1;
                if(abs(R2(i,j)) < Tol)
                    noPivotColIndex = [noPivotColIndex; j];
                else
                    nPivots = nPivots + 1;
                    pivotColIndex = [pivotColIndex; j];
                end
            end

            % permutation matrix. First: Pivots Cols. After noPivot Cols
            E = speye(nCols_R2);
            E2 = E(:,[pivotColIndex' noPivotColIndex']);


            % 5.b) Calculate R3 = R2 * E2 ***********************************************
            % After reordering we have R2*E2 = [ R3 | D3 ]
            %                                  [  0 |  0 ]
            R3 = R2 * E2;


            % 5.c) Extract vars R3, D3, b4 from R3 and b3 *******************************
            R4 = R3(1:nPivots,1:nPivots);
            %D4 = R3(nPivots+1:end,nPivots+1:end); % never used
            d4 = d3(1:nPivots);
            
            % 5.d) Calculate lambda *****************************************************
            Lambda2_free = zeros(nRows_R2-nPivots,1);
            Lambda2_ind = R4\d4;
            Lambda2 = [Lambda2_ind; Lambda2_free];
            % undo the column(variable) reordering
            Lambda = E2 * Lambda2;
            
            
            % 6) From linear equation [2] Deltaq = R1\(d1-F1*Lambda) --------------------------------------------
            Deltaq = full(R1\(d1-F1*Lambda));
            
            
            % 7) Rank of the coefficient matrix  ----------------------------------------------------------------
            % We suppose that E has full column rank otherwise and error is displayed.
            % We know that E and F are L.I. then:
            %     rank(C) = rank([E F]) = rank(E) + rank(F) = nCoords + rank(F) = nCoords + rank(R3)
            rankC = nCoords + nPivots;
            
            % 7) Check if the problem is under-guided -----------------------------------------------------------------
            % if ReordMeth == 2 %#### 0
            nVars = nRows_R2/2;
            nVarsUnder = 0;
            ValueVariable = [];
            for i=1:length(noPivotColIndex)
                if noPivotColIndex(i) <= nVars
                    nVarsUnder = nVarsUnder + 1;
                    ValueVariable = [ValueVariable;noPivotColIndex(i)];
                end
            end
            if nVarsUnder > 0
                str1 = sprintf(['      WARNING: There could be ',num2str(nVarsUnder),' variables under-guided !!!\n']);
                         str2 = '               The underguided variables are probably: ';
                for i=1:nVarsUnder
                    str2 = [str2,S.Reconstruction.Subject.q{ValueVariable(i)},' '];
                end
                str3 = [str1,str2];
                if(S.Settings.Display == 1 || S.Settings.Display == 2)
                    disp(str3);
                end
                fprintf(S.ExperLogFileId, '%s\n', str3);
            end
        end
        function weightedOTM(S)
            % sizes
            [NFrames, NInputs] = size(S.Reconstruction.g_t);
            NCoords       = size(S.Reconstruction.Subject.q,1);
            NDrivenCoords = size(S.Reconstruction.z,1);
            NWeights_Wm   = size(S.Reconstruction.Wm,1);
            NWeights_Ws   = size(S.Reconstruction.Ws,1);
            NInputs_Phi   = size(S.Reconstruction.gs_t,2);
            NCtrs         = size(S.Reconstruction.Subject.Phi,1);
            velApprox = 1;
            %check sizes
            if NCoords ~= length(S.Reconstruction.q0)
                error(sprintf(['\n------------------------------------------------\n'...
                    ' q and q0 have different length!!'...
                    '\n------------------------------------------------\n']));
            end
            if NDrivenCoords ~= NInputs
                error(sprintf(['\n--------------------------------------------------------------\n'...
                    'z(driven coords) and g_t(Inputs) have different length!!'...
                    '\n--------------------------------------------------------------\n']));
            end
            if NDrivenCoords ~= NWeights_Wm
                error(sprintf(['\n------------------------------------------------\n'...
                    ' z and wm have different length!!'...
                    '\n------------------------------------------------\n']));
            end
            % ----------------------------------------------------------------------------------------
            % Newton - Rapson iteration
            % ----------------------------------------------------------------------------------------
            % display info
            str = sprintf(['   Solving problem with Optimal Tracking Method. Problem sizes:\n'...
                '    ',num2str(NCoords),' generalized coordinates\n',...
                '    ',num2str(NInputs),' driven generalized coordinates\n',...
                '    ',num2str(NCtrs),' kinematic constraints in Phi(q)\n',...
                '    Phiq has size (',num2str(NCtrs),' x ',num2str(NCoords),')\n',...
                '    Phiq''*Phiq has size (',num2str(NCoords),' x ',num2str(NCoords),')\n',...
                '    The complete matrix of the problem has size (',num2str(2*NCoords),' x ',num2str(2*NCoords),')\n']);
            if(S.Settings.Display == 1 || S.Settings.Display == 2)
                disp(str);
            end
            fprintf(S.ExperLogFileId, '%s\n', str);
            % initialize
            S.q_t = zeros(NFrames, NCoords);
                        
            % variable initialization
            StepData = cell(NFrames,4);
            Results  = zeros(NFrames,2);
            q0       = S.Reconstruction.q0;
            ws0      = zeros(NWeights_Ws,1);
            wm0_all  = zeros(NFrames, NWeights_Wm);
            ws0_all  = zeros(NFrames, NWeights_Ws);
            nIterations = 25;
            errorMin    = 10;
            DivergCount = 0;
            StartWeighting = 0;
            
            % Parameter definition
            Xframes     = 25;  % Check after 'Xframes' if convergence is ok otherwise EXIT
            SIM_END_OK  = 1;  % if simulation is finished before time this var is changed to 0
            DivergThreshold = 1e+3;
            MaxDivergFrames = 3;
            
%             % ---------------------------------
%             % Pre-process weighting definition
%             % ---------------------------------
%             wm2 = preProcWeightDef(wm, q, z);
            % --------------------------------------------------------
            % Iteration loop for each frame
            % --------------------------------------------------------
            for i=1:NFrames
                
                % Initialise CPU timer
                Tini = cputime;
                
                % info in command window
                strStep = ['    Step Nº = ',num2str(i),'/',num2str(NFrames)];
                fprintf(S.ExperLogFileId, '%s', strStep);
                % --------------------------------------------------------
                % load inputs for the current step
                % --------------------------------------------------------
                g0  = S.Reconstruction.g_t(i,:)';
                if isempty(S.Reconstruction.gs_t)
                    gs0 = [];
                else
                    gs0 = S.Reconstruction.gs_t(i,:)';
                end
                % --------------------------------------------------------
                % distance error for each marker in previous frame
                % --------------------------------------------------------
%                 if errorMin <= SOL.Epsilon
%                     StartWeighting = 1;
%                 end
%                 
%                 if i > 5 || StartWeighting
%                 end
                % ----------------------------------------------------------------------
                % Setting weghting factors for driver constraints in objective function
                % ----------------------------------------------------------------------
                [wm0,wm0_all] = S.setWeightingFactors(g0,i,wm0_all);
                % ----------------------------------------------------------------------
                % Setting weghting factors for driver constraints in equality constraints
                % ----------------------------------------------------------------------
                for j = 1:NWeights_Ws
                    ws_j = S.Reconstruction.ws{j};
                    if isscalar(ws_j)
                        ws0(j) = ws_j;
                    elseif isvector(ws_j) && ~iscell(ws_j) % This is a double array (nFrames x 1)
                        ws0(j) = ws_j(i); % weight for the current sample time
                    end
                end
                
                % store the value of the weight for each frame
                ws0_all(i,:) = ws0';
                % --------------------------------------------------------
                % Weighted Newton-Rapson iteration
                % --------------------------------------------------------
                % 1) Select iteration MODE
                %    LsearchON = 2; % 0-OFF, 1-Bang/Bang 2-meritFunction
                %    ReordMeth = 0; % 0-as is, 1-colperm, 2-symrcm, 3-symamd
                LSearchMASTER = 2;
                if (nIterations == 25)  &&  (errorMin > S.Epsilon)
                    % SAFE MODE
                    % ----------------------
                    % 1) For the first frame (i=1) nIteration is initialised to 25 and
                    %    errorMin to 10. Then we enter always this option for i=1
                    % 2) The second condition is to avoid entering here when convergence
                    %    is reached in the last iteration, i.e. nIterations = 25
                    
                    LsearchON = LSearchMASTER; % 0-OFF, 1-Bang/Bang 2-meritFunction
                    ReordMeth = 0;             % 0-as is, 1-colperm, 2-symrcm, 3-symamd
                    
                    % info in command window
                    strMode = ' - SAFE MODE';
                    if(S.Settings.Display == 2)
                        str = [strStep,strMode];
                        disp(str);
                    end
                    fprintf(S.ExperLogFileId, '%s\n', strMode);
                    
                else
                 % SPEED MODE
%                  if ~isempty(SOLVER)
%                      LSearch = SOLVER(2);
%                      if strcmp(LSearch,'1')
%                          LsearchON = LSearchMASTER;
%                          ReordMeth = 0;
%                          
%                      elseif strcmp(LSearch,'2')
%                          LsearchON = 0;
%                          ReordMeth = 0;
%                          
%                      elseif strcmp(LSearch,'3')
%                          LsearchON = LSearchMASTER;
%                          ReordMeth = 2;
%                          
%                      elseif strcmp(LSearch,'4')
                         LsearchON = 2;
                         ReordMeth = 2;
                         
%                      end
%                  else
%                      % if SOLVER is not defined the default LSearch strategy is
%                      LsearchON = 0; % before: LSearchMASTER
%                      ReordMeth = 2; % before: 0
                     
                 
                 % info in command window
                 strMode = ' - SPEED MODE';
                 if(S.Settings.Display == 2)
                     str = [strStep,strMode];
                     disp(str);
                 end
                 fprintf(S.ExperLogFileId, '%s\n', strMode);
                end
                %2) Perform iteration
%                 [q1, errorMin, nIterations, IterData] = stepweightedOTM(q0, g0, wm0, gs0, ws0, PoszInq, FilenamePhi, FilenamePhiq, Epsilon, LsearchON, ReordMeth);
                [q1, errorMin, nIterations, IterData] = S.stepWeightedOTM(q0, g0, wm0, gs0, ws0,LsearchON, ReordMeth);
                 % ---------------------------------------------------------------
                % Display information
                % ---------------------------------------------------------------
                if LsearchON                  
                    str = ['      Iter Nº = ',num2str(nIterations,'%2d'),',  NormPhi = ',num2str(IterData(nIterations,1),'%6.3e'), ...
                        ',  objfun = ',num2str(IterData(nIterations,2),'%6.3e'),',  normDeltaq = ',num2str(IterData(nIterations,3),'%6.3e'), ...
                        ',  beta = ',num2str(IterData(nIterations,4)),',  rank(C) = ',num2str(IterData(nIterations,5),'%2d')];
                else
                    % store data for efficiency analysis of the method
                    IterData = [IterData; [NormPhi, objfun, normDeltaq, rankC, Tend]];
                    
                    str = ['      Iter Nº = ',num2str(nIterations,'%2d'),',  NormPhi = ',num2str(IterData(nIterations,1),'%6.6e'), ...
                        ',  objfun = ',num2str(IterData(nIterations,2),'%6.6e'),',  normDeltaq = ',num2str(IterData(nIterations,3),'%6.3e'), ...
                        ',  rank(C) = ',num2str(IterData(nIterations,4),'%2d')];
                end
%                 disp(str);
                % Measure step time
                Tend = cputime - Tini;
                % store data of the current step for later study
                S.q_t(i,:) = q1;
                Results(i,1) = nIterations;
                Results(i,2) = errorMin;
                StepData(i,:)  = {nIterations, errorMin, IterData, Tend};
                
                % Initial approximation for next step
                if velApprox == 1
                    q0 = velapprox(S.q_t, i);
                else
                    q0 = S.q_t(i,:)';
                end
                % --------------------------------------------------------
                % Check for Divergence. If divergence ocurrs EXIT
                % --------------------------------------------------------
                % After 3 consecutive frames of divergence EXIT reconstruction
                if errorMin > DivergThreshold
                    DivergCount = DivergCount + 1;
                    if DivergCount == MaxDivergFrames
                        % Resize the variables previously allocated
                        S.q_t      = S.q_t(1:i,:);
                        StepData = StepData(1:i,:);
                        Results  = Results(1:i,:);
                        % Set variable SIM_END_OK to 0 to warn other functions or procedures
                        SIM_END_OK  = 0;
                        % send a warning message
                        str = sprintf(['    WARNING: Divergence detected during ',num2str(MaxDivergFrames),' consecutive frames. Simulation stopped']);
                        disp(str);
                        fprintf(S.ExperLogFileId, '%s\n', str);
                        % exit from the FOR-loop
                        break
                    end
                else
                    DivergCount = 0;
                end
                
            end % end of for nFrames
            % ----------------------------------------------------------------------------------------
            % Results convergence clasification
            % ----------------------------------------------------------------------------------------
            % initialize
            frame=zeros(9,1);
            iterplus=0;
            
            % sizes
            nFrames=size(Results,1);
            
            for i=1:nFrames
                if     (   0 < Results(i,2)) && (Results(i,2) <= 1e-7)
                    frame(1,1) = frame(1,1)+1;
                elseif (1e-7 < Results(i,2)) && (Results(i,2) <= 1e-6)
                    frame(2,1) = frame(2)+1;
                elseif (1e-6 < Results(i,2)) && (Results(i,2) <= 1e-5)
                    frame(3,1) = frame(3,1)+1;
                elseif (1e-5 < Results(i,2)) && (Results(i,2) <= 1e-4)
                    frame(4,1) = frame(4,1)+1;
                elseif (1e-4 < Results(i,2)) && (Results(i,2) <= 1e-3)
                    frame(5,1) = frame(5,1)+1;
                elseif (1e-3 < Results(i,2)) && (Results(i,2) <= 1e-2)
                    frame(6,1) = frame(6,1)+1;
                elseif (1e-2 < Results(i,2)) && (Results(i,2) <= 1e-1)
                    frame(7,1) = frame(7,1)+1;
                elseif (1e-1 < Results(i,2)) && (Results(i,2) <=  1.0)
                    frame(8,1) = frame(8,1)+1;
                elseif Results(i,2) > 1
                    frame(9,1) = frame(9,1)+1;
                end
                iterplus = iterplus+Results(i,1);
            end
            
            %Number of blank spaces in frames
            for i =1:length(frame)
                if     (frame(i,1) >=    0) && (frame(i,1) <    10)
                    frame(i,2) = 3;
                elseif (frame(i,1) >=   10) && (frame(i,1) <   100)
                    frame(i,2) = 2;
                elseif (frame(i,1) >=  100) && (frame(i,1) <  1000)
                    frame(i,2) = 1;
                elseif (frame(i,1) >= 1000) && (frame(i,1) < 10000)
                    frame (i,2) = 0;
                end
            end
            
            percentageframe(:,1) = frame(:,1)/nFrames*100;
            
            %Number of blank spaces in percentage
            for i =1:length(percentageframe)
                if     (percentageframe(i,1) >=  0) && (percentageframe(i,1) <  10)
                    percentageframe(i,2) = 2;
                elseif (percentageframe(i,1) >= 10) && (percentageframe(i,1) < 100)
                    percentageframe(i,2) = 1;
                elseif percentageframe(i,1) == 100
                    percentageframe (i,2) = 0;
                end
            end
            
            % ---------------------------------------------------------------
            % Display convergence clasification
            % ---------------------------------------------------------------
            str = sprintf(['\n    --------------------------------------------------\n', ...
                '     Solver convergence report. TOL = ',num2str(S.Epsilon),'\n', ... % Sergio tenía TOL no se para que
                '    --------------------------------------------------\n',...
                '      Error Range      Frames       Percentages\n',...
                '     -------------   ----------   ---------------\n'...
                '       < 1e-7           ',blanks(frame(1,2)),num2str(frame (1,1),'%4d'),'          ', ...
                blanks(percentageframe(1,2)),num2str(percentageframe (1,1),'%4.2f'),'\n', ...
                '     [1e-7, 1e-6]       ',blanks(frame(2,2)),num2str(frame (2,1),'%4d'),'          ', ...
                blanks(percentageframe(2,2)),num2str(percentageframe (2,1),'%4.2f'),'\n', ...
                '     [1e-6, 1e-5]       ',blanks(frame(3,2)),num2str(frame (3,1),'%4d'),'          ', ...
                blanks(percentageframe(3,2)),num2str(percentageframe (3,1),'%4.2f'),'\n', ...
                '     [1e-5, 1e-4]       ',blanks(frame(4,2)),num2str(frame (4,1),'%4d'),'          ', ...
                blanks(percentageframe(4,2)),num2str(percentageframe (4,1),'%4.2f'),'\n', ...
                '     [1e-4, 1e-3]       ',blanks(frame(5,2)),num2str(frame (5,1),'%4d'),'          ', ...
                blanks(percentageframe(5,2)),num2str(percentageframe (5,1),'%4.2f'),'\n', ...
                '     [1e-3, 1e-2]       ',blanks(frame(6,2)),num2str(frame (6,1),'%4d'),'          ', ...
                blanks(percentageframe(6,2)),num2str(percentageframe (6,1),'%4.2f'),'\n', ...
                '     [1e-2, 1e-1]       ',blanks(frame(7,2)),num2str(frame (7,1),'%4d'),'          ', ...
                blanks(percentageframe(7,2)),num2str(percentageframe (7,1),'%4.2f'),'\n', ...
                '     [1e-1, 1.0 ]       ',blanks(frame(8,2)),num2str(frame (8,1),'%4d'),'          ', ...
                blanks(percentageframe(8,2)),num2str(percentageframe (8,1),'%4.2f'),'\n', ...
                '       > 1.0            ',blanks(frame(9,2)),num2str(frame (9,1),'%4d'),'          ', ...
                blanks(percentageframe(9,2)),num2str(percentageframe (9,1),'%4.2f'),'\n', ...
                '     Number of total iterations: ', num2str(iterplus)]);
            if(S.Settings.Display == 1 || S.Settings.Display == 2)
                disp(str);
            end
            fprintf(S.ExperLogFileId, '%s\n', str);

            % ---------------------------------------------------------------
            % Check for number of frames with convergence above TOL value
            % ---------------------------------------------------------------
            ConvFramePerTol = frame(:,1); % array(9x1)
            if     (   0 < S.Epsilon) && (S.Epsilon <= 1e-7)
                nFramesAboveTol = sum(ConvFramePerTol(2:end));
            elseif (1e-7 < S.Epsilon) && (S.Epsilon <= 1e-6)
                nFramesAboveTol = sum(ConvFramePerTol(3:end));
            elseif (1e-6 < S.Epsilon) && (S.Epsilon <= 1e-5)
                nFramesAboveTol = sum(ConvFramePerTol(4:end));
            elseif (1e-5 < S.Epsilon) && (S.Epsilon <= 1e-4)
                nFramesAboveTol = sum(ConvFramePerTol(5:end));
            elseif (1e-4 < S.Epsilon) && (S.Epsilon <= 1e-3)
                nFramesAboveTol = sum(ConvFramePerTol(6:end));
            elseif (1e-3 < S.Epsilon) && (S.Epsilon <= 1e-2)
                nFramesAboveTol = sum(ConvFramePerTol(7:end));
            elseif (1e-2 < S.Epsilon) && (S.Epsilon <= 1e-1)
                nFramesAboveTol = sum(ConvFramePerTol(8:end));
            elseif (1e-1 < S.Epsilon)
                nFramesAboveTol = sum(ConvFramePerTol(9:end));
            end
            
            PercentFramesAboveTol = 100*nFramesAboveTol/nFrames;
            
            % If simulation finished ok (SIM_END_OK=1) and % of frames that not converged > 5% display warning.
            if PercentFramesAboveTol > 5 && SIM_END_OK
                str = sprintf(['    WARNING: ',num2str(PercentFramesAboveTol,'%3.1f'),'%% of the frames have not converged to TOL=',num2str(S.Epsilon)]);
                disp(str);
                fprintf(S.ExperLogFileId, '%s\n', str);
            end
            S.Reconstruction.q_t = S.q_t;
            
%             S.calcVelAccel();
%             S.Reconstruction.qdot_t = S.qdot_t;
%             S.Reconstruction.qdot2_t = S.qdot2_t;
        end
        
 
    end
    
end

