function [R, d] = calcOptPose(Glob_Pos_Markers, TCS_Pos_Markers)
            % CALCOPTPOSE Algorithm based on SODERKVIST and WEDIN (1993) that calculates
            % the optimal pose of a body given the global and local coordinates of n points (markers) of the body.
            %
            %   [R, d] = calcOptPose(Glob_Pos_Markers, TCS_Pos_Markers)
            %   Inputs:
            %     + Glob_Pos_Markers are the global coordinates of the markers located in
            %       the body - double array (3 x NMarkers)
            %     + TCS_Pos_Markers are the local coordinates of the markers located in
            %       the body referred to TCS - double array(3 x NMarkers)
            %   Outputs:
            %     + R is the rotation matrix between the global coordinate system and
            %       TCS - double array(3x3)
            %     + d is the position vector of the TCS origin referred to global coordinate
            %       system - double array(3x1)

            globMark = Glob_Pos_Markers;
            locMark  = TCS_Pos_Markers;
            
            NMarkers = size(globMark,2);
            
            % 1) -----------------------
            if isempty(locMark)
                meanGlobMark = zeros(3, 1);
                meanLocMark = zeros(3, 1);
            else
                meanGlobMark = (1/NMarkers) * sum(globMark, 2); % sum by rows
                meanLocMark  = (1/NMarkers) * sum(locMark, 2); % sum by rows
            end
            
            % 2) -----------------------
            A = zeros(3, NMarkers);
            for i=1 : NMarkers
                A(:, i) =  locMark(:, i) - meanLocMark;
            end
            B = zeros(3, NMarkers);
            for i=1 : NMarkers
                B(:, i) =  globMark(:, i) - meanGlobMark;
            end
            
            % 3) -----------------------
            C = B * A';
            
            % 4) -----------------------
            [P, Gamma, Q] = svd(C);
            
            % 5) -----------------------
            R = P * diag([1; 1; det(P*Q')]) * Q';
            
            % 6) -----------------------
            d = meanGlobMark - R * meanLocMark;
end