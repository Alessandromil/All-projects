function [Glob_Pos_OrCoordSist,Glob_R_BodyCoordSist, BodyCoordSist_Pos_Marker]=calcMarkersTCS(Glob_Pos_Markers)
            NMarkers = size(Glob_Pos_Markers,2); %(S.LocalMarkers,1);
            % Centroid
            Glob_Pos_OrCoordSist = (1 / NMarkers) * sum(Glob_Pos_Markers, 2); % sum by rows
            % Markers local coordinates in TCS
            for i=1:NMarkers
                TCS_Pos_Markers(:, i) = Glob_Pos_Markers(:, i) - Glob_Pos_OrCoordSist;
            end
            % -------------------------------------------------------------------------------------
            % Body fixed coordinate system defined by the principal axes of inertia
            % -------------------------------------------------------------------------------------
            % Inertia tensor of the marker cluster
            I11 =  TCS_Pos_Markers(2, :) * TCS_Pos_Markers(2, :)' ...
                +  TCS_Pos_Markers(3, :) * TCS_Pos_Markers(3, :)';
            I12 = - TCS_Pos_Markers(1, :) * TCS_Pos_Markers(2, :)';
            I13 = - TCS_Pos_Markers(1, :) * TCS_Pos_Markers(3, :)';
            I21 = I12;
            I22 =  TCS_Pos_Markers(1, :) * TCS_Pos_Markers(1, :)' ...
                +  TCS_Pos_Markers(3, :) * TCS_Pos_Markers(3, :)';
            I23 = - TCS_Pos_Markers(2, :) * TCS_Pos_Markers(3, :)';
            I31 = I13;
            I32 = I23;
            I33 =  TCS_Pos_Markers(1, :) * TCS_Pos_Markers(1, :)' ...
                + TCS_Pos_Markers(2, :) * TCS_Pos_Markers(2, :)';
            
            I = [I11, I12, I13; I21, I22, I23; I31, I32, I33];
            % Eigenvalues and Eigenvectors
%             disp([S.Name]);
            [eigvec, eigval] = eig(I);
            
            
            % Sort the axes with the rule: x-axis maximum moment of inertia, y-axis medium, z-axis minimum
            [AscOrd_eigval, AscOrd_Index] = sort([eigval(1, 1); eigval(2, 2); eigval(3, 3)]);
            DesOrd_Index  = flipud(AscOrd_Index);
            
            % Build the rotation matrix. The third vector of the rotation matrix is obtained as the cross
            % product of the two first. This is to avoid a left hand coordinate system that we could obtain from
            % changing the order of the columns of the rotation matrix.
            TCS_R_BodyCoordSist = [eigvec(:, DesOrd_Index(1)), eigvec(: , DesOrd_Index(2)), cross(eigvec(: , DesOrd_Index(1)), ...
                eigvec(: , DesOrd_Index(2)))];
            
            % Markers new local coordinates defined by the eigenvectors.
            BodyCoordSist_R_TCS = TCS_R_BodyCoordSist';
            BodyCoordSist_Pos_Marker = BodyCoordSist_R_TCS * TCS_Pos_Markers;
            
            % Note: Glob_R_TCS is the identity matrix. Both coord syst. are "parallel"
            Glob_R_BodyCoordSist = TCS_R_BodyCoordSist;
        end