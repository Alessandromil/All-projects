classdef SPH_SENSOR < SENSOR
    %SPH_SENSOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Segment1 % Handle to a Segment 1                            SEGMENT
        Segment2 % Handle to a Segment 2                            SEGMENT
        RotSeq   % Rotation secuence to measure angles              double
        Perm1x   % Vector x to rotate from Previos Segment to Joint double(3x1)
        Perm1y   % Vector y to rotate from Previos Segment to Joint double(3x1)
        Perm2x   % Vector x to rotate from Joint to Next Segment    double(3x1)
        Perm2y   % Vector y to rotate from Joint to Next Segment    double(3x1)
        Parent_R_Joint  % Fixed rotation between Parent and Joint   double(3x3)
        Child_R_Joint   % Fixed rotation between Child and Joint    double(3x3)
    end
    
    methods
        function SS = SPH_SENSOR(Name,Type,Segment1,Segment2,RotSeq,Perm1x,Perm1y,Perm2x,Perm2y)
            SS = SS@SENSOR(Name,Type);
            SS.Segment1 = Segment1;
            SS.Segment2 = Segment2;
            SS.RotSeq = RotSeq;
            SS.Perm1x = Perm1x;
            SS.Perm1y = Perm1y;
            SS.Perm2x = Perm2x;
            SS.Perm2y = Perm2y;
        end
        function checkSubjectParsData(SS,File)
            if ~isempty(SS.Perm1x)
                if ischar(SS.Perm1x)
                    error(['The value of ',SS.Perm1x,' for Perm1 in sensor ',SS.Name,' is not in File: ',File])
                end
            end
            if ~isempty(SS.Perm2x)
                if ischar(SS.Perm2x)
                    error(['The value of ',SS.Perm2x,' for Perm2 in sensor ',SS.Name,' is not in File: ',File])
                end
            end
        end
        function Sensor = getSensorData(SS,q_t)
            % 1- Calc Permutation Matrices
            % calculate axes coordinates
            if ~isempty(SS.Perm1x)
                Xaxis = SS.Perm1x;
                Yaxis = SS.Perm1y;
                Xaxis = Xaxis/norm(Xaxis);
                Yaxis = Yaxis/norm(Yaxis);
                Zaxis = cross(Xaxis,Yaxis);
                Zaxis = Zaxis/norm(Zaxis);
                Yaxis = cross(Zaxis,Xaxis);
                Yaxis = Yaxis/norm(Yaxis);
                Xaxis = cross(Yaxis,Zaxis);
                Xaxis = Xaxis/norm(Xaxis);
                % Fixed Seg1_R_Joint
                Perm1 = [Xaxis,Yaxis,Zaxis];    
            else
                Perm1 = eye(3);
                SS.Perm1x = [1;0;0];
                SS.Perm1y = [0;1;0];
            end
            SS.Parent_R_Joint = Perm1;
            if ~isempty(SS.Perm2x)
                Xaxis = SS.Perm2x;
                Yaxis = SS.Perm2y;
                Xaxis = Xaxis/norm(Xaxis);
                Yaxis = Yaxis/norm(Yaxis);
                Zaxis = cross(Xaxis,Yaxis);
                Zaxis = Zaxis/norm(Zaxis);
                Yaxis = cross(Zaxis,Xaxis);
                Yaxis = Yaxis/norm(Yaxis);
                Xaxis = cross(Yaxis,Zaxis);
                Xaxis = Xaxis/norm(Xaxis);
                % Fixed Seg2_R_Joint
                Perm2 = [Xaxis,Yaxis,Zaxis];
            else
                Perm2 = eye(3);
                SS.Perm2x = [1;0;0];
                SS.Perm2y = [0;1;0];
            end
            SS.Child_R_Joint = Perm2;
            if strcmp(SS.RotSeq,'123')
                Sensor.Name1 = [SS.Name,'_Alpha'];
                Sensor.Name2 =  [SS.Name,'_Beta'];
                Sensor.Name3 = [SS.Name,'_Gamma'];
            elseif strcmp(SS.RotSeq,'321')
                Sensor.Name1 = [SS.Name,'_Gamma'];
                Sensor.Name2 =  [SS.Name,'_Beta'];
                Sensor.Name3 = [SS.Name,'_Alpha'];
            else
                error([SS.RotSeq,' is not a type of rotation sequence in ',SS.Name,' sensor.'])
            end
            NSamples = size(q_t,1);
            for i=1:NSamples
                % 2- Calc Glob_R_Seg1 & Glob_R_Seg2
                Glob_R_Seg1 = SS.Segment1.getRd(q_t(i,:));
                Glob_R_Seg2 = SS.Segment2.getRd(q_t(i,:));
                % 3- Calc relative rotation matrix
                Joint_R_Joint = Perm1'* Glob_R_Seg1' * Glob_R_Seg2 * Perm2;
                % 4 Calculate the euler angles from de relative rotation matrix
                if strcmp(SS.RotSeq,'123')
                    if i == 1
                        angelref = 0;
                        Sensor.Val(i,1) =  0;
                        Sensor.Val(i,2) =  0;
                        Sensor.Val(i,3) =  0;
                    else
                       angelref = Sensor.Val(i-1,3); 
%                         angref = [Sensor.Val(i-1,1);Sensor.Val(i-1,2);Sensor.Val(i-1,3)];  % values in previous frame
%                     end
%                     [Sensor.Val(i,1),Sensor.Val(i,2),Sensor.Val(i,3)] = SS.weuler_xyz(Joint_R_Joint,angref);
                    end
                    [Sensor.Val(i,1),Sensor.Val(i,2),Sensor.Val(i,3)]= rot123s(Joint_R_Joint,angelref);
                elseif strcmp(SS.RotSeq,'321')
%                     [Sensor.Val(i,1),Sensor.Val(i,2),Sensor.Val(i,3)] = rot321s(Joint_R_Joint);
                    if i == 1
                        angref = [0; 0; 0];
                    else
                        angref = [Sensor.Val(i-1,1);Sensor.Val(i-1,2);Sensor.Val(i-1,3)];  % values in previous frame
                    end
                    [Sensor.Val(i,1),Sensor.Val(i,2),Sensor.Val(i,3)] = SS.weuler_zyx(Joint_R_Joint,angref);
                end
            end
            % 5 Pass angles from radians to grades
            for i=1:NSamples
                for j=1:3
                    Sensor.Val(i,j) = Sensor.Val(i,j)*180/pi;
                end
            end
            
            
        end
        function [psi,teta,phi] = weuler_xyz(SS,r,angref)
            % [psi,teta,phi] = weuler_xyz(r,angref);
            %
            % fonction permettant de calculer les trois angles d'euler suivant la
            % sequence de rotation rx(psi), ry(teta), rz(phi) (axes mobiles)
            
            %   Gilles Monnier, INRETS, Altran pour Renault
            %   Implementation:
            %       Matlab v6.5 r13
            %   Toolboxes:
            %       general
            %		Dépendance externe:
            %     aucune
            %   Historique des revisions majeures:
            %       Rev 1.0 15/10/2004 : creation
            
            epsi=1.e-6;
            
            % temporarily disconnects atan2 warning. The warning is reconnected at the
            % end of the function.
%             warning('off','MATLAB:atan2:complexArgument')
            
            
            %---valeurs possibles pour teta
            st2=r(1,3);
            if abs(st2)<epsi
                st2=0;
                teta2(1)=0;
                teta2(2)=pi;
            elseif abs(abs(st2)-1)<epsi
                teta2(1)=pi/2;
                teta2(2)=-pi/2;
            else
                teta2(1)=asin(st2);
                q=abs(st2);
                aa=st2/q;
                teta2(2)=pi*aa-teta2(1);
            end
            
            %---on calcule les solutions qui en decoulent
            for i=1:2
                ct2=cos(teta2(i));
                if abs(ct2)<epsi %indétermination entre psi et phi
                    phi2(i)=0;
                    st1=r(3,2);
                    ct1=r(2,2);
                    if abs(st1)<epsi
                        psi2(i)=(1-ct1)*pi/2;
                    elseif abs(ct1)<epsi
                        psi2(i)=st1*pi/2;
                    else
                        psi2(i)=atan2(st1,ct1);
                    end
                else
                    % ------------       pour psi
                    st1= -r(2,3)/ct2;
                    ct1=r(3,3)/ct2;
                    if abs(st1)<=epsi
                        psi2(i)=(1-ct1)*pi/2;
                    elseif abs(ct1)<=epsi
                        psi2(i)=st1*pi/2;
                    else
                        psi2(i)=atan2(st1,ct1);
                    end
                    %------------        et pour phi
                    st3=-r(1,2)/ct2;
                    ct3=r(1,1)/ct2;
                    if abs(st3)<epsi
                        phi2(i)=(1-ct3)*pi/2;
                    elseif abs(ct3)<epsi
                        phi2(i)=st3*pi/2;
                    else
                        phi2(i)=atan2(st3,ct3);
                    end
                end
            end
            
            %----------------------on retient la solution du moindre
            %                         deplacement angulaire global
            for i=1:2
                D3=min([abs(phi2(i)-angref(3)) abs(phi2(i)-angref(3)-2*pi) abs(phi2(i)-angref(3)+2*pi)]);
                D2=min([abs(teta2(i)-angref(2)) abs(teta2(i)-angref(2)-2*pi) abs(teta2(i)-angref(2)+2*pi)]);
                D1=min([abs(psi2(i)-angref(1)) abs(psi2(i)-angref(1)-2*pi) abs(psi2(i)-angref(1)+2*pi)]);
                
                som(i)=D1+D2+D3;
                %abs(phi2(i)-angref(3))+abs(teta2(i)-angref(2))+...
                %abs(psi2(i)-angref(1));
            end
            isol=1;
            if som(2)<som(1)
                isol=2;
            end
            phi=phi2(isol);
            teta=teta2(isol);
            psi=psi2(isol);
            
            % Re-connects atan2 warning for complex numbers.
%             warning('on','MATLAB:atan2:complexArgument')
        end
        function [psi,teta,phi]=  weuler_zyx(SS,r,angref)
            % [psi,teta,phi] = weuler_xyz(r,angref);
            %
            % fonction permettant de calculer les trois angles d'euler suivant la
            % sequence de rotation rz(psi), ry(teta), rx(phi) (axes mobiles)
            %   Gilles Monnier, Altran pour Renault
            %   Implementation:
            %       Matlab v6.5 r13
            %   Toolboxes:
            %       general
            %		Dépendance externe:
            %     aucune
            %   Historique des revisions majeures:
            %       Rev 1.0 23/08/2006 : creation
            epsi=1.e-12;
%             epsi=1.e-6; Problemas de precisión
            
            % temporarily disconnects atan2 warning. The warning is reconnected at the
            % end of the function.
%             warning('off','MATLAB:atan2:complexArgument')
            
            
            %---valeurs possibles pour teta
            st2=-r(3,1);
            if abs(st2)<epsi
                st2=0;
                teta2(1)=0;
                teta2(2)=pi;
            elseif abs(abs(st2)-1)<epsi
                teta2(1)=pi/2;
                teta2(2)=-pi/2;
            else
                teta2(1)=asin(st2);
                q=abs(st2);
                aa=st2/q;
                teta2(2)=pi*aa-teta2(1);
            end
            
            %---on calcule les solutions qui en decoulent
            for i=1:2
                ct2=cos(teta2(i));
                if abs(ct2)<epsi
                    phi2(i)=0;
                    st1=-r(1,2);
                    ct1=r(2,2);
                    if abs(st1)<epsi
                        psi2(i)=(1-ct1)*pi/2;
                    elseif abs(ct1)<epsi
                        psi2(i)=st1*pi/2;
                    else
                        psi2(i)=atan2(st1,ct1);
                    end
                else
                    % ------------       pour psi
                    st1= r(2,1)/ct2;
                    ct1=r(1,1)/ct2;
                    if abs(st1)<=epsi
                        psi2(i)=(1-ct1)*pi/2;
                    elseif abs(ct1)<=epsi
                        psi2(i)=st1*pi/2;
                    else
                        psi2(i)=atan2(st1,ct1);
                    end
                    %------------        et pour phi
                    st3=r(3,2)/ct2;
                    ct3=r(3,3)/ct2;
                    if abs(st3)<epsi
                        phi2(i)=(1-ct3)*pi/2;
                    elseif abs(ct3)<epsi
                        phi2(i)=st3*pi/2;
                    else
                        phi2(i)=atan2(st3,ct3);
                    end
                end
            end
            
            %----------------------on retient la solution du moindre
            %                         deplacement angulaire global
            for i=1:2
                % 	% ORIGINAL
                %     D3=min([abs(phi2(i)-angref(3)) abs(phi2(i)-angref(3)-2*pi) abs(phi2(i)-angref(3)+2*pi)]);
                % 	D2=min([abs(teta2(i)-angref(2)) abs(teta2(i)-angref(2)-2*pi) abs(teta2(i)-angref(2)+2*pi)]);
                % 	D1=min([abs(psi2(i)-angref(1)) abs(psi2(i)-angref(1)-2*pi) abs(psi2(i)-angref(1)+2*pi)]);
                % 	som(i)=D1+D2+D3;
                % MODIFIED
                som(i)= abs(phi2(i)-angref(3)) + abs(teta2(i)-angref(2)) + abs(psi2(i)-angref(1));
            end
            isol=1;
            if som(2)<som(1)
                isol=2;
            end
            
            phi=phi2(isol);
            teta=teta2(isol);
            psi=psi2(isol);
            
            % Re-connects atan2 warning for complex numbers.
%             warning('on','MATLAB:atan2:complexArgument')
%             
%             return,
        end
    end
    
end

