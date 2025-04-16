classdef TRS_SENSOR < SENSOR
    %TRS_SENSOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        InitPoint % handle to POINT
        EndPoint  % handle to POINT 
    end
    
    methods
        function TS = TRS_SENSOR(Name,Type,InitPoint,EndPoint)
            TS = TS@SENSOR(Name,Type);
            TS.InitPoint = InitPoint;
            TS.EndPoint = EndPoint;
        end
        function Sensor = getSensorData(TS,q_t)
            NSamples = size(q_t,1);
            Sensor.Name1 = [TS.Name,'x'];
            Sensor.Name2 = [TS.Name,'y'];
            Sensor.Name3 = [TS.Name,'z'];
            for i=1:NSamples
                if TS.InitPoint.Fixed == 1
                    PIni = TS.InitPoint.GlobalCoord;
                else
                    PIni = [q_t(i,TS.InitPoint.PosInq);q_t(i,TS.InitPoint.PosInq+1);q_t(i,TS.InitPoint.PosInq+2)];
                end
                if TS.EndPoint.Fixed == 1
                    PEnd = TS.EndPoint.GlobalCoord;
                else
                    PEnd = [q_t(i,TS.EndPoint.PosInq);q_t(i,TS.EndPoint.PosInq+1);q_t(i,TS.EndPoint.PosInq+2)];
                end
                Sensor.Val(i,1) = PEnd(1) - PIni(1);
                Sensor.Val(i,2) = PEnd(2) - PIni(2);
                Sensor.Val(i,3) = PEnd(3) - PIni(3);
            end
            
        end
            
        
    end
    
end

