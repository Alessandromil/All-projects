function SensorIndex = getSensorIndex(Name,Sensors)

NSensors = size(Sensors,1);
for i=1:NSensors
    if strcmpi(Sensors{i}.Type,'SPH')
        if strcmpi(Name,Sensors{i}.Name)
            SensorIndex = i;
            return
        end
    end
end
error(['The sensor ',Name,' is not defined.'])