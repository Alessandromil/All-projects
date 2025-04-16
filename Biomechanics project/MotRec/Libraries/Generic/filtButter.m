function Data_f = filtButter(Data, sampleFreq, cutFreq)

if(~isempty(find(isnan(Data))))
    warning('NaN values not accepted in a trajectory!!');
    Data_f=Data;
    return
end

% size
nSamples = size(Data,1);
if nSamples <= 3*3  % Three times the filter order!!
    %warning('nSamples must be > 9 !');
    Data_f=Data;
    return
end

% Get parameters of Butterworth filter order 3
normFreq = cutFreq/(sampleFreq/2); % Shanon theorem
[B,A]    = butter(3,normFreq); 

% filtering in double direction to avoid desfase angular
Data_f = filtfilt(B,A,Data); % Double direction filtering. 

