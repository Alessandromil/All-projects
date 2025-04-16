% plot results using function compget y plot

ResultsPath = '\Experiments\ExpRightLeg\Results\'; % in Windows only

ResultsFile = 'ExpRightLeg_RightLeg_Group10_ankle1.sen';
%ResultsFile = 'ExpRightLeg_RightLeg_Group10_knee1.sen';

SensorName = 'Ankle_sensor'; % Ankle
Ankle_sensorData = compget([pwd,ResultsPath], ResultsFile, SensorName);
figure('Name','Ankle','NumberTitle','off','WindowStyle','docked'), hold on,
title(ResultsFile),
plot(Ankle_sensorData), 
legend('Dorsiflexion/Plantarflexion','Fixed parameter','Inversion/Eversion')
xlabel('Photograms');
ylabel('Degree');

maxValues = [];
minValues = [];
maxIndices = [];
minIndices = [];

[numRows, numCols] = size(Ankle_sensorData); 

for col = 1:numCols
    [maxValue, maxIndex] = max(Ankle_sensorData(:, col));  
    [minValue, minIndex] = min(Ankle_sensorData(:, col));  
    maxValues = [maxValues; maxValue];
    minValues = [minValues; minValue];
    maxIndices = [maxIndices; maxIndex];
    minIndices = [minIndices; minIndex];
   
    plot(maxIndex, maxValue, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5); 
    plot(minIndex, minValue, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);  
end

legend('Dorsiflexion/Plantarflexion', 'Fixed parameter', 'Inversion/Eversion', 'Maximum', 'Minimum');

SensorName = 'Knee_sensor'; % Knee
Knee_sensorData = compget([pwd,ResultsPath], ResultsFile, SensorName);
figure('Name','Knee','NumberTitle','off','WindowStyle','docked'), hold on,
title(ResultsFile),  
plot(Knee_sensorData),

xlabel('Photograms');
ylabel('Degree');

[minValue, minIndex] = min(Knee_sensorData); 
if ismatrix(Knee_sensorData)
    [minValue, minIndex] = min(Knee_sensorData(:));
    [row, col] = ind2sub(size(Knee_sensorData), minIndex); 
end
plot(minIndex, minValue, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5); 
legend('Flexion/Extension','Fixed parameter','Fixed parameter','minimum');



