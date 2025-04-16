clear, clc, close all 

media_base = [];
for i = 1:size(Data, 1)
    base = Data.baseFHR24bpm(i); 
    media_base = [media_base, mean(base)];
end