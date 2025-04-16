clc 
clear 
close all 

malati = load('IUGR.mat'); 
Data_m = malati.IUGR; 
par = load('PARAMETRI_struct.mat'); 
NamesPar
par = par.par; 
parametri = cell2mat(struct2cell(par)); %cell2mat converte un cell array che contiene dati dello stesso tipo in una unica matrice.struct2cell converte una struct array in una cell array
[Data_m.baseFHR24bpm, Data_m.intacc24bpm, Data_m.intdec24bpm, Data_m.FHR120bpm, Data_m.FHR24bpm, Data_m.FHR120ms, Data_m.FHR24ms, Data_m.base120bpm, Data_m.intacc120bpm, Data_m.intdec120bpm]...
    = cellfun(@(x, y) PreProc([x', y'], parametri, 1), Data_m.FHR, Data_m.QUALITA, 'uni', 0);
%idx = cellfun(@(x) isempty(x), data.baseFHR24bpm); 
%data(idx, :) = []; 
[Data_m.grandiAcc, Data_m.piccoleAcc, Data_m.grandiDec, Data_m.piccoleDec] = cellfun(@(x, y) accDec(x, y), Data_m.intacc24bpm, Data_m.intdec24bpm, 'uni', 0); 


[Data_m.FHR24bpmsenza, Data_m.FHR24mssenza] = cellfun(@(a, b, c, d, e, f) excludeAccDec(a, b, c, d, e, f), Data_m.FHR24bpm, Data_m.FHR24ms, Data_m.grandiAcc, Data_m.piccoleAcc, Data_m.grandiDec, Data_m.piccoleDec, 'uni', 0); 