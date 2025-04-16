clc 
clear 
close all 
sani = load('SANI (1).mat'); 
Data = sani.SANI; 
sani = load('SANI (1).mat'); 
Data = sani.SANI; 
par = load("PARAMETRI_struct.mat");
NamesPar
par = par.par; 
parametri = cell2mat(struct2cell(par)); %cell2mat converte un cell array che contiene dati dello stesso tipo in una unica matrice.struct2cell converte una struct array in una cell array
[Data.baseFHR24bpm, Data.intacc24bpm, Data.intdec24bpm, Data.FHR120bpm, Data.FHR24bpm, Data.FHR120ms, Data.FHR24ms, Data.base120bpm, Data.intacc120bpm, Data.intdec120bpm]...
    = cellfun(@(x, y) PreProc([x', y'], parametri, 1), Data.FHR, Data.QUALITA, 'uni', 0);
idx = cellfun(@(x) isempty(x), Data.baseFHR24bpm); 
Data(idx, :) = []; 
[Data.grandiAcc, Data.piccoleAcc, Data.grandiDec, Data.piccoleDec] = cellfun(@(x, y) accDec(x, y), Data.intacc24bpm, Data.intdec24bpm, 'uni', 0); 


[Data.FHR24bpmsenza, Data.FHR24mssenza] = cellfun(@(a, b, c, d, e, f) excludeAccDec(a, b, c, d, e, f), Data.FHR24bpm, Data.FHR24ms, Data.grandiAcc, Data.piccoleAcc, Data.grandiDec, Data.piccoleDec, 'uni', 0); 

malati = load('IUGR.mat'); 
Data_m = malati.IUGR;
[Data_m.baseFHR24bpm, Data_m.intacc24bpm, Data_m.intdec24bpm, Data_m.FHR120bpm, Data_m.FHR24bpm, Data_m.FHR120ms, Data_m.FHR24ms, Data_m.base120bpm, Data_m.intacc120bpm, Data_m.intdec120bpm]...
    = cellfun(@(x, y) PreProc([x', y'], parametri, 1), Data_m.FHR, Data_m.QUALITA, 'uni', 0);
idx = cellfun(@(x) isempty(x), Data_m.baseFHR24bpm); 
Data_m(idx, :) = []; 
[Data_m.grandiAcc, Data_m.piccoleAcc, Data_m.grandiDec, Data_m.piccoleDec] = cellfun(@(x, y) accDec(x, y), Data_m.intacc24bpm, Data_m.intdec24bpm, 'uni', 0); 


[Data_m.FHR24bpmsenza, Data_m.FHR24mssenza] = cellfun(@(a, b, c, d, e, f) excludeAccDec(a, b, c, d, e, f), Data_m.FHR24bpm, Data_m.FHR24ms, Data_m.grandiAcc, Data_m.piccoleAcc, Data_m.grandiDec, Data_m.piccoleDec, 'uni', 0);

baseline = Data{1,"base120bpm"}{1,1};
% 
% [nr_campioni_1, ~] = size(baseline);
% nr_secondi = nr_campioni_1 / 2; 
% fs = 2;
% 
% t=0:1/fs:(nr_campioni_1 - 1)/fs; %definisco il vettore tempo così: va da 0, con passo 1/fs, al numero di secondi meno uno perchè è necessario che il numero di campioni e il vettore tempo abbiano la stessa dimensione.

% figure(1)
% plot(t, baseline(:, 1))
% ylabel('bpm')
% xlabel('tempo [s]')
% %title('FHR 120 bpm')
% box on 

segnale_FHR = Data{1, "FHR120bpm"}{1,1}(:,1);
[nr_campioni, ~] = size(segnale_FHR); %dipende dalla durata dell'acquisizione e dalla fs  
fs = 2;

t_1=0:1/fs:(nr_campioni - 1)/fs;

%SEGNALE FHR 120 BPM
figure(1) 
plot(t_1, segnale_FHR(:, 1)) %unione di valori discreti consecutivi tra di loro
xlabel('tempo [s]')
ylabel('bpm')
%title('FHR 120 bpm')
box on 

%SEGNALE FHR + BASELINE
figure(2)
hold on 
plot(t_1, segnale_FHR(:, 1)) %unione di valori discreti consecutivi tra di loro
plot(t_1, baseline(:, 1), 'LineWidth', 2)
xlabel('tempo [s]')
ylabel('bpm')
hold off
box on

%FINESTRA TRE MINUTI, SEGNALE 120 bpm
N120 = 3*120; 
N3min = floor(length(segnale_FHR)/N120);
spezzone = segnale_FHR(1:N120); 
t_spezzone = 0:1/fs:(length(spezzone) - 1)/fs; 

figure(3)
hold on 
plot(t_1, segnale_FHR)
plot(t_spezzone, spezzone)
xlabel('tempo [s]')
ylabel('bpm')
hold off
box on 
 
%FINESTRA 1 MINUTO
FHR24bpm = Data{1, "FHR24bpm"}{1,1}(:,1);
[nr_campioni_1, ~] = size(FHR24bpm); %dipende dalla durata dell'acquisizione e dalla fs  
fs = 0.4;
N1min = floor(length(FHR24bpm)/24);
t=0:1/fs:(nr_campioni_1 - 1)/fs;

spezzone = FHR24bpm(1:24);
t_spezzone = 0:1/fs:(length(spezzone) - 1)/fs; 

figure(4)
hold on 
plot(t, FHR24bpm)
plot(t_spezzone,spezzone)
xlabel('tempo [s]')
ylabel('bpm')
box on 

segnale_FHR = Data{1, "FHR120ms"}{1,1}(:,1);
[nr_campioni, ~] = size(segnale_FHR); %dipende dalla durata dell'acquisizione e dalla fs  
fs = 2;

t_1=0:1/fs:(nr_campioni - 1)/fs;
figure(1) 
plot(t_1, segnale_FHR(:, 1)) %unione di valori discreti consecutivi tra di loro
xlabel('tempo [s]')
ylabel('bpm')
%title('FHR 120 bpm')
box on 