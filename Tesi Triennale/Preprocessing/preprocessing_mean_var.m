clc 
clear 
close all 

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

% Prendiamo la quinta colonna di Data e la diciassettesima colonna
% di Data. Ci muoviamo lungo la quinta. Ogni volta che l'elemento
% è compreso tra 30 e 31 escluso 31 aggiungiamo a un box tutti gli 
% elementi contenuti nell'elemento corrispondente alla diciassettesima 
% colonna, poi alla fine calcoliamo la media di tutti i numeri
% contenuti in questo box, e facciamo così anche per tutte le altre
% settimane

%baseFHR120bpm
media_sani = zeros(1, 3);
varianza_sani = zeros(1, 3);
for j = 1:10
bbox = []; 
for i = 1:size(table2array(Data(:,5)))
    if (table2array(Data(i,5)) >= 30 + 1*(j-1) && table2array(Data(i,5)) < 31 + 1*(j-1))
       bbox = [bbox; cell2mat(table2array(Data(i, 24)))]; 
    end
end
media_sani(j) = mean(bbox);
varianza_sani(j) = var(bbox);
% Apen_sani(j) = Apen; 
% SampEn_sani(j) = SampEn; 
% LZ2_sani(j) = LZ2; 
% LZ3_sani(j)= LZ3; 
end

% FATTO DA GIANLUCA
% Stessa cosa per Data_m
media_m = zeros(1, 3);
varianza_m = zeros(1, 3);

for j = 1:10
bbox = []; 
for i = 1:size(table2array(Data_m(:,5)))
    if (table2array(Data_m(i,5)) >= 30 + 1*(j-1) && table2array(Data_m(i,5)) < 31 + 1*(j-1))
        bbox = [bbox; cell2mat(table2array(Data_m(i, 24)))]; 
    end
end
media_m(j) = mean(bbox);
varianza_m(j) = var(bbox);
end

% FATTO DA GIANLUCA
figure(1)
hold on
plot(1:10, abs(media_sani - media_m), 'k-', 'LineWidth', 1.5)
plot(1:10, abs(media_sani - media_m), 'ko', 'MarkerSize', 15)
hold off
box on

%close all 

figure(2)
hold on
plot(1:10, media_sani, 'y-', 'LineWidth', 1.5, 'DisplayName', 'sani')
plot(1:10, media_sani, 'y.', 'MarkerSize', 10, 'DisplayName', 'sani')
plot(1:10, media_sani, 'ko', 'MarkerSize', 5, 'DisplayName', 'sani')
plot(1:10, media_m, 'b-',   'LineWidth', 1.5, 'DisplayName', 'IUGR')
plot(1:10, media_m, 'b.',   'MarkerSize', 10, 'DisplayName', 'IUGR')
plot(1:10, media_m, 'ko',   'MarkerSize', 5, 'DisplayName', 'IUGR')
hold off
xticks([1:10])
xticklabels({'30','31','32','33','34','35', '36', '37', '38', '39'})
box on

%close all 

figure(3)
hold on
plot(1:10, abs(varianza_sani - varianza_m), 'k-', 'LineWidth', 1.5)
plot(1:10, abs(varianza_sani - varianza_m), 'ko', 'MarkerSize', 15)
xticks([1:10])
xticklabels({'30','31','32','33','34','35', '36', '37', '38', '39'})
hold off
box on
%close all

%FHR120ms
media_sani = zeros(1,3);
varianza_sani = zeros(1, 3);
for j = 1:10
bbox = []; 
for i = 1:size(table2array(Data(:,5)))
    if (table2array(Data(i,5)) >= 30 + 1*(j-1) && table2array(Data(i,5)) < 31 + 1*(j-1))
        bbox = [bbox; cell2mat(table2array(Data(i, 22)))]; 
    end
end
media_sani(j) =  mean(bbox);
varianza_sani(j) = var(bbox);
end

media_m = zeros(1, 3);
varianza_m = zeros(1, 3);
for j = 1:10
bbox = []; 
for i = 1:size(table2array(Data_m(:,5)))
    if (table2array(Data_m(i,5)) >= 30 + 1*(j-1) && table2array(Data_m(i,5)) < 31 + 1*(j-1))
        bbox = [bbox; cell2mat(table2array(Data_m(i, 22)))];
    end
end
media_m(j) = mean(bbox);
varianza_m(j) = var(bbox);
end

figure(4)
hold on
plot(1:10, media_sani, 'y-', 'LineWidth', 1.5, 'DisplayName', 'sani')
plot(1:10, media_sani, 'y.', 'MarkerSize', 10, 'DisplayName', 'sani')
plot(1:10, media_sani, 'ko', 'MarkerSize', 5, 'DisplayName', 'sani')
plot(1:10, media_m, 'b-',   'LineWidth', 1.5, 'DisplayName', 'IUGR')
plot(1:10, media_m, 'b.',   'MarkerSize', 10, 'DisplayName', 'IUGR')
plot(1:10, media_m, 'ko',   'MarkerSize', 5, 'DisplayName', 'IUGR')
hold off
xticks([1:10])
xticklabels({'30','31','32','33','34','35', '36', '37', '38', '39'})
box on

%close all

figure(5)
hold on
plot(1:10, abs(media_sani - media_m), 'k-', 'LineWidth', 1.5)
plot(1:10, abs(media_sani - media_m), 'ko', 'MarkerSize', 15)
xticks([1:10])
xticklabels({'30','31','32','33','34','35', '36', '37', '38', '39'})
hold off
box on

%close all 

figure(6)
hold on
plot(1:10, abs(varianza_sani - varianza_m), 'k-', 'LineWidth', 1.5)
plot(1:10, abs(varianza_sani - varianza_m), 'ko', 'MarkerSize', 15)
xticks([1:10])
xticklabels({'30','31','32','33','34','35', '36', '37', '38', '39'})
hold off
box on
%close all

%FHR24mssenza
media_sani = zeros(1,3);
varianza_sani = zeros(1, 3);
for j = 1:10
bbox = []; 
for i = 1:size(table2array(Data(:,5)))
    if (table2array(Data(i,5)) >= 30 + 1*(j-1) && table2array(Data(i,5)) < 31 + 1*(j-1))
        bbox = [bbox; cell2mat(table2array(Data(i, 32)))]; 
    end
end
media_sani(j) =  mean(bbox);
varianza_sani(j) = var(bbox);
end

media_m = zeros(1, 3);
varianza_m = zeros(1, 3);
j = 0;
for j = 1:10
bbox = []; 
for i = 1:size(table2array(Data_m(:,5)))
    if (table2array(Data_m(i,5)) >= 30 + 1*(j-1) && table2array(Data_m(i,5)) < 31 + 1*(j-1))
        bbox = [bbox; cell2mat(table2array(Data_m(i, 32)))]; 
    end
end
media_m(j) = mean(bbox);
varianza_m(j) = var(bbox);
end

figure(7)
hold on
plot(1:10, media_sani, 'y-', 'LineWidth', 1.5, 'DisplayName', 'sani')
plot(1:10, media_sani, 'y.', 'MarkerSize', 10, 'DisplayName', 'sani')
plot(1:10, media_sani, 'ko', 'MarkerSize', 5, 'DisplayName', 'sani')
plot(1:10, media_m, 'b-',   'LineWidth', 1.5, 'DisplayName', 'IUGR')
plot(1:10, media_m, 'b.',   'MarkerSize', 10, 'DisplayName', 'IUGR')
plot(1:10, media_m, 'ko',   'MarkerSize', 5, 'DisplayName', 'IUGR')
hold off
xticks([1:10])
xticklabels({'30','31','32','33','34','35', '36', '37', '38', '39'})
box on

%close all

figure(8)
hold on
plot(1:10, abs(media_sani - media_m), 'k-', 'LineWidth', 1.5)
plot(1:10, abs(media_sani - media_m), 'ko', 'MarkerSize', 15)
xticks([1:10])
xticklabels({'30','31','32','33','34','35', '36', '37', '38', '39'})
hold off
box on

% % close all 

figure(9)
hold on
plot(1:10, abs(varianza_sani - varianza_m), 'k-', 'LineWidth', 1.5)
plot(1:10, abs(varianza_sani - varianza_m), 'ko', 'MarkerSize', 15)
xticks([1:10])
xticklabels({'30','31','32','33','34','35', '36', '37', '38', '39'})
hold off
box on
% % close all

% % % %CHIEDERE A GIULIO!! il problema qua è che le intacc e intdecc sono nx4
% % % %perciò continuo ad avere dei problemi. probabilmente devo riuscire a
% % % %convertire le "matrici" in vettori ma provando sia con ROWS2VECTORS che
% % % %con reshape ho avuto dei problemi.
%intaccelerazioni 
% % media_sani = zeros(1,3);
% % varianza_sani = zeros(1, 3);
% % for j = 1:10
% % bbox = []; 
% % for1 i = 1:size(table2array(Data(:,5)))
% %     if (table2array(Data(i,5)) >= 30 + 1*(j-1) && table2array(Data(i,5)) < 31 + 1*(j-1))
% %         bbox = [bbox; cell2mat(table2array(Data_m(i, 25)))]; 
% %     end
% % end
% % media_sani(j) =  mean(bbox);
% % varianza_sani(j) = var(bbox);
% % end
% % 
% % media_m = zeros(1, 3);
% % varianza_m = zeros(1, 3);
% % j = 0;
% % for j = 1:10
% % bbox = []; 
% % for i = 1:size(table2array(Data_m(:,5)))
% %     if (table2array(Data_m(i,5)) >= 30 + 1*(j-1) && table2array(Data_m(i,5)) < 31 + 1*(j-1))
% %         bbox = [bbox; cell2mat(table2array(Data_m(i, 25)))]; 
% %     end
% % end
% % media_m(j) = mean(bbox);
% % varianza_m(j) = var(bbox);
% % end


 
% % % %intdec
% % % media_sani = zeros(1,3);
% % % varianza_sani = zeros(1, 3);
% % % for j = 1:10
% % % bbox = []; 
% % % for i = 1:size(table2array(Data(:,5)))
% % %     if (table2array(Data(i,5)) >= 30 + 1*(j-1) && table2array(Data(i,5)) < 31 + 1*(j-1))
% % %         bbox = [bbox; cell2mat(table2array(Data_m(i, 25)))]; 
% % %     end
% % % end
% % % media_sani(j) =  mean(bbox);
% % % varianza_sani(j) = var(bbox);
% % % end
% % % 
% % % media_m = zeros(1, 3);
% % % varianza_m = zeros(1, 3);
% % % j = 0;
% % % for j = 1:10
% % % bbox = []; 
% % % for i = 1:size(table2array(Data_m(:,5)))
% % %     if (table2array(Data_m(i,5)) >= 30 + 1*(j-1) && table2array(Data_m(i,5)) < 31 + 1*(j-1))
% % %         bbox = [bbox; cell2mat(table2array(Data_m(i, 25)))]; 
% % %     end
% % % end
% % % media_m(j) = mean(bbox);
% % % varianza_m(j) = var(bbox);
% % % end
% % % 
% % % %intdec
% % % media_sani = zeros(1,3);
% % % varianza_sani = zeros(1, 3);
% % % for j = 1:10
% % % bbox = []; 
% % % for i = 1:size(table2array(Data(:,5)))
% % %     if (table2array(Data(i,5)) >= 30 + 1*(j-1) && table2array(Data(i,5)) < 31 + 1*(j-1))
% % %         bbox = [bbox; cell2mat(table2array(Data_m(i, 26)))]; 
% % %     end
% % % end
% % % media_sani(j) =  mean(bbox);
% % % varianza_sani(j) = var(bbox);
% % % end
% % % 
% % % media_m = zeros(1, 3);
% % % varianza_m = zeros(1, 3);
% % % j = 0;
% % % for j = 1:10
% % % bbox = []; 
% % % for i = 1:size(table2array(Data_m(:,5)))
% % %     if (table2array(Data_m(i,5)) >= 30 + 1*(j-1) && table2array(Data_m(i,5)) < 31 + 1*(j-1))
% % %         bbox = [bbox; cell2mat(table2array(Data_m(i, 26)))]; 
% % %     end
% % % end
% % % media_m(j) = mean(bbox);
% % % varianza_m(j) = var(bbox);
% % % end




