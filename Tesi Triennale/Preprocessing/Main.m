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

N24 = 24*3; 
N120 = 120*3;

% Prendiamo la quinta colonna di Data e la diciassettesima colonna
% di Data. Ci muoviamo lungo la quinta. Ogni volta che l'elemento
% è compreso tra 30 e 31 escluso 31 aggiungiamo a un box tutti gli 
% elementi contenuti nell'elemento corrispondente alla diciassettesima 
% colonna, poi alla fine calcoliamo la media di tutti i numeri
% contenuti in questo box, e facciamo così anche per tutte le altre
% settimane
%% costruisco table dei parametri
ParametriSani = Data(:, [1 3 5]);
ParametriMalati = Data_m(:, [1 3 5]); 

for i = 1:size(Data,1)
    segnale = Data{i,"FHR120bpm"}{1,1}(:,1);
    baseline = Data{i,"base120bpm"}{1,1};
    qualita = Data{i,"QUALITA"}{1,1}';
    FHR24bpmsenza = Data{i,"FHR24bpmsenza"}{1,1};
  
    N3min = floor(length(segnale)/N120);
    N1min = floor(length(FHR24bpmsenza)/24);
    Apen = [];
    SampEn = [];
    LZ2 = [];
    LZ3 = [];
    STV = [];
    for n = 1:N3min
        spezzone = segnale(1+(n-1)*N120:n*N120);
        spezzoneQualita = qualita(1+(n-1)*N120:n*N120);   
        %controllo qualità
        if sum(spezzoneQualita>64)<(0.05*N120) %meno del 5% di interpolati
            [Apen(n),SampEn(n)] = apsampen(spezzone,2,0.15,1);
            [LZ2(n),LZ3(n)] = ComputeLZ(spezzone);
            %stima spettro
            spezzone = spezzone-mean(spezzone);
            xdft = fft(spezzone); %mappiamo il segnale nel dominio nel tempo nel dominio delle frequenze
            xdft = xdft(1:N120/2+1);
            psdx = (1/(2*N120)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1); %il primo campione è in continua. 
            freq = 0:2/length(spezzone):1;
            
            PTot = var(spezzone);
            LF(n) = (0.5*sum(psdx(freq>0.03 & freq<0.15)))/PTot;
            MF(n) = (0.5*sum(psdx(freq>0.15 & freq<0.5)))/PTot;
            HF(n) = (0.5*sum(psdx(freq>0.5 & freq<1)))/PTot;

        else
            [Apen(n),SampEn(n),LZ2(n),LZ3(n),LF(n),MF(n),HF(n)] = deal(nan);
        end
    end

    for n = 1:N1min
    spezzone24 = FHR24bpmsenza(1+(n-1)*24:n*24);
    spezzoneQualita = qualita(1+(n-1)*120:n*120); %occhio che qualità in120
    %controllo qualità
    if sum(spezzoneQualita>64)<(0.05*120) %meno del 5% di interpolati
        [STV(n), ~, ~, ~]=STV_II_m(spezzone24);
    else
        STV(n) = nan;
    end

end

    ParametriSani{i,"media"} = mean(baseline);
    ParametriSani{i,"Apen"} = mean(Apen(~isnan(Apen)));
    ParametriSani{i,"SampEn"} = mean(SampEn(~isnan(SampEn)));
    
    ParametriSani{i,"STV"} = mean(STV(~isnan(STV)));
    ParametriSani{i,"LZ2"} = mean(LZ2(~isnan(LZ2)));
    ParametriSani{i,"LZ3"} = mean(LZ3(~isnan(LZ3)));

    ParametriSani{i,"LF"} = mean(LF(~isnan(LF)));
    ParametriSani{i,"MF"} = mean(MF(~isnan(MF)));
    ParametriSani{i,"HF"} = mean(HF(~isnan(HF)));
end

for i = 1:size(Data_m,1)
    segnale = Data_m{i,"FHR120bpm"}{1,1}(:,1);
    baseline = Data_m{i,"base120bpm"}{1,1};
    qualita = Data_m{i,"QUALITA"}{1,1}';
    FHR24bpmsenza = Data_m{i,"FHR24bpmsenza"}{1,1};
  
    N3min = floor(length(segnale)/N120);
    N1min = floor(length(FHR24bpmsenza)/24);
    Apen = [];
    SampEn = [];
    LZ2 = [];
    LZ3 = [];
    STV = [];
    for n = 1:N3min
        spezzone = segnale(1+(n-1)*N120:n*N120);
        spezzoneQualita = qualita(1+(n-1)*N120:n*N120);   
        %controllo qualità
        if sum(spezzoneQualita>64)<(0.05*N120) %meno del 5% di interpolati
            [Apen(n),SampEn(n)] = apsampen(spezzone,2,0.15,1);
            [LZ2(n),LZ3(n)] = ComputeLZ(spezzone);
            %stima spettro
            spezzone = spezzone-mean(spezzone);
            xdft = fft(spezzone);
            xdft = xdft(1:N120/2+1);
            psdx = (1/(2*N120)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            freq = 0:2/length(spezzone):1;
            
            PTot = var(spezzone);
            LF(n) = (0.5*sum(psdx(freq>0.03 & freq<0.15)))/PTot;
            MF(n) = (0.5*sum(psdx(freq>0.15 & freq<0.5)))/PTot;
            HF(n) = (0.5*sum(psdx(freq>0.5 & freq<1)))/PTot;

        else
            [Apen(n),SampEn(n),LZ2(n),LZ3(n),LF(n),MF(n),HF(n)] = deal(nan);
        end
    end

    for n = 1:N1min
    spezzone24 = FHR24bpmsenza(1+(n-1)*24:n*24);
    spezzoneQualita = qualita(1+(n-1)*120:n*120); %occhio che qualità in120
    %controllo qualità
    if sum(spezzoneQualita>64)<(0.05*120) %meno del 5% di interpolati
        [STV(n), ~, ~, ~]=STV_II_m(spezzone24);
    else
        STV(n) = nan;
    end

end

    ParametriMalati{i,"media"} = mean(baseline);
    ParametriMalati{i,"Apen"} = mean(Apen(~isnan(Apen)));
    ParametriMalati{i,"SampEn"} = mean(SampEn(~isnan(SampEn)));
    
    ParametriMalati{i,"STV"} = mean(STV(~isnan(STV)));
    ParametriMalati{i,"LZ2"} = mean(LZ2(~isnan(LZ2)));
    ParametriMalati{i,"LZ3"} = mean(LZ3(~isnan(LZ3)));

    ParametriMalati{i,"LF"} = mean(LF(~isnan(LF)));
    ParametriMalati{i,"MF"} = mean(MF(~isnan(MF)));
    ParametriMalati{i,"HF"} = mean(HF(~isnan(HF)));
end
        %ho preso i patnum e li ho messi in un vettore
        %ho confrontato ogni patnum con i patnum della tabella che avevamo già completato. 
        %se il patnum della tabella è uguale a quello del vettore (il vettore è ottenuto tramite unique quindi tutti i 
        % patnum compaiono una sola volta)
        %si controlla ogni data di monitoraggio: in base all'intervallo in cui cade 
        % ho inserito la media del paziente per la specifica data in una scatola
        %che si svuota ad ogni nuovo patnum. quando il patnum corrente è
        %stato controllato su tutta la tabella e dopo che ho inserito tutte
        %le medie all'interno delle scatole relative all'intervallo, osservo se il contenuto della
        %scatola è maggiore di 1: in questo caso il contenuto della scatola
        %viene mediato. ogni scatola contiene le medie della specifica
        %settimana e ogni volta che il patnum è stato controllato su tutta
        %la tabella, svuoto la scatola in quella non corrente, ovvero in
        %quella che alla fine inserisco nella tabella. Sono partita
        %dall'assunto che per ogni paziente c'è almeno una registrazione
        %in ogni bin. 
    numero_identificativoS = unique(ParametriSani.patnum); 
    media_30_33 = []; 
    media_34_36 = []; 
    media_37_38 = [];
    for i = 1:length(numero_identificativoS) 
        pat_num = numero_identificativoS(i); 
        scatola_1 = []; 
        scatola_2 = []; 
        scatola_3 = []; 
        for n = 1:size(ParametriSani, 1) 
            if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 33 && ParametriSani.sett_gestazione(n) >= 30
               scatola_1 = [scatola_1, ParametriSani.media(n)]; 
            end
            if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 36 && ParametriSani.sett_gestazione(n) >= 34
               scatola_2 = [scatola_2, ParametriSani.media(n)];
            end
            if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 38 && ParametriSani.sett_gestazione(n) >= 37
               scatola_3 = [scatola_3, ParametriSani.media(n)];
            end
        end
            if length(scatola_1) > 1
               media_30_33 = [media_30_33, mean(scatola_1)]; 
            else 
                media_30_33 = [media_30_33, scatola_1];
            end
            if length(scatola_2) > 1
               media_34_36 = [media_34_36, mean(scatola_2)]; 
            else 
                media_34_36 = [media_34_36, scatola_2];
            end
            if length(scatola_3) > 1
               media_37_38 = [media_37_38, mean(scatola_3)]; 
            else 
                media_37_38 = [media_37_38, scatola_3];
            end
         end
         media_30_33 = media_30_33'; 
         media_34_36 = media_34_36'; 
         media_37_38 = media_37_38'; 
         media_sani = table(numero_identificativoS, media_30_33, media_34_36, media_37_38);
   
         %malati
         numero_identificativoM = unique(ParametriMalati.patnum);
    media_30_33 = []; 
    media_34_36 = []; 
    media_37_38 = [];
    for i = 1:length(numero_identificativoM) 
        pat_num = numero_identificativoM(i); 
        scatola_1 = []; 
        scatola_2 = []; 
        scatola_3 = []; 
        for n = 1:size(ParametriMalati, 1) 
            if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 33 && ParametriMalati.sett_gestazione(n) >= 30
               scatola_1 = [scatola_1, ParametriMalati.media(n)]; 
            end
            if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 36 && ParametriMalati.sett_gestazione(n) >= 34
               scatola_2 = [scatola_2, ParametriMalati.media(n)];
            end
            if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 38 && ParametriMalati.sett_gestazione(n) >= 37
               scatola_3 = [scatola_3, ParametriMalati.media(n)];
            end
        end
            if length(scatola_1) > 1
               media_30_33 = [media_30_33, mean(scatola_1)]; 
            else 
                media_30_33 = [media_30_33, scatola_1];
            end
            if length(scatola_2) > 1
               media_34_36 = [media_34_36, mean(scatola_2)]; 
            else 
                media_34_36 = [media_34_36, scatola_2];
            end
            if length(scatola_3) > 1
               media_37_38 = [media_37_38, mean(scatola_3)]; 
            else 
                media_37_38 = [media_37_38, scatola_3];
            end
         end
         media_30_33 = media_30_33'; 
         media_34_36 = media_34_36'; 
         media_37_38 = media_37_38'; 
         media_malati = table(numero_identificativoM, media_30_33, media_34_36, media_37_38);
 
 %faccio la stessa cosa con gli altri parametri 
 %STV
 STV_30_33 = []; 
 STV_34_36 = []; 
 STV_37_38 = []; 
 for i = 1:length(numero_identificativoS) 
        pat_num = numero_identificativoS(i); 
        scatola_1 = []; 
        scatola_2 = []; 
        scatola_3 = []; 
        for n = 1:size(ParametriSani, 1) 
            if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 33 && ParametriSani.sett_gestazione(n) >= 30
               scatola_1 = [scatola_1, ParametriSani.STV(n)]; 
            end
            if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 36 && ParametriSani.sett_gestazione(n) >= 34
               scatola_2 = [scatola_2, ParametriSani.STV(n)];
            end
            if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 38 && ParametriSani.sett_gestazione(n) >= 37
               scatola_3 = [scatola_3, ParametriSani.STV(n)];
            end
        end
            if length(scatola_1) > 1
               STV_30_33 = [STV_30_33, mean(scatola_1)]; 
            else 
                STV_30_33 = [STV_30_33, scatola_1];
            end
            if length(scatola_2) > 1
               STV_34_36 = [STV_34_36, mean(scatola_2)]; 
            else 
                STV_34_36 = [STV_34_36, scatola_2];
            end
            if length(scatola_3) > 1
               STV_37_38 = [STV_37_38, mean(scatola_3)]; 
            else 
                STV_37_38 = [STV_37_38, scatola_3];
            end
         end
       
         STV_30_33 = STV_30_33'; 
         STV_34_36 = STV_34_36'; 
         STV_37_38 = STV_37_38'; 
         STV_sani = table(numero_identificativoS, STV_30_33, STV_34_36, STV_37_38);

         STV_30_33 = []; 
         STV_34_36 = []; 
         STV_37_38 = []; 
        for i = 1:length(numero_identificativoM) 
            pat_num = numero_identificativoM(i); 
            scatola_1 = []; 
            scatola_2 = []; 
            scatola_3 = []; 
        for n = 1:size(ParametriMalati, 1) 
            if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 33 && ParametriMalati.sett_gestazione(n) >= 30
               scatola_1 = [scatola_1, ParametriMalati.STV(n)]; 
            end
            if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 36 && ParametriMalati.sett_gestazione(n) >= 34
               scatola_2 = [scatola_2, ParametriMalati.STV(n)];
            end
            if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 38 && ParametriMalati.sett_gestazione(n) >= 37
               scatola_3 = [scatola_3, ParametriMalati.STV(n)];
            end
        end
            if length(scatola_1) > 1
               STV_30_33 = [STV_30_33, mean(scatola_1)]; 
            else 
                STV_30_33 = [STV_30_33, scatola_1];
            end
            if length(scatola_2) > 1
               STV_34_36 = [STV_34_36, mean(scatola_2)]; 
            else 
                STV_34_36 = [STV_34_36, scatola_2];
            end
            if length(scatola_3) > 1
               STV_37_38 = [STV_37_38, mean(scatola_3)]; 
            else 
                STV_37_38 = [STV_37_38, scatola_3];
            end
         end
         STV_30_33 = STV_30_33'; 
         STV_34_36 = STV_34_36'; 
         STV_37_38 = STV_37_38';
         STV_malati = table(numero_identificativoM, STV_30_33, STV_34_36, STV_37_38);

         %ApEn
         ApEn_30_33 = []; 
         ApEn_34_36 = []; 
         ApEn_37_38 = []; 
         for i = 1:length(numero_identificativoS) 
            pat_num = numero_identificativoS(i); 
            scatola_1 = []; 
            scatola_2 = []; 
            scatola_3 = []; 
            for n = 1:size(ParametriSani, 1) 
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 33 && ParametriSani.sett_gestazione(n) >= 30
                   scatola_1 = [scatola_1, ParametriSani.Apen(n)]; 
                end
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 36 && ParametriSani.sett_gestazione(n) >= 34
                   scatola_2 = [scatola_2, ParametriSani.Apen(n)];
                end
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 38 && ParametriSani.sett_gestazione(n) >= 37
                   scatola_3 = [scatola_3, ParametriSani.Apen(n)];
                end
            end
                if length(scatola_1) > 1
                   ApEn_30_33 = [ApEn_30_33, mean(scatola_1)]; 
                else 
                    ApEn_30_33 = [ApEn_30_33, scatola_1];
                end
                if length(scatola_2) > 1
                   ApEn_34_36 = [ApEn_34_36, mean(scatola_2)]; 
                else 
                    ApEn_34_36 = [ApEn_34_36, scatola_2];
                end
                if length(scatola_3) > 1
                   ApEn_37_38 = [ApEn_37_38, mean(scatola_3)]; 
                else 
                    ApEn_37_38 = [ApEn_37_38, scatola_3];
                end
             end
           
             ApEn_30_33 = ApEn_30_33'; 
             ApEn_34_36 = ApEn_34_36'; 
             ApEn_37_38 = ApEn_37_38'; 
             ApEn_sani = table(numero_identificativoS, ApEn_30_33, ApEn_34_36, ApEn_37_38);

             ApEn_30_33 = []; 
             ApEn_34_36 = []; 
             ApEn_37_38 = []; 
             for i = 1:length(numero_identificativoM) 
                pat_num = numero_identificativoM(i); 
                scatola_1 = []; 
                scatola_2 = []; 
                scatola_3 = []; 
                for n = 1:size(ParametriMalati, 1) 
                    if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 33 && ParametriMalati.sett_gestazione(n) >= 30
                       scatola_1 = [scatola_1, ParametriMalati.Apen(n)]; 
                    end
                    if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 36 && ParametriMalati.sett_gestazione(n) >= 34
                       scatola_2 = [scatola_2, ParametriMalati.Apen(n)];
                    end
                    if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 38 && ParametriMalati.sett_gestazione(n) >= 37
                       scatola_3 = [scatola_3, ParametriMalati.Apen(n)];
                    end
                end
                    if length(scatola_1) > 1
                       ApEn_30_33 = [ApEn_30_33, mean(scatola_1)]; 
                    else 
                        ApEn_30_33 = [ApEn_30_33, scatola_1];
                    end
                    if length(scatola_2) > 1
                       ApEn_34_36 = [ApEn_34_36, mean(scatola_2)]; 
                    else 
                        ApEn_34_36 = [ApEn_34_36, scatola_2];
                    end
                    if length(scatola_3) > 1
                       ApEn_37_38 = [ApEn_37_38, mean(scatola_3)]; 
                    else 
                        ApEn_37_38 = [ApEn_37_38, scatola_3];
                    end
                 end
                 ApEn_30_33 = ApEn_30_33'; 
                 ApEn_34_36 = ApEn_34_36'; 
                 ApEn_37_38 = ApEn_37_38'; 
                 ApEn_malati = table(numero_identificativoM, ApEn_30_33, ApEn_34_36, ApEn_37_38);
         
         %SampEn        
         SampEn_30_33 = []; 
         SampEn_34_36 = []; 
         SampEn_37_38 = []; 
         for i = 1:length(numero_identificativoS) 
            pat_num = numero_identificativoS(i); 
            scatola_1 = []; 
            scatola_2 = []; 
            scatola_3 = []; 
            for n = 1:size(ParametriSani, 1) 
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 33 && ParametriSani.sett_gestazione(n) >= 30
                   scatola_1 = [scatola_1, ParametriSani.SampEn(n)]; 
                end
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 36 && ParametriSani.sett_gestazione(n) >= 34
                   scatola_2 = [scatola_2, ParametriSani.SampEn(n)];
                end
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 38 && ParametriSani.sett_gestazione(n) >= 37
                   scatola_3 = [scatola_3, ParametriSani.SampEn(n)];
                end
            end
                if length(scatola_1) > 1
                   SampEn_30_33 = [SampEn_30_33, mean(scatola_1)]; 
                else 
                    SampEn_30_33 = [SampEn_30_33, scatola_1];
                end
                if length(scatola_2) > 1
                   SampEn_34_36 = [SampEn_34_36, mean(scatola_2)]; 
                else 
                    SampEn_34_36 = [SampEn_34_36, scatola_2];
                end
                if length(scatola_3) > 1
                   SampEn_37_38 = [SampEn_37_38, mean(scatola_3)]; 
                else 
                    SampEn_37_38 = [SampEn_37_38, scatola_3];
                end
         end 
                 SampEn_30_33 = SampEn_30_33'; 
                 SampEn_34_36 = SampEn_34_36'; 
                 SampEn_37_38 = SampEn_37_38'; 
                 
                 SampEn_sani = table(numero_identificativoS, SampEn_30_33, SampEn_34_36, SampEn_37_38);         
             
         SampEn_30_33 = []; 
         SampEn_34_36 = []; 
         SampEn_37_38 = []; 
         
         for i = 1:length(numero_identificativoM) 
            pat_num = numero_identificativoM(i); 
            scatola_1 = []; 
            scatola_2 = []; 
            scatola_3 = []; 
                for n = 1:size(ParametriMalati, 1) 
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 33 && ParametriMalati.sett_gestazione(n) >= 30
                   scatola_1 = [scatola_1, ParametriMalati.SampEn(n)]; 
                end
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 36 && ParametriMalati.sett_gestazione(n) >= 34
                   scatola_2 = [scatola_2, ParametriMalati.SampEn(n)];
                end
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 38 && ParametriMalati.sett_gestazione(n) >= 37
                   scatola_3 = [scatola_3, ParametriMalati.SampEn(n)];
                end
            end
                if length(scatola_1) > 1
                   SampEn_30_33 = [SampEn_30_33, mean(scatola_1)]; 
                else 
                    SampEn_30_33 = [SampEn_30_33, scatola_1];
                end
                if length(scatola_2) > 1
                   SampEn_34_36 = [SampEn_34_36, mean(scatola_2)]; 
                else 
                    SampEn_34_36 = [SampEn_34_36, scatola_2];
                end
                if length(scatola_3) > 1
                   SampEn_37_38 = [SampEn_37_38, mean(scatola_3)]; 
                else 
                    SampEn_37_38 = [SampEn_37_38, scatola_3];
                end
         end 
                 SampEn_30_33 = SampEn_30_33'; 
                 SampEn_34_36 = SampEn_34_36'; 
                 SampEn_37_38 = SampEn_37_38'; 
                 
                 SampEn_malati = table(numero_identificativoM, SampEn_30_33, SampEn_34_36, SampEn_37_38);         
%LZ2
         LZ2_30_33 = []; 
         LZ2_34_36 = []; 
         LZ2_37_38 = []; 
         for i = 1:length(numero_identificativoS) 
            pat_num = numero_identificativoS(i); 
            scatola_1 = []; 
            scatola_2 = []; 
            scatola_3 = []; 
            for n = 1:size(ParametriSani, 1) 
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 33 && ParametriSani.sett_gestazione(n) >= 30
                   scatola_1 = [scatola_1, ParametriSani.LZ2(n)]; 
                end
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 36 && ParametriSani.sett_gestazione(n) >= 34
                   scatola_2 = [scatola_2, ParametriSani.LZ2(n)];
                end
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 38 && ParametriSani.sett_gestazione(n) >= 37
                   scatola_3 = [scatola_3, ParametriSani.LZ2(n)];
                end
            end
                if length(scatola_1) > 1
                   LZ2_30_33 = [LZ2_30_33, mean(scatola_1)]; 
                else 
                    LZ2_30_33 = [LZ2_30_33, scatola_1];
                end
                if length(scatola_2) > 1
                   LZ2_34_36 = [LZ2_34_36, mean(scatola_2)]; 
                else 
                    LZ2_34_36 = [LZ2_34_36, scatola_2];
                end
                if length(scatola_3) > 1
                   LZ2_37_38 = [LZ2_37_38, mean(scatola_3)]; 
                else 
                    LZ2_37_38 = [LZ2_37_38, scatola_3];
                end
         end 

                 LZ2_30_33 = LZ2_30_33'; 
                 LZ2_34_36 = LZ2_34_36'; 
                 LZ2_37_38 = LZ2_37_38'; 
                 
                 LZ2_sani = table(numero_identificativoS, LZ2_30_33, LZ2_34_36, LZ2_37_38);         
             
         LZ2_30_33 = []; 
         LZ2_34_36 = []; 
         LZ2_37_38 = []; 
         
         for i = 1:length(numero_identificativoM) 
            pat_num = numero_identificativoM(i); 
            scatola_1 = []; 
            scatola_2 = []; 
            scatola_3 = []; 
                for n = 1:size(ParametriMalati, 1) 
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 33 && ParametriMalati.sett_gestazione(n) >= 30
                   scatola_1 = [scatola_1, ParametriMalati.LZ2(n)]; 
                end
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 36 && ParametriMalati.sett_gestazione(n) >= 34
                   scatola_2 = [scatola_2, ParametriMalati.LZ2(n)];
                end
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 38 && ParametriMalati.sett_gestazione(n) >= 37
                   scatola_3 = [scatola_3, ParametriMalati.LZ2(n)];
                end
            end
                if length(scatola_1) > 1
                   LZ2_30_33 = [LZ2_30_33, mean(scatola_1)]; 
                else 
                    LZ2_30_33 = [LZ2_30_33, scatola_1];
                end
                if length(scatola_2) > 1
                   LZ2_34_36 = [LZ2_34_36, mean(scatola_2)]; 
                else 
                    LZ2_34_36 = [LZ2_34_36, scatola_2];
                end
                if length(scatola_3) > 1
                   LZ2_37_38 = [LZ2_37_38, mean(scatola_3)]; 
                else 
                    LZ2_37_38 = [LZ2_37_38, scatola_3];
                end
         end 
                 LZ2_30_33 = LZ2_30_33'; 
                 LZ2_34_36 = LZ2_34_36'; 
                 LZ2_37_38 = LZ2_37_38'; 
                 
                 LZ2_malati = table(numero_identificativoM, LZ2_30_33, LZ2_34_36, LZ2_37_38); 
         %LZ3 
         LZ3_30_33 = []; 
         LZ3_34_36 = []; 
         LZ3_37_38 = []; 
         for i = 1:length(numero_identificativoS) 
            pat_num = numero_identificativoS(i); 
            scatola_1 = []; 
            scatola_2 = []; 
            scatola_3 = []; 
            for n = 1:size(ParametriSani, 1) 
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 33 && ParametriSani.sett_gestazione(n) >= 30
                   scatola_1 = [scatola_1, ParametriSani.LZ3(n)]; 
                end
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 36 && ParametriSani.sett_gestazione(n) >= 34
                   scatola_2 = [scatola_2, ParametriSani.LZ3(n)];
                end
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 38 && ParametriSani.sett_gestazione(n) >= 37
                   scatola_3 = [scatola_3, ParametriSani.LZ3(n)];
                end
            end
                if length(scatola_1) > 1
                   LZ3_30_33 = [LZ3_30_33, mean(scatola_1)]; 
                else 
                    LZ3_30_33 = [LZ3_30_33, scatola_1];
                end
                if length(scatola_2) > 1
                   LZ3_34_36 = [LZ3_34_36, mean(scatola_2)]; 
                else 
                    LZ3_34_36 = [LZ3_34_36, scatola_2];
                end
                if length(scatola_3) > 1
                   LZ3_37_38 = [LZ3_37_38, mean(scatola_3)]; 
                else 
                    LZ3_37_38 = [LZ3_37_38, scatola_3];
                end
         end 

                 LZ3_30_33 = LZ3_30_33'; 
                 LZ3_34_36 = LZ3_34_36'; 
                 LZ3_37_38 = LZ3_37_38'; 
                 
                 LZ3_sani = table(numero_identificativoS, LZ3_30_33, LZ3_34_36, LZ3_37_38);         
             
         LZ3_30_33 = []; 
         LZ3_34_36 = []; 
         LZ3_37_38 = []; 
         
         for i = 1:length(numero_identificativoM) 
            pat_num = numero_identificativoM(i); 
            scatola_1 = []; 
            scatola_2 = []; 
            scatola_3 = []; 
                for n = 1:size(ParametriMalati, 1) 
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 33 && ParametriMalati.sett_gestazione(n) >= 30
                   scatola_1 = [scatola_1, ParametriMalati.LZ3(n)]; 
                end
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 36 && ParametriMalati.sett_gestazione(n) >= 34
                   scatola_2 = [scatola_2, ParametriMalati.LZ3(n)];
                end
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 38 && ParametriMalati.sett_gestazione(n) >= 37
                   scatola_3 = [scatola_3, ParametriMalati.LZ3(n)];
                end
            end
                if length(scatola_1) > 1
                   LZ3_30_33 = [LZ3_30_33, mean(scatola_1)]; 
                else 
                    LZ3_30_33 = [LZ3_30_33, scatola_1];
                end
                if length(scatola_2) > 1
                   LZ3_34_36 = [LZ3_34_36, mean(scatola_2)]; 
                else 
                    LZ3_34_36 = [LZ3_34_36, scatola_2];
                end
                if length(scatola_3) > 1
                   LZ3_37_38 = [LZ3_37_38, mean(scatola_3)]; 
                else 
                    LZ3_37_38 = [LZ3_37_38, scatola_3];
                end
         end 
                 LZ3_30_33 = LZ3_30_33'; 
                 LZ3_34_36 = LZ3_34_36'; 
                 LZ3_37_38 = LZ3_37_38'; 
                 
                 LZ3_malati = table(numero_identificativoM, LZ3_30_33, LZ3_34_36, LZ3_37_38); 
      %LF  
         LF_30_33 = []; 
         LF_34_36 = []; 
         LF_37_38 = []; 
         for i = 1:length(numero_identificativoS) 
            pat_num = numero_identificativoS(i); 
            scatola_1 = []; 
            scatola_2 = []; 
            scatola_3 = []; 
            for n = 1:size(ParametriSani, 1) 
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 33 && ParametriSani.sett_gestazione(n) >= 30
                   scatola_1 = [scatola_1, ParametriSani.LF(n)]; 
                end
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 36 && ParametriSani.sett_gestazione(n) >= 34
                   scatola_2 = [scatola_2, ParametriSani.LF(n)];
                end
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 38 && ParametriSani.sett_gestazione(n) >= 37
                   scatola_3 = [scatola_3, ParametriSani.LF(n)];
                end
            end
                if length(scatola_1) > 1
                   LF_30_33 = [LF_30_33, mean(scatola_1)]; 
                else 
                    LF_30_33 = [LF_30_33, scatola_1];
                end
                if length(scatola_2) > 1
                   LF_34_36 = [LF_34_36, mean(scatola_2)]; 
                else 
                    LF_34_36 = [LF_34_36, scatola_2];
                end
                if length(scatola_3) > 1
                   LF_37_38 = [LF_37_38, mean(scatola_3)]; 
                else 
                    LF_37_38 = [LF_37_38, scatola_3];
                end
         end 

                 LF_30_33 = LF_30_33'; 
                 LF_34_36 = LF_34_36'; 
                 LF_37_38 = LF_37_38'; 
                 
                 LF_sani = table(numero_identificativoS, LF_30_33, LF_34_36, LF_37_38);         
             
         LF_30_33 = []; 
         LF_34_36 = []; 
         LF_37_38 = []; 
         
         for i = 1:length(numero_identificativoM) 
            pat_num = numero_identificativoM(i); 
            scatola_1 = []; 
            scatola_2 = []; 
            scatola_3 = []; 
                for n = 1:size(ParametriMalati, 1) 
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 33 && ParametriMalati.sett_gestazione(n) >= 30
                   scatola_1 = [scatola_1, ParametriMalati.LF(n)]; 
                end
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 36 && ParametriMalati.sett_gestazione(n) >= 34
                   scatola_2 = [scatola_2, ParametriMalati.LF(n)];
                end
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 38 && ParametriMalati.sett_gestazione(n) >= 37
                   scatola_3 = [scatola_3, ParametriMalati.LF(n)];
                end
            end
                if length(scatola_1) > 1
                   LF_30_33 = [LF_30_33, mean(scatola_1)]; 
                else 
                    LF_30_33 = [LF_30_33, scatola_1];
                end
                if length(scatola_2) > 1
                   LF_34_36 = [LF_34_36, mean(scatola_2)]; 
                else 
                    LF_34_36 = [LF_34_36, scatola_2];
                end
                if length(scatola_3) > 1
                   LF_37_38 = [LF_37_38, mean(scatola_3)]; 
                else 
                    LF_37_38 = [LF_37_38, scatola_3];
                end
         end 
                 LF_30_33 = LF_30_33'; 
                 LF_34_36 = LF_34_36'; 
                 LF_37_38 = LF_37_38'; 
                 
                 LF_malati = table(numero_identificativoM, LF_30_33, LF_34_36, LF_37_38); 
 %MF  
         MF_30_33 = []; 
         MF_34_36 = []; 
         MF_37_38 = []; 
         for i = 1:length(numero_identificativoS) 
            pat_num = numero_identificativoS(i); 
            scatola_1 = []; 
            scatola_2 = []; 
            scatola_3 = []; 
            for n = 1:size(ParametriSani, 1) 
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 33 && ParametriSani.sett_gestazione(n) >= 30
                   scatola_1 = [scatola_1, ParametriSani.MF(n)]; 
                end
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 36 && ParametriSani.sett_gestazione(n) >= 34
                   scatola_2 = [scatola_2, ParametriSani.MF(n)];
                end
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 38 && ParametriSani.sett_gestazione(n) >= 37
                   scatola_3 = [scatola_3, ParametriSani.MF(n)];
                end
            end
                if length(scatola_1) > 1
                   MF_30_33 = [MF_30_33, mean(scatola_1)]; 
                else 
                    MF_30_33 = [MF_30_33, scatola_1];
                end
                if length(scatola_2) > 1
                   MF_34_36 = [MF_34_36, mean(scatola_2)]; 
                else 
                    MF_34_36 = [MF_34_36, scatola_2];
                end
                if length(scatola_3) > 1
                   MF_37_38 = [MF_37_38, mean(scatola_3)]; 
                else 
                    MF_37_38 = [MF_37_38, scatola_3];
                end
         end 

                 MF_30_33 = MF_30_33'; 
                 MF_34_36 = MF_34_36'; 
                 MF_37_38 = MF_37_38'; 
                 
                 MF_sani = table(numero_identificativoS, MF_30_33, MF_34_36, MF_37_38);         
             
         MF_30_33 = []; 
         MF_34_36 = []; 
         MF_37_38 = []; 
         
         for i = 1:length(numero_identificativoM) 
            pat_num = numero_identificativoM(i); 
            scatola_1 = []; 
            scatola_2 = []; 
            scatola_3 = []; 
                for n = 1:size(ParametriMalati, 1) 
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 33 && ParametriMalati.sett_gestazione(n) >= 30
                   scatola_1 = [scatola_1, ParametriMalati.MF(n)]; 
                end
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 36 && ParametriMalati.sett_gestazione(n) >= 34
                   scatola_2 = [scatola_2, ParametriMalati.MF(n)];
                end
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 38 && ParametriMalati.sett_gestazione(n) >= 37
                   scatola_3 = [scatola_3, ParametriMalati.MF(n)];
                end
            end
                if length(scatola_1) > 1
                   MF_30_33 = [MF_30_33, mean(scatola_1)]; 
                else 
                    MF_30_33 = [MF_30_33, scatola_1];
                end
                if length(scatola_2) > 1
                   MF_34_36 = [MF_34_36, mean(scatola_2)]; 
                else 
                   MF_34_36 = [MF_34_36, scatola_2];
                end
                if length(scatola_3) > 1
                   MF_37_38 = [MF_37_38, mean(scatola_3)]; 
                else 
                    MF_37_38 = [MF_37_38, scatola_3];
                end
         end 
                 MF_30_33 = MF_30_33'; 
                 MF_34_36 = MF_34_36'; 
                 MF_37_38 = MF_37_38'; 
                 
                 MF_malati = table(numero_identificativoM, MF_30_33, MF_34_36, MF_37_38); 
                 
                   %HF  
         HF_30_33 = []; 
         HF_34_36 = []; 
         HF_37_38 = []; 
         for i = 1:length(numero_identificativoS) 
            pat_num = numero_identificativoS(i); 
            scatola_1 = []; 
            scatola_2 = []; 
            scatola_3 = []; 
            for n = 1:size(ParametriSani, 1) 
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 33 && ParametriSani.sett_gestazione(n) >= 30
                   scatola_1 = [scatola_1, ParametriSani.HF(n)]; 
                end
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 36 && ParametriSani.sett_gestazione(n) >= 34
                   scatola_2 = [scatola_2, ParametriSani.HF(n)];
                end
                if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 38 && ParametriSani.sett_gestazione(n) >= 37
                   scatola_3 = [scatola_3, ParametriSani.HF(n)];
                end
            end
                if length(scatola_1) > 1
                   HF_30_33 = [HF_30_33, mean(scatola_1)]; 
                else 
                    HF_30_33 = [HF_30_33, scatola_1];
                end
                if length(scatola_2) > 1
                   HF_34_36 = [HF_34_36, mean(scatola_2)]; 
                else 
                    HF_34_36 = [HF_34_36, scatola_2];
                end
                if length(scatola_3) > 1
                   HF_37_38 = [HF_37_38, mean(scatola_3)]; 
                else 
                    HF_37_38 = [HF_37_38, scatola_3];
                end
         end 

                 HF_30_33 = HF_30_33'; 
                 HF_34_36 = HF_34_36'; 
                 HF_37_38 = HF_37_38'; 
                 
                 HF_sani = table(numero_identificativoS, HF_30_33, HF_34_36, HF_37_38);         
             
         HF_30_33 = []; 
         HF_34_36 = []; 
         HF_37_38 = []; 
         
         for i = 1:length(numero_identificativoM) 
            pat_num = numero_identificativoM(i); 
            scatola_1 = []; 
            scatola_2 = []; 
            scatola_3 = []; 
                for n = 1:size(ParametriMalati, 1) 
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 33 && ParametriMalati.sett_gestazione(n) >= 30
                   scatola_1 = [scatola_1, ParametriMalati.HF(n)]; 
                end
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 36 && ParametriMalati.sett_gestazione(n) >= 34
                   scatola_2 = [scatola_2, ParametriMalati.HF(n)];
                end
                if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 38 && ParametriMalati.sett_gestazione(n) >= 37
                   scatola_3 = [scatola_3, ParametriMalati.HF(n)];
                end
            end
                if length(scatola_1) > 1
                   HF_30_33 = [HF_30_33, mean(scatola_1)]; 
                else 
                    HF_30_33 = [HF_30_33, scatola_1];
                end
                if length(scatola_2) > 1
                   HF_34_36 = [HF_34_36, mean(scatola_2)]; 
                else 
                    HF_34_36 = [HF_34_36, scatola_2];
                end
                if length(scatola_3) > 1
                   HF_37_38 = [HF_37_38, mean(scatola_3)]; 
                else 
                    HF_37_38 = [HF_37_38, scatola_3];
                end
         end 
                 HF_30_33 = HF_30_33'; 
                 HF_34_36 = HF_34_36'; 
                 HF_37_38 = HF_37_38'; 
                 
                 HF_malati = table(numero_identificativoM, HF_30_33, HF_34_36, HF_37_38); 
              
