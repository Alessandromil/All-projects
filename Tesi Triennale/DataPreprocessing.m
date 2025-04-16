clear
close
clc

%caricate i dati in matlab (qualcosa tipo:
%SANI = load(SANI.mat);
%Data = SANI.SANI;) %nel caso lo importasse come struct non sono sicuro

NamesPar 
parametri=cell2mat(struct2cell(par));

%% pre processing del segnale
%restituisce:
% linea di base (24 punti per minuto)
% accelerazioni (indice inizio,indice fine, altezza picco, inidice picco - 24 punti per minuto)
% decelerazioni (come sopra)
% FHR con interpolazioni in bpm (120 punti per minuto - come l'originale)
% FHR con interpolazioni in bpm (24 punti per minuto  - sottocampionato)
% FHR con interpolazioni in ms (120 punti per minuto - come l'originale)
% FHR con interpolazioni in ms (24 punti per minuto  - sottocampionato)
% linea di base (120 punti per minuto)
% accelerazioni (indice inizio,indice fine, altezza picco, inidice picco - 120 punti per minuto)
% decelerazioni (come sopra)
[Data.baseFHR24bpm, Data.intacc24bpm, Data.intdec24bpm, Data.FHR120bpm, Data.FHR24bpm,Data.FHR120ms,Data.FHR24ms,Data.base120bpm,Data.intacc120bpm,Data.intdec120bpm] = cellfun(@(x,y) PreProc([x',y'],parametri,1),Data.FHR,Data.QUALITA,'uni',0);
%tolgo le osservazion in cui è fallito il pre processing
idx = cellfun(@(x) isempty(x), Data.baseFHR24bpm);
Data(idx,:)=[];
[Data.grandiAcc,Data.piccoleAcc,Data.grandiDec,Data.piccoleDec] = cellfun(@(x,y) accDec(x,y), Data.intacc24bpm, Data.intdec24bpm,'uni',0);
[Data.FHR24bpmsenza,Data.FHR24mssenza]=cellfun(@(a,b,c,d,e,f) excludeAccDec(a,b,c,d,e,f),Data.FHR24bpm,Data.FHR24ms,Data.grandiAcc,Data.piccoleAcc,Data.grandiDec,Data.piccoleDec,'uni',0);



%save(cartella,'Data','-v7.3');



