% function [base, intacc, intdec, dati, datir,
%           base_interp, intacc_interp, intdec_interp,
%           buco_iniziale]=mantelnew(parametri, esame)
%
% Calcola linea di base, accelerazioni e decelerazioni, attraverso
% l'algoritmo di [Mantel et al.].
% Ref: International Journal of Bio-Medical Computing
% Volume 25 No. 4 May 1990
% "Computer analysis of antepartum fetal heart rate:
%  1. baseline determination. (pp. 261-272)
%  2. detection of accelerations and decelerations. (pp. 273-286)"
% 
% Parametri è un vettore riga di 64 elementi; solo alcuni sono rilevanti 
% per questa routine:
% PARAMETRI(6)=VERDE;
% PARAMETRI(7)=GIALLO;
% PARAMETRI(8)=ROSSO;
% PARAMETRI(9)=INTERPOLATO;
% PARAMETRI(21)=SOGLIA_INF_ACC;
% PARAMETRI(22)=SOGLIA_SUP_ACC;
% PARAMETRI(23)=LUNGH_MIN_ACC;
% PARAMETRI(24)=SOGLIA_INF_DEC;
% PARAMETRI(25)=SOGLIA_SUP_DEC;
% PARAMETRI(26)=LUNGH_MIN_DEC;
% PARAMETRI(27)=LUNGH_MAX_GAP;
% PARAMETRI(28)=PROB_MIN;
% PARAMETRI(29)=FLAG_FASE_2;
%
% In uscita:
% base           -> baseline (24 punti/minuto);
% dati           -> esame(:,1:2), corretto (120 punti/minuto);
% datiR          -> esame(:,1:2), corretto, filtrato e sottocampionato
%                   (24 punti/minuto);
% intacc, intdec -> accelerazioni e decelerazioni nel seguente formato:
%
% [(colonna 1)  (colonna2)  (colonna3)  (colonna4)]
%   inizio       fine        valoremax   posizionemax
%
% Ulteriori uscite opzionali sono:
%
% base_interp    -> baseline (120 punti/minuto);
% intacc_interp,
% intdec_interp  -> accelerazioni e decelerazioni (120 punti/minuto).
% buco_iniziale ->  numero di campioni ROSSI con cui inizia la
%                   variabile esame (120 punti/minuto).
%
% TUTTI i valori sono in battiti per minuto (bpm)
%
% ATTENZIONE: base, intacc e intdec si riferiscono al segnale sotto-
% campionato, come suggerito da Mantel (24 punti/minuto). Utilizzare
% base_interp, intacc_interp e intdec_interp per valori riferiti al
% segnale iniziale (120 punti/minuto).
%
% SEE ALSO: mantel, mantelC, mantelDLL, mantelFULL.

% NLG, Dipartimento di Bioingegneria, Politecnico di Milano, Italia.
% Revisione 6; 23 Giugno 2000
% Codice C/Matlab scritto da Roberto Sassi
% Copyright (c) 1998-99 2000 by Roberto Sassi
