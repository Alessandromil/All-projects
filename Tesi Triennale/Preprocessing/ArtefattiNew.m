% function [esame_out, num_correzioni]=ArtefattiNew(parametri, esame_in)
%
% Localizza eventali artefatti sul segnale cardiografico.
% Prende spunto dall'algoritmo descritto in:
%
% van Geijn H P, Jongsma H W, de Hann J and Eskes T K A B
% "Analysis of heart rate and beat-to-beat variability: 
% Interval difference index"
% Am J Obstet Gynecol, 138 (1980) 246-252.
%
% Parametri è un vettore riga di 64 elementi; solo alcuni sono rilevanti 
% per questa routine:
% PARAMETRI(8)=ROSSO;
% PARAMETRI(41)=FCF_MAX;
% PARAMETRI(42)=LUNGH_MAX_GAP_INAFFIDABILI;
%
% SEE ALSO: MantelNew.

% NLG, Dipartimento di Bioingegneria, Politecnico di Milano, Italia.
% versione 6; 23 Giugno 2000
% Codice C/Matlab scritto da Roberto Sassi
