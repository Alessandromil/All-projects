In linea di massima le funzioni si chiamano su finestre di 1 o 3 minuti e poi si fa la media. 
Alcuni file sono in c, vi ho passato il mex tramite cui da Windows dovreste poter chiamare le funzioni come se fossero delle normali 
funzioni di MATLAB. Escludete dal calcolo le finesre in cui ci sono interpolazioni 
(guardate il vettore qualita, sono interpolati se la qualità è maggiore di 64) 


STV_II calcola Short term variability e Interval Index. 
Passare alla funzione finestre di 1 minuto del segnale campionato 
a 24 punti per minuto escluse grandi accelerazioni e decelerazioni (FHR24mssenza) (potete ignorare gli altri input)

apsampen è precompliato. la funzione si chiama così:
[Apen, SampEn]=apsampen(spezzzone,m,r,1); 
m e r si possono scegliere (2 e 0.15 sono valori sensati, per esempio). 
1 significa che r è già normalizzato rispetto alla deviazione standard (implicitamente lo moltiplica per std(spez)). 
Spezzone è in 120 punti per minuto. è indifferente se in ms o bpm
ho creato la funzione ComputeLZ che chiama delle funzioni precomplicate. 
Comunque basta passare lo spezzone del segnale (in 120 punti per minuto. è indifferente se in ms o in bpm) e vi ritorna LZ binario e ternario.
