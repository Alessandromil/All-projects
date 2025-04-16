%esempio per chiamare le funzioni

%genero segnale simulato per prova (ignorate agilmente la procedura)
segnale120 = randn(10000,1);
segnale24 = randn(2000,1); %per STV
qualita = ones(10000,1)*32; qualita(abs(randn(10000,1))>1) = 64; qualita(abs(randn(10000,1))>2) = 128;

%scorro sul segnale e calcolo i parametri
N120 = 120*3;
N3min = floor(length(segnale120)/N120); %numero finestre di 3 minuti
N1min = floor(length(segnale120)/120);

%calcolo parametri da 3 minuti
for i = 1:N3min
    spezzone = segnale120(1+(i-1)*N120:i*N120);
    spezzoneQualita = qualita(1+(i-1)*N120:i*N120);   
    %controllo qualità
    if sum(spezzoneQualita>64)<(0.05*N120) %meno del 5% di interpolati
        [Apen(i),SampEn(i)] = apsampen(spezzone,2,0.15,1);
        [LZ2(i),LZ3(i)] = ComputeLZ(spezzone);
    else
        [Apen(i),SampEn(i),LZ2(i),LZ3(i)] = deal(nan);
    end
end

%è solo STV, che è calcolata in 24 punti per minuto
for i = 1:N1min
    spezzone24 = segnale24(1+(i-1)*24:i*24);
    spezzoneQualita = qualita(1+(i-1)*120:i*120); %occhio che qualità in120
    %controllo qualità
    if sum(spezzoneQualita>64)<(0.05*120) %meno del 5% di interpolati
        [STV(i), ~, ~, ~]=STV_II_m(spezzone24);
    else
        STV(i) = nan;
    end

end

%rimuovo outlier e calcolo media
STV(abs(STV)>4*std(STV)) = nan; 
STV = nanmean(STV);

%etcetc