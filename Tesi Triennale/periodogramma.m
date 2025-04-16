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
            
%             figure
%             plot(freq/pi,pow2db(psdx))
%             grid on
%             title("Periodogram Using FFT")
%             xlabel("Normalized Frequency (\times\pi rad/sample)")
%             ylabel("Power/Frequency (dB/(rad/sample))")
            
            PTot = var(spezzone);
            
            LF(n) = (0.5*sum(psdx(freq>0.03 & freq<0.15)))/PTot;
            MF(n) = (0.5*sum(psdx(freq>0.15 & freq<0.5)))/PTot;
            HF(n) = (0.5*sum(psdx(freq>0.5 & freq<1)))/PTot;

        else
            [Apen(n),SampEn(n),LZ2(n),LZ3(n),LF(n),MF(n),HF(n)] = deal(nan);
        end
    end
end 
            figure
            plot(freq/pi,pow2db(psdx))
            grid on
            title("Periodogram Using FFT")
            xlabel("Normalized Frequency (\times\pi rad/sample)")
            ylabel("Power/Frequency (dB/(rad/sample))")
      