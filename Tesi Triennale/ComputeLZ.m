function [LZ2,LZ3] = ComputeLZ(spezzone)
    %normalizzo
    ss=spezzone./std(spezzone);
    snorm=ss-mean(ss)*ones(size(ss));
    %codifico
    ssn = codeB(snorm, 2, 0.02);
    ssn1 = codeB(snorm, 3, 0.01);
    % LZ spezzone normalizzato
    LZ2 = LZ(string(ssn),2);   %LZ k=2 simboli con soglia p = 0.02
    LZ3 = LZ(string(ssn1),3);   %LZ k=3 simboli con soglia p = 0.01
end