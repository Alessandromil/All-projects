function [dati24bpmsenza,dati24mssenza]=excludeAccDec(dati24bpm,dati24ms,grandiAcc,piccoleAcc,grandiDec,piccoleDec,escludi,ESCLUSO,INTERPOLATO)

if nargin ==6
    
    ESCLUDI_ACC_GRANDI = 1;
    ESCLUDI_ACC_PICCOLE = 0;
    ESCLUDI_DEC_GRANDI= 1;
    ESCLUDI_DEC_PICCOLE= 1;
    ESCLUDI_DATI_ROSSI = 1;
    ESCLUSO = -1;
    INTERPOLATO = 96;
    
else
    ESCLUDI_ACC_GRANDI=escludi(1);
    ESCLUDI_ACC_PICCOLE=escludi(2);
    ESCLUDI_DEC_GRANDI=escludi(3);
    ESCLUDI_DEC_PICCOLE=escludi(4);
    ESCLUDI_DATI_ROSSI =escludi(5);
    
end

dati24bpmsenza=dati24bpm(:,1);
dati24mssenza=dati24ms;
qualita=dati24bpm(:,2);

if ESCLUDI_ACC_GRANDI
    for jj=1:size(grandiAcc,1)
        dati24bpmsenza(grandiAcc(jj,1):grandiAcc(jj,2))=ESCLUSO;
        dati24mssenza(grandiAcc(jj,1):grandiAcc(jj,2))=ESCLUSO;
    end
end

if ESCLUDI_ACC_PICCOLE
    for jj=1:size(piccoleAcc,1)
        dati24bpmsenza(piccoleAcc(jj,1):piccoleAcc(jj,2))=ESCLUSO;
        dati24mssenza(piccoleAcc(jj,1):piccoleAcc(jj,2))=ESCLUSO;
    end
end

if ESCLUDI_DEC_GRANDI
    for jj=1:size(grandiDec,1)
        dati24bpmsenza(grandiDec(jj,1):grandiDec(jj,2))=ESCLUSO;
        dati24mssenza(grandiDec(jj,1):grandiDec(jj,2))=ESCLUSO;
    end
end

if ESCLUDI_DEC_PICCOLE
    for jj=1:size(piccoleDec,1)
        dati24bpmsenza(piccoleDec(jj,1):piccoleDec(jj,2))=ESCLUSO;
        dati24mssenza(piccoleDec(jj,1):piccoleDec(jj,2))=ESCLUSO;
    end
end

if ESCLUDI_DATI_ROSSI
    pos_dati_rossi=find(qualita==INTERPOLATO);
    if ~isempty(pos_dati_rossi)
        dati24bpmsenza(pos_dati_rossi,1)=ESCLUSO;
        dati24mssenza(pos_dati_rossi,1)=ESCLUSO;
    end
end

end
