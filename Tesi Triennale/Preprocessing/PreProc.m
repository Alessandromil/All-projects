function  [base, intacc, intdec, dati120bpm, dati24bpm,dati120ms,dati24ms,varargout] = PreProc(esame,PARAMETRI,varargin)

if size(esame,1)<size(esame,2)
    esame = esame';
end

[esame_corretto, ~]=ArtefattiNew(PARAMETRI, esame);

if nargin ==2 %versione base
    [base, intacc, intdec, dati120bpm, dati24bpm] = ...
        MantelNew(PARAMETRI, esame_corretto);

%aggiungi interpolazione accelerazioni e linea di base
elseif nargin ==3 %linea di base e interpolazioni 120
    try
        [base, intacc, intdec, dati120bpm, dati24bpm, baseI, intaccI, intdecI,~] = ...
            MantelNew(PARAMETRI, esame_corretto);
    
        varargout{1} = baseI;
        varargout{2} = intaccI;
        varargout{3} = intdecI;
    catch
        warning('Mantel Non eseguito - capisci perché')
        [base, intacc, intdec, dati120bpm, dati24bpm,dati120ms,dati24ms,varargout{1},varargout{2},varargout{3}] = deal([]);
        return
    end
end

dati120ms = 60000./dati120bpm(:,1);
dati24ms  = 60000./dati24bpm(:,1);

end
