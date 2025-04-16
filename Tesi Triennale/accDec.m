function [grandiAcc,piccoleAcc,grandiDec,piccoleDec] = accDec(intacc,intdec)

%%  aggiustamento delle accelerazioni e decelerazioni secondo la nuova definizione

try 
    grandiAcc=intacc(intacc(:,3)>15,:);
    piccoleAcc=intacc(intacc(:,3)<=15,:);
    grandiDec=intdec(intdec(:,3)>20,:);
    piccoleDec=intdec(intdec(:,3)<=20,:);
catch
    1
end


% dec=[grandiDec' piccoleDec']';

end