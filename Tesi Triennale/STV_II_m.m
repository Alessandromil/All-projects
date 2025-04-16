function [valorestv, valorestv_dalton, valoreii, valoreii_2ctg]=STV_II_m(minuto,ESCLUSO,SOGLIA)
% NLG, Dipartimento di Bioingegneria, 1999.  
% Revisione 5; 15 Maggio 2000, Politecnico di Milano, Italia.
% Codice Matlab scritto da Roberto Sassi

if nargin > 1,
   %% I PUNTI POSTI AL VALORE "ESCLUSO" VENGONO TRASCURATI
   intnonNULLI=cercaintervalli(minuto~=ESCLUSO);
   NUMintnonNULLI=size(intnonNULLI);
   differenza=[];
   spezzone_nonNULLI_max=0;
   if NUMintnonNULLI~=0,
      for ii=1:NUMintnonNULLI,
         inizio=intnonNULLI(ii,1);
         fine=intnonNULLI(ii,2);
         if fine > inizio,
            differenza=[differenza; abs(diff(minuto(inizio:fine,1)))];
            if spezzone_nonNULLI_max < fine-inizio+1,
               spezzone_nonNULLI_max=fine-inizio+1;
            end,
         end,
      end,
   end,
   if spezzone_nonNULLI_max < SOGLIA,
      valorestv=NaN;
      valorestv_dalton=NaN;
      valoreii=NaN;
      valoreii_2ctg=NaN;
   else,
      valorestv=mean(differenza);
      %% VIENE CORRETTO PER EVITARE CHE SIA DISTORTO DALLA EVENTUALE MANCANZA DI PUNTI
      numeropunti=size(differenza,1);
      valorestv_dalton=(sum(differenza)/2)*((size(minuto,1)-1)/numeropunti);
      valoreii_2ctg=std(differenza);
      valoreii=valoreii_2ctg/valorestv;
   end,
else
   %% TUTTI I PUNTI ENTRANO NEL CALCOLO
   differenza=abs(diff(minuto(:,1)));
   valorestv=mean(differenza);
   valorestv_dalton=sum(differenza)/2;
   valoreii_2ctg=std(differenza);
   valoreii=valoreii_2ctg/valorestv;
end,

return,
