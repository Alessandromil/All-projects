cont_1 = []; 
cont_2 = []; 
cont_3 = [];
cont_paz = []; 
for i = 1:length(numero_identificativoM) 
        pat_num = numero_identificativoM(i);  
        for n = 1:size(ParametriMalati, 1) 
            if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 33 && ParametriMalati.sett_gestazione(n) >= 30
              cont_1 = [cont_1, pat_num]; 
            end
            if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 36 && ParametriMalati.sett_gestazione(n) >= 34
              cont_2 = [cont_2, pat_num]; 
            end
            if ParametriMalati.patnum(n) == pat_num && ParametriMalati.sett_gestazione(n) <= 38 && ParametriMalati.sett_gestazione(n) >= 37
               cont_3 = [cont_3, pat_num];
            end
        end
end 
cont_paz = unique([cont_1, cont_2, cont_3]); 
