cont_1 = []; 
cont_2 = []; 
cont_3 = [];
cont_paz = []; 
for i = 1:length(numero_identificativoS) 
        pat_num = numero_identificativoS(i);  
        for n = 1:size(ParametriSani, 1) 
            if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 33 && ParametriSani.sett_gestazione(n) >= 30
              cont_1 = [cont_1, pat_num]; 
            end
            if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 36 && ParametriSani.sett_gestazione(n) >= 34
              cont_2 = [cont_2, pat_num]; 
            end
            if ParametriSani.patnum(n) == pat_num && ParametriSani.sett_gestazione(n) <= 38 && ParametriSani.sett_gestazione(n) >= 37
               cont_3 = [cont_3, pat_num];
            end
        end
end 
cont_paz = unique([cont_1, cont_2, cont_3]); 
