settimane = 30:38; 
cont_s = zeros(size(settimane)); 
for j = 1:length(settimane)
    settimana = settimane(j); 
    corrente = 0; 
        for n = 1:size(Data(:, 1), 1)
            if Data.sett_gestazione(n) == settimana
                corrente = corrente + 1; 
            end
        end
    cont_s(j) = corrente; 
end

cont_m = zeros(size(settimane)); 
for j = 1:length(settimane)
    settimana = settimane(j); 
    corrente = 0; 
        for n = 1:size(Data_m(:, 1), 1)
            if Data_m.sett_gestazione(n) == settimana
                corrente = corrente + 1; 
            end
        end
    cont_m(j) = corrente; 
end