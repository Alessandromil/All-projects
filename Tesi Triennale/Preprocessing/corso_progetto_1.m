campioni_settimana_1 = []; 
campioni_settimana_2 = []; 
campioni_settimana_3 = []; 
counter = 0; 
campioni = []; 

for i=1:size(Data, 1)
    if Data.sett_gestazione(i) <= 33 && Data.sett_gestazione(i) >= 30 
        campioni_settimana_1 = [campioni_settimana_1, Data.patnum(i)];
    end
    campioni_settimana_1 = unique(campioni_settimana_1);
end


for i=1:size(Data, 1)
    if Data.sett_gestazione(i) <= 36 && Data.sett_gestazione(i) >= 34 
        campioni_settimana_2 = [campioni_settimana_2, Data.patnum(i)];
    end
    campioni_settimana_2 = unique(campioni_settimana_2);
end

for i=1:size(Data, 1)
    if Data.sett_gestazione(i) <= 38 && Data.sett_gestazione(i) >= 37 
        campioni_settimana_3 = [campioni_settimana_3, Data.patnum(i)];
    end
    campioni_settimana_3 = unique(campioni_settimana_3);
end

campioni = [campioni_settimana_1, campioni_settimana_2, campioni_settimana_3];
campioni_finali = [];

for i = 1:size(Data, 1)
for n = 1: size(campioni)
    if Data.patnum(i) == campioni(n)
        counter = counter + 1; 
    end
end

if counter >= 3
    campioni_finali = [campioni_finali, Data.patnum(i)];
end 
counter = 0; 
end 
