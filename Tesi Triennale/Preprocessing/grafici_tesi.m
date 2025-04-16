clc 
% clear
close all 

Y = [media_sani.media_30_33, media_malati.media_30_33, media_sani.media_34_36, media_malati.media_34_36, media_sani.media_37_38, media_malati.media_37_38]; 
labels = {'30-33', '30-33', '34-36', '34-36', '37-38', '37-38'};
figure(1)
h = boxplot(Y, 'labels', labels, 'Colors', 'k', 'BoxStyle','outline');
viola = [0.4940 0.1840 0.5560]; 
fucsia = [0.6350 0.0780 0.1840]; 
colors = [viola; fucsia; viola; fucsia; viola; fucsia];
m = findobj(gca,'Tag','Box');
for j=1:length(m)
    patch(get(m(j),'XData'),get(m(j),'YData'),colors(j, :),'FaceAlpha',.5);
end
set(h,{'linew'},{1})
%title('Media della linea di base'); 
legend('Malati', 'Sani')
box on 

figure(2)
Y_1 = [LF_sani.LF_30_33, LF_malati.LF_30_33, LF_sani.LF_34_36, LF_malati.LF_34_36, LF_sani.LF_37_38, ...
    LF_malati.LF_37_38];
h = boxplot(Y_1, 'labels', labels, 'Colors', 'k', 'BoxStyle','outline');
labels = {'30-33', '30-33', '34-36', '34-36', '37-38', '37-38'};
xlabel('Settimana di gestazione'); 
ylabel('Densità di potenza spettrale sulle LF%')
verde = [0.4660 0.6740 0.1880]; 
giallo = [0.9290 0.6940 0.1250]; 
colors = [giallo; verde; giallo; verde; giallo; verde];
m = findobj(gca,'Tag','Box');
for j=1:length(m)
    patch(get(m(j),'XData'),get(m(j),'YData'),colors(j, :),'FaceAlpha',.5);
end
set(h,{'linew'},{1})
%title('Andamento LF nelle settimane'); 
legend('Malati', 'Sani')
box on 

Y = [STV_sani.STV_30_33, STV_malati.STV_30_33]; 
labels = {'30-33', '30-33'};

figure(3)
h = boxplot(Y, 'labels', labels, 'Colors', 'k', 'BoxStyle','outline');
viola = [0.4940 0.1840 0.5560]; 
fucsia = [0.6350 0.0780 0.1840]; 
colors = [viola; fucsia];
m = findobj(gca,'Tag','Box');
for j=1:length(m)
    patch(get(m(j),'XData'),get(m(j),'YData'),colors(j, :),'FaceAlpha',.5);
end
set(h,{'linew'},{1})
title('STV della FHR'); 
legend('Malati', 'Sani')
box on 



