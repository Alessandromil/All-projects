function CreaBoxplot(Y)

labels = {'30-33', '30-33', '34-36', '34-36', '37-38', '37-38'};

figure(1)
h = boxplot(Y, 'labels', labels, 'Colors', 'k', 'BoxStyle','outline');
blue = [0.9290 0.6940 0.1250]; 
red = [0.4660 0.6740 0.1880]; 
colors = [blue; red; blue; red; blue; red];
m = findobj(gca,'Tag','Box');
for j=1:length(m)
    patch(get(m(j),'XData'),get(m(j),'YData'),colors(j, :),'FaceAlpha',.5);
end
set(h,{'linew'},{1})
% legend('Malati', 'Sani')
xlabel('Settimane')
%ylabel('Media [bpm]')
box on 

end