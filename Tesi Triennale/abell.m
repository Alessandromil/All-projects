clc 
clear
close all

a = -pi : pi/2 : pi;                                % Define Corners
ph = pi/4;                                          % Define Angular Orientation (‘Phase’)
x = [cos(a+ph); cos(a+ph)]/cos(ph);
y = [sin(a+ph); sin(a+ph)]/sin(ph);
z = [-ones(size(a)); ones(size(a))];
figure
surf(x, y, z, 'FaceColor',[0.8500 0.3250 0.0980])                      % Plot Cube
hold on
patch(x', y', z', [0.9290 0.6940 0.1250])                              % Make Cube Appear Solid
hold off
axis([ -1  1    -1  1    -1  1]*1.5)
xlabel('Chi')
xticks(-1:0.5:1.5)
xticklabels({'', 'Ospedali', 'Studi medici','Gruppi industriali', 'Grossisti', 'Farmacie'})
ylabel('Come')
yticks(-1:0.5:1.5)
yticklabels({'Acido ialuronico', 'Ricerca e sviluppo'})
zlabel('Cosa')
zticks(-1:0.5:1.5)
%zticklabels({'', 'Salute articolare', 'Riparazione tissutale', 'Oftalmologia', 'Ginecologia', 'Medicina Estetica', 'Neuroscienze', 'Cardiovascolare', '', '', '', ''})
title('Modello di Abell')
grid on
