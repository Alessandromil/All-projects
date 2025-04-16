% demo
clear all %clears all the variables in the workspace
close all %closes all the figures that were open
clc %clears the command window

% load data
load data.mat %loads the data

% plot
figure(1) %opens the figure
subplot(2,1,1) %activates the first figure in a 2-row 1-column figure
plot(x,y1) %plot function 1
xlabel('x(m)') %writes the x-label
ylabel('y(m)') %writes the y-label
axis([0 20 0 400]) %sets the values for the x- and y-axes.
title('Y1 function') %writes the title for the figure
subplot(2,1,2)
plot(x,y2)
xlabel('x(m)')
ylabel('y(m)')
axis([0 20 0 400])
title('Y2 function')

% compute integrals
integral1 = trapz(x,y1); %computes the integral of function 1.
integral2 = trapz(x,y2);
disp(['Intregral(Y1) = ' num2str(integral1) ' m2, Integral(Y2) = ' num2str(integral2) ' m2.']) %displays the results in the command window

% compute average values
average1 = mean(y1); %computes the average value of function 1.
average2 = mean(y2);
disp(['Mean(Y1) = ' num2str(average1) ' m, Mean(Y2) = ' num2str(average2) ' m.'])