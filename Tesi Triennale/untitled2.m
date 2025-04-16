% prova pi greco

clear; clc; close all;

counter = 0;
counter2 = 0;
approx = 1000;
delta =  1000;

while delta > 0.01
    x = rand(1, 2) - 0.5;
    counter = counter + 1;
    if sqrt((x(1))^2 + (x(2))^2) <= 1
        counter2 = counter2 + 1;
    end
        approx = 0.25*counter2/counter;
        delta = abs(approx - pi)
end
