function [H, MI, S] = second2HMS (elapsedtime)

% SECONDS2HMS returns the input time in seconds
% in format hours:minutes:seconds
%
%   [H, MI, S] = second2HMS (elapsedtime)
%   Inputs:
%     + elapsedtime is the time in seconds
%   Outputs:
%     + H is double variable with the hour. If its value 
%       is zero this means that the time in less than 1 hour
%     + MI is double variable with the minutes. If its value 
%       is zero this means that the time in less than 1 min.
%     + S is double variable with the seconds.


H = elapsedtime/3600;
if H<1
    H=0;
    if elapsedtime>60
        S  = rem(elapsedtime,60);  
        MI = (elapsedtime -S)/60;
        S  = round(S); 
    elseif elapsedtime > 10
        MI = 0;
        S  = round(elapsedtime);
    else
        MI = 0;
        S  = elapsedtime;
    end
else
    H=(elapsedtime-rem(elapsedtime,3600))/3600;
    S = rem(elapsedtime,3600);
    if S<60
        MI=0;
    else
        Sprima=S;
        S  = rem(Sprima,60);
        MI =(Sprima-S)/60;
        S  = round(S); 
    end
end