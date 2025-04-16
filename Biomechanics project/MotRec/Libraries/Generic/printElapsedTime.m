function printElapsedTime(H, MI, S, FileID, MessageString)

% PRINTELAPSEDTIME prints a message and the time (hours:min:second)
% in the Matlab command window and in a text file (*.log)
%
%   printElapsedTime(H, MI, S, FileID, MessageString)
%   Inputs:
%     + H is double variable with the hour. If its value 
%       is zero this means that the time in less than 1 hour
%     + MI is double variable with the minutes. If its value 
%       is zero this means that the time in less than 1 min.
%     + S is double variable with the seconds. 
%     + FileID is a scalar MATLAB integer, called a file 
%       identifier. You use the FileID as the first argument 
%       to other file input/output routines.
%     + MessageString is the string to be printed on the file
%       and the Matlab command window together with the time.
%   Outputs:
%     NONE


% defines quantity of information to be displayed


if MI==0 & H==0
    str = sprintf([MessageString, num2str(S,'%2.0f'), ' second(s)']);
%     fprintf(FileID, '%s\n', str);
%     if INFO_TYPE == 1
         disp(str);
%     elseif INFO_TYPE == 0
%         disp(['  ',str]);
%      end

elseif H==0
    str = sprintf([MessageString, num2str(MI,'%2.0f'), ' minute(s)  ', num2str(S), ' second(s)']);
%     fprintf(FileID, '%s\n', str);
%     if INFO_TYPE == 1
         disp(str);
%     elseif INFO_TYPE == 0
%         disp(['  ',str]);
%     end

else
    str = sprintf([MessageString, num2str(H,'%2.0f'), ' hour(s)  ', num2str(MI), ' minute(s)  ', num2str(S), ' second(s)']);
%     fprintf(FileID, '%s\n', str);
%     if INFO_TYPE == 1
         disp(str);
%     elseif INFO_TYPE == 0
%         disp(['  ',str]);
%     end

end
