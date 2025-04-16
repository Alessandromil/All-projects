function [data,index] = compread(file_name)
%
% [data,index] = compread(file_name)
%
% Reads the data and the name index of a COMPAMM result file.
%
f = fopen(file_name,'r');

% Get number of variables from first line:
ln = fgetl(f);
[t,r] = strtok(ln);
n = str2num(r);

% Skip next ( = INDEX_BEGIN)
ln = fgetl(f);

% Skip next ( = TIME)
ln = fgetl(f);

index = 'Time';
for i=1:n-1,
  ln = fgetl(f);
  [s1,s2] = strtok(ln,'"');
  index = str2mat(index,s1);
end

% Skip next ( = INDEX_END)
ln = fgetl(f);

data = fscanf(f,'%f',[n,inf]);

fclose(f);

