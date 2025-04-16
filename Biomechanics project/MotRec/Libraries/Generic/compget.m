function data = compget(Path, File, varname)
%
% data = compget(filename, varname)
%
% Print names and get data for all variables with varname as a substring.
% Argument "filename" is the name of a COMPAMM result file.
% Signals are returned as columns. 
%


[~, Extension] = getFilenameAndExt(File);
checkFileAndPath(Path,File,{Extension});

[d,i] = compread([Path,File]);
inds = fstrings(varname,i);
data = [];

[m,n] = size(d);

for j=1:length(inds),
  if n >= 3,
    %disp([i(inds(j),:), num2str(d(inds(j),1:3)), ' ...']);
  elseif n >= 2,
    %disp([i(inds(j),:), num2str(d(inds(j),1:2)), ' ...']);
  else
    %disp([i(inds(j),:), num2str(d(inds(j),1))]);
  end
  data = [data, d(inds(j),:)'];
end

