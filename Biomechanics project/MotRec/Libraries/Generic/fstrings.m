function inds = fstrings(s,strlist)
%
% inds = fstrings(s,strlist)
%
% Return indices of strings in strlist that contains substring s.
%

[n,m] = size(strlist);
inds = [];
for i=1:n,
  ind = findstr(strlist(i,:),s);
  if ~isempty(ind),
    inds = [inds, i];
  end
end
