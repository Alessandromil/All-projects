function Index = getCellIndex(Props, Name)

% GETINDEX returns the row index of the string 'Name'
% in the cell 'Props'. The first column of the cell must
% contain strings. The function searches for the row index of
% the string 'Name' in this first column of strings in
% the cell.
%
%   Index = getCellIndex(Props, Name)
%   Inputs:
%     + Props is a cell (nElements x nCols) where nElements
%       and nCols can have any size. The first column of the cell
%       Props{:,1} must contains strings. This column is 
%       where the function looks for the string 'Name'
%     + Name is the string that the function searches in 
%       the cell
%   Outputs:
%     + Indexis the index of the string 'Name' in
%       the cell 'Props' - integer


% Find body properties
nProps = size(Props, 1);
Index = [];

for i = 1 : nProps
    if strcmp(Name, Props{i, 1}) == 1
           Index = i;
    end
end
if isempty(Index)
    disp(['  The  "',Name,'" is not contained in the cell']);    
end
