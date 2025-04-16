function [Filename, Extension] = getFilenameAndExt(FilenameExt)

% GETFILENAMEANDEXT returns Filename and Extension of a file.
% The filename may have several points. The extension is consider
% the last letters (any number) after the last point.
%   e.g. "This.is.a.long.name.final"
%        Filename = "This.is.a.long.name"
%        Extension = "final"
%

% Flip name from left to right. Then extension is first from left
FilenameExt_flip = fliplr(FilenameExt);

% Get Extension length
Ext_flip       = strtok(FilenameExt_flip,'.');
Ext_Length     = length(Ext_flip); % Extension number of characters

% Get filename and extension
Filename  = FilenameExt(1:end-(Ext_Length+1));
Extension = FilenameExt(end-(Ext_Length-1):end);

