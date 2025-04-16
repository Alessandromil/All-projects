function [apen, sampen]=apsampen(varargin)
%APSAMPEN  Computes Approximate Entropy as defined by Steven Pincus in:
%   S. M. Pincus, "Approximate entropy as a measure of system complexity",
%   in Proc. Natl. Acad. Sci. USA, vol. 88, pp. 2297-2301 (1991).
%   S. M. Pincus and W. M. Huang, "Approximate entropy: statistical
%   properties and applications", Community Statistic Theory and Methods 21,
%   3061-3077 (1992).
%   S. M. Pincus, "Approximate entropy (ApEn) as a complexity measure",
%   Chaos 5, 110-117 (1995).
%
%   Also computes Sample Entropy as defined by
%   Richman JS and Moorman JR, "Physiological time-series analysis using 
%   approximate entropy and sample entropy", Am J Physiol Heart Circ Physiol. 
%   2000 Jun;278(6):H2039-49.  
%   
%   APEN=APSAMPEN(X,M,R) computes approximate entropy of the vector X;
%   M and R are two parameters. M can be considered as a Taken's embedding 
%   dimension, R as a noise rejection level.
%
%   [APEN, SAMPEN]=APSAMPEN(X,M,R) computes both approximate entropy and sample entropy
%   of the vector X;
%
%   [APENE, SAMPEN]=APSAMPEN(X,M,R,FLAG_STD)
%   If FLAG_STD==1, R=R*std(data_vector).
%   If FLAG_STD==0, R=R.
%   By default FLAG_STD is set to 1.
%
% Copyright (C) 2003 Roberto Sassi
% Copyright (C) 2003 Elena Fabani and Roberto Sassi
% 
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version. 
% 
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
% General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License along 
% with this program; if not, write to the Free Software Foundation, Inc., 
% 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. 
% 
% You may contact the author by e-mail (sassi@biomed.polimi.it) or postal 
% mail (Dipartimento di Bioingegneria, Politecnico di Milano, 
% Piazza Leonardo da Vinci 32, 20133 Milano, ITALY).  
%
% Matlab (C) wrapper papen.c
% Revision: 2, (2003) Roberto Sassi <sassi@biomed.polimi.it>
% Revision: 1, (1998) Davide Torti and Roberto Sassi
% Dipartimento di Bioingegneria
% Politecnico di Milano, ITALY
% Last modified: 11 July 2003
% 
% routine apen()
% Revision: 3
% Copyright (c) 1996 Elena Fabani and Roberto Sassi
% Dipartimento di Bioingegneria
% Politecnico di Milano, ITALY
% Last modified: 12 April 1996
% 
% The algorithm apen() exploits a circular array to reduce the number of
% comparisons to N*(N-1)/2-m*(m-1)/2 (being N the number of points in 
% data_vector). It was published in:
% E. Fabani and R. Sassi, "Studio dell'entropia approssimata per la
% classificazione di serie temporali: applicazioni al segnale di variabilita 
% cardiaca", Tesi di Laurea, Politecnico di Milano (1996). (IN ITALIAN).

if( exist('apsampen') ~= 3 ),
    if( exist('apsampen.c') ~= 2 ),
        error('File apsampen.c is missing.');
    else,
        mex apsampen.c
        [apen, sampen]=apsampen(varargin{:});
    end,
end,

