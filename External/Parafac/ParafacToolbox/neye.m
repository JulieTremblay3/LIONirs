function G=neye(Fac);
% NEYE  Produces a super-diagonal array
%
%function G=neye(Fac);
%
% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 1.00 $ Date 5. Aug. 1998 $ Not compiled $
%
% This algorithm requires access to:
% 'getindxn'
%
% See also:
% 'parafac' 'maxvar3' 'maxdia3'
%
% ---------------------------------------------------------
%             Produces a super-diagonal array
% ---------------------------------------------------------
%	
% G=neye(Fac);
%
% Fac      : A row-vector describing the number of factors
%            in each of the N modes. Fac must be a 1-by-N vector. 
%            Ex. [3 3 3] or [2 2 2 2]



% Copyright, 1998 - 
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. Furthermore, the
% code can not be made part of anything but the 'N-way Toolbox'.
% In case of doubt, contact the holder of the copyrights.
%
% Claus A. Andersson
% Chemometrics Group, Food Technology
% Department of Food and Dairy Science
% Royal Veterinary and Agricultutal University
% Rolighedsvej 30, DK-1958 Frederiksberg, Denmark
% E-mail claus@andersson.dk

N=size(Fac,2);
if N==1,
   fprintf('Specify ''Fac'' as e vector to define the order of the core, e.g.,.\n')
   fprintf('G=eyecore([2 2 2 2])\n')
end;

G=zeros(Fac(1),prod(Fac(2:N)));

for i=1:Fac(1),
   [gi,gj]=getindxn(Fac,ones(1,N)*i);
   G(gi,gj)=1;
end;

G = reshape(G,Fac);