function [i,j]=GetIndxn(R,Idx);
%GETINDXN
%
%[i,j]=GetIndxn(R,Idx)
%
% Copyright
% Claus A. Andersson 1995-
% Chemometrics Group, Food Technology
% Department of Food and Dairy Science
% Royal Veterinary and Agricultutal University
% Rolighedsvej 30, T254
% DK-1958 Frederiksberg
% Denmark
% E-mail: claus@andersson.dk

l=size(Idx,2);

i=Idx(1);
j=Idx(2);

if l==3,
  j = j + R(2)*(Idx(3)-1);
 else
  for q = 3:l,
    j = j + prod(R(2:(q-1)))*(Idx(q)-1);
  end;
end;
