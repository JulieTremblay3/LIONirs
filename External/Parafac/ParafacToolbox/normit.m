function [Y,D]=normit(X)
%NORMIT normalize
%
%[Y,D]=normit(X);
%Normalizes X, such that X  = Y *D
%It also applies that    X' = D'*Y'
%Always give X column wise


% Copyright
% Claus A. Andersson 1996-
% Chemometrics Group, Food Technology
% Department of Food and Dairy Science
% Royal Veterinary and Agricultutal University
% Rolighedsvej 30, DK-1958 Frederiksberg, Denmark
% E-mail: claus@andersson.dk

[a b]=size(X);

Y=zeros(a,b);

SS=sqrt(sum(X.^2));
for i=1:b,
  Y(:,i)=X(:,i)./SS(i);
end;

if nargout==2,
  D=(Y'*Y)\(Y'*X);
end;

