function [MX,SX] = SumMat2(X,v)


MX = sum(X(v,:,:),1);
MX = squeeze(MX)';

% Esta fue implementada  para testar media contra cero
%function [MX,SX] = SumMat(X,v)

%[sx,sy,sz]=size(X);
%ind=find(X(:) < 0);
%X(ind) = 0;
%MX = squeeze(sum(X,1));

% MX = zeros(sy, sz);
% for j=1:sy
%   for k=1:sz
%     MX(j,k)=0;  
%     for i=1:sx
%      if (X(i,j,k)>0) 
%        MX (j,k)= MX (j,k)+ X(i,j,k);
%       end;
%     end;
%   end;
% end,
% MX = squeeze(MX)';

