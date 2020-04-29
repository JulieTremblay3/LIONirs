function st=stdnan(X);

%STDNAN estimate std with NaN's
%
% Estimates the standard deviation of each column of X
% when there are NaN's in X.
%
% Columns with only NaN's get a standard deviation of zero


% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
% 
%
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
% E-mail: claus@andersson.dk

[I,J]=size(X);

st=[];
for j=1:J
  id=find(~isnan(X(:,j)));
  if length(id)
    st=[st std(X(id,j))];
  else
    st=[st 0];
  end
end