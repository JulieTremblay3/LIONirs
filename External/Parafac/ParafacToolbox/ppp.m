function AB=ppp(A,B);

% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
%
% Copyright, 1998 - 
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. Furthermore, the
% code can not be made part of anything but the 'N-way Toolbox'.
% In case of doubt, contact the holder of the copyrights.
%
% Rasmus Bro
% Chemometrics Group, Food Technology
% Department of Food and Dairy Science
% Royal Veterinary and Agricultutal University
% Rolighedsvej 30, DK-1958 Frederiksberg, Denmark
% Phone  +45 35283296
% Fax    +45 35283245
% E-mail rb@kvl.dk
%
% The parallel proportional profiles product - triple-P product
% For two matrices with similar column dimension the triple-P product
% is ppp(A,B) = [kron(B(:,1),A(:,1) .... kron(B(:,F),A(:,F)]
% 
% AB = ppp(A,B);
%
% NB. This file is obsolete. Use kr.m instead but not that it takes
% inputs oppositely

%
% Copyright 1998
% Rasmus Bro
% KVL,DK
% rb@kvl.dk


disp('PPP.M is obsolete and will be removed in future versions. ')
disp('use KR.M instead. Note that kr(B,A) = ppp(A,B)')

[I,F]=size(A);
[J,F1]=size(B);

if F~=F1
   error(' Error in ppp.m - The matrices must have the same number of columns')
end

AB=zeros(I*J,F);
for f=1:F
   ab=A(:,f)*B(:,f).';
   AB(:,f)=ab(:);
end