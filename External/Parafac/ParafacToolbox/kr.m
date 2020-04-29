function AB = kr(A,B);
%KR Khatri-Rao product
%
% The Khatri - Rao product
% For two matrices with similar column dimension the khatri-Rao product
% is kr(A,B) = [kron(B(:,1),A(:,1) .... kron(B(:,F),A(:,F)]
% 
% I/O AB = ppp(A,B);
%
% kr(A,B) equals ppp(B,A) - where ppp is the triple-P product = 
% the parallel proportional profiles product which was originally 
% suggested in Bro, Ph.D. thesism, 1998

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
% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $

[I,F]=size(A);
[J,F1]=size(B);

if F~=F1
   error(' Error in kr.m - The matrices must have the same number of columns')
end

AB=zeros(I*J,F);
for f=1:F
   ab=B(:,f)*A(:,f).';
   AB(:,f)=ab(:);
end