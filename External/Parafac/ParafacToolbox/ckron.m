function C=ckron(A,B)
%CKRON
% C=ckron(A,B)
%
% Claus Andersson, Jan. 1996

% Should not be compiled to overwrite ckron.mex

[mA,nA] = size(A);
[mB,nB] = size(B);

C = zeros(mA*mB,nA*nB);
if mA*nA <= mB*nB
  for i = 1:mA
  iC = 1+(i-1)*mB:i*mB;
    for j = 1:nA
      jC = 1+(j-1)*nB:j*nB;
      C(iC,jC) = A(i,j)*B;
    end
  end
else
  for i = 1:mB
    iC = i:mB:(mA-1)*mB+i;
    for j = 1:nB
      jC = j:nB:(nA-1)*nB+j;
      C(iC,jC) = B(i,j)*A;
    end
  end
end

%In order to avoid compiling
%break
