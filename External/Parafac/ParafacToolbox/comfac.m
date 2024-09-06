function [A,B,C,FIT,IT]=comfac(X,Fac,Options,DoComp,CompIt,Init);

%COMFAC Algorithm for fitting the complex-valued PARAFAC model
% 
% See e.g. Rasmus Bro, N. D. Sidiropoulos, and G. B. Giannakis. A Fast 
%          Least Squares Algorithm for Separating Trilinear Mixtures. 
%          Proc. ICA99 ­ Int. Workshop on Independent Component Analysis 
%          and Blind Signal Separation. Jan. 11­15, Aussois, France:289-294, 1999.
%
%	        N. D. Sidiropoulos, G. B. Giannakis, and Rasmus Bro. Blind PARAFAC 
%          Receivers for DS-CDMA Systems. IEEE Transactions on Signal Processing March, 2000.
%
% The algorithm works by first compressing the data using a Tucker3 models. Subsequently the 
% PARAFAC model is fitted to the compressed array, either initialized with DTLD (~ESPRIT) or 
% with PARAFAC-ALS. The solution is de-compressed to the original space, and a few safe-guard
% PARAFAC-ALS steps are performed. Further optimization of this algorithm is possible, e.g.
% if the data are very large, if the profiles are very correlated or if the noise is huge. This
% has not been pursued here, in order to have a generally applicable algorithm
% 
% INPUT
% X       : I x J x K data three-way array 
% Fac     : Number of factors in PARAFAC model
%
% OPTIONAL INPUTS
% Options : A 1x2 vector
% (1)     : Mode to compress to dimension two in
%           DTLD if smallest mode-dimension is more
%           than two. Default is the smallest dimension
% (2)     : Number of extra components in Tucker3
%           compression model. Default if no given 
%           is one
% 
% ADVANCED OPTIONS
% DoComp  : 0      => No compression
%           1      => Compression
% CompIt  : CompIt => Max number of iterations in Compression
% Init    : 0      => GRAM/DTLD
%           1      => Ten x PARAFAC with 10 iterations started from RandOrth
%
% 
% I/O
% [A,B,C,FIT,IT]=comfac(X,Fac,Options,DoComp,CompIt,Init);
% 
% Or short : [A,B,C]=comfac(X,Fac);
%
%
% Copyright, 1998 - 
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. In case of doubt, 
% contact the holder of the copyrights.
%
% Copyright 1998
% Rasmus Bro
% KVL, Denmark, rb@kvl.dk
%      & 
% Nikos Sidiropoulos
% Univ. Minnesota, nikos@ece.umn.edu
% 

% 31/3/00 RB fixed bug when dimensions not suitable for DTLD

crit = 1e-5; % Criterion used for fitting with Alternating LS algorithm

if length(size(X))~=3
   error(' The data must be held in a three-way array')
end
   
DimX = size(X);
X = reshape(X,DimX(1),prod(DimX(2:3)));

dbg=1;

DeafultExtraCompInTucker=2;
if exist('DoComp')~=1|isempty(DoComp)
   DoComp = 1;
end
if exist('CompIt')~=1|isempty(CompIt)
   CompIt = 3;
end
if exist('Init')~=1|isempty(Init)
   Init   = 1;
end
if length(find(DimX>Fac))<2 % GRAM/ESPRIT cannot be used
   ReplaceTLDwithRand = 1;
else
   ReplaceTLDwithRand = 0;
end   

if dbg
%   home
   disp(' ')
   disp([' Fitting ',num2str(Fac),'-comp. PARAFAC model using COMFAC ...'])
   disp([' Array size: ',num2str(DimX(1)),' x ',num2str(DimX(2)),' x ',num2str(DimX(3))])
end

if nargin<3|isempty(Options)
  Options=[0 DeafultExtraCompInTucker];
end
if Options(1)==0
  [out,CompToTwo] = min(DimX);
else
  CompToTwo = Options(1);
end
if length(Options)<2
  Options(2)=DeafultExtraCompInTucker; % Number of additional components in compression as compared to model
end
W = Options(2)+Fac;
W = [W W W];
for i=1:3
  if W(i)>DimX(i)
    W(i)=DimX(i);
  end
end

% Display things
if DoComp & dbg
  disp([' Array compressed to size ',num2str(W(1)),' x ',num2str(W(2)),' x ',num2str(W(3))])
end

% Compression

if DoComp
   [Ut,Vt,Zt,Gt,fit]=tucker3(X,DimX,W,CompIt);
else
   Gt = X;
   W = DimX;  
   Ut=speye(DimX(1));
   Vt=speye(DimX(2));
   Zt=speye(DimX(3));
end

% Initialize PARAFAC using DTLD
if Init == 0
   if ~ReplaceTLDwithRand
      disp(' Initializing using direct trilinear decomposition')
      [Adtld,Bdtld,Cdtld]=cdtld(Gt,W,Fac,CompToTwo);
   else
      disp(' Initializing using random values because tld cannot be used due to the sizes')
      [Adtld,Bdtld,Cdtld,fit]=cparafac(Gt,W,Fac,crit,0,0,0,10);
   end
   
elseif Init == 1
   disp(' Initializing using best of 11 initial small runs')
   if ~ReplaceTLDwithRand
      [Adtld,Bdtld,Cdtld]=cdtld(Gt,W,Fac,CompToTwo);
   else
      [Adtld,Bdtld,Cdtld,fit]=cparafac(Gt,W,Fac,crit,0,0,0,10);
   end
  fitout=sum(sum(abs(  Gt-Adtld*ppp(Bdtld,Cdtld).').^2));
  for rep=1:10
     [adtld,bdtld,cdtldvar,fit]=cparafac(Gt,W,Fac,crit,0,0,0,10);
     if fit<fitout
       Adtld=adtld;Bdtld=bdtld;Cdtld=cdtldvar;fitout=fit;
     end
  end
end

% Fit PARAFAC model in compressed space
if DoComp
  disp(' Fitting PARAFAC in compressed space')
  [Acomp,Bcomp,Ccomp,fit2,it2]=cparafac(Gt,W,Fac,crit,Adtld,Bdtld,Cdtld);
  % Decompress to original domain
  disp(' Transforming interim solution to original space')
  Ainit = Ut*Acomp;
  Binit = Vt*Bcomp;
  Cinit = Zt*Ccomp;
else
  Ainit = Adtld;
  Binit = Bdtld;
  Cinit = Cdtld;
  it2=0;
  fit2=NaN;
end

% Fit PARAFAC model in raw space
disp(' Fitting PARAFAC in original space')
[A,B,C,FIT,IT]=cparafac(X,DimX,Fac,crit,Ainit,Binit,Cinit);
disp(' Algorithm converged')


function [A,B,C,fit,it]=cparafac(X,DimX,Fac,crit,A,B,C,maxit,DoLineSearch);

% Complex PARAFAC-ALS
% Fits the PARAFAC model Xk = A*Dk*B.' + E
% where Dk is a diagonal matrix holding the k'th
% row of C.
%
% Uses on-the-fly projection-compression to speed up 
% the computations. This requires that the first mode 
% is the largest to be effective
% 
% INPUT
% X       : Data
% DimX    : Dimension of X
% Fac     : Number of factors
% OPTIONAL INPUT
% crit    : Convergence criterion (default 1e-6)
% A,B,C   : Initial parameter values
%
% I/O
% [A,B,C,fit,it]=parafac(X,DimX,Fac,crit,A,B,C);
%
% Copyright 1998
% Rasmus Bro
% KVL, Denmark, rb@kvl.dk


% Initialization
if nargin<8
  maxit   = 2500;      % Maximal number of iterations
end
showfit = pi;         % Show fit every 'showfit'th iteration (set to pi to avoid)

if nargin<4
  crit=1e-6;
end

if crit==0
  crit=1e-6;
end

I = DimX(1);
J = DimX(2);
K = DimX(3);

InitWithRandom=0;
if nargin<7
   InitWithRandom=1;
end
if nargin>6 & size(A,1)~=I
  InitWithRandom=1;
end

if InitWithRandom

  if I<Fac
    A = rand(I,Fac);
  else
    A = orth(rand(I,Fac));
  end
  if J<Fac
    B = rand(J,Fac);
  else
    B = orth(rand(J,Fac));
  end
  if K<Fac
    C = rand(K,Fac);
  else
    C = orth(rand(K,Fac));
  end
end

SumSqX = sum(sum(abs(X).^2));
fit    = sum(sum(abs(X-A*ppp(B,C).')));
fit0   = fit;
fitold = 2*fit;
it     = 0;
Delta  = 5;

while abs((fit-fitold)/fitold)>crit & it<maxit & fit>10*eps
   it=it+1;
   fitold=fit;

   % Do line-search
   if rem(it+2,2)==-1
      [A,B,C,Delta]=linesrch(X,DimX,A,B,C,Ao,Bo,Co,Delta);
   end

   Ao=A;Bo=B;Co=C;
   % Update A
   Xbc=0;

   for k=1:K
     Xbc = Xbc + X(:,(k-1)*J+1:k*J)*conj(B*diag(C(k,:)));
   end
   A = Xbc*inv((B'*B).*(C'*C)).';

   % Project X down on orth(A) - saves time if first mode is large
   [Qa,Ra]=qr(A,0);
   x=Qa'*X;

   % Update B
   Xac=0;
   for k=1:K
     Xac = Xac + x(:,(k-1)*J+1:k*J).'*conj(Ra*diag(C(k,:)));
   end
   B = Xac*inv((Ra'*Ra).*(C'*C)).';
   
   % Update C
   ab=inv((Ra'*Ra).*(B'*B));
   for k=1:K
     C(k,:) = (ab*diag(Ra'* x(:,(k-1)*J+1:k*J)*conj(B))).';
   end

   % Calculating fit. Using orthogonalization instead
   %fit=0;for k=1:K,residual=X(:,(k-1)*J+1:k*J)-A*diag(C(k,:))*B.';fit=fit+sum(sum((abs(residual).^2)));end
   [Qb,Rb]=qr(B,0);
   [Z,Rc]=qr(C,0);
   fit=SumSqX-sum(sum(abs(Ra*ppp(Rb,Rc).').^2));

   if rem(it,showfit)==0
     fprintf(' %12.10f       %g        %3.4f \n',fit,it,100*(1-fit/fit0));
   end
end

%fprintf(' %12.10f       %g        %3.4f \n',fit,it,100*(1-fit/fit0));


function [A,B,C]=cgram(X1,X2,F);

% cGRAM - Complex Generalized Rank Annihilation Method
% Fits the PARAFAC model directly for the case of a 
% three-way array with only two frontal slabs.
% For noise-free trilinear data the algorithm is exact.
% 
% INPUTS:
% X1    : I x J matrix of data from observation one
% X2    : I x J matrix of data from observation two
% Fac   : Number of factors
% 
% OUTPUTS:
% A     : Components in the row mode (I x F)
% B     : Components in the column mode (J x F)
% C     : Weights for each slab; C(1,:) are the component 
%         weights for first slab such that the approximation
%         of X1 is equivalent to X1 = A*diag(C(1,:))*B.'
%
% Copyright 1998
% Rasmus Bro, KVL, DK
% rb@kvl.dk


  DontShowOutput = 1;

  [U,s,V]=svd(X1+X2);
  U=U(:,1:F);
  V=V(:,1:F);

%  S2=S1*b
%  b*k = k*l =>
%  S1*b*k = S1*k*l =>
%  S2*k = S1*k*l =>
%  inv(S1*k)*S2*k = l = diagonal

  S1=U'*X1*V;
  S2=U'*X2*V;
  [v,d]=eig(S1\S2);

  ddd=d;
  d=diag(d);out=abs(d)>eps;v=v(:,out);d=d(out); %only significant terms
  [dd,out]=sort(abs(d));out=flipud(out);
  d=d(out);d=diag(d);v(:,out);v=v/norm(v);% sort them

  A = U*S1*v;
  B=V/v';
  B=(B.')';
  C=(pinv(ppp(A,B))*[X1(:) X2(:)]).';

  if ~DontShowOutput
    fit = sum(sum(abs([X1 X2] - [A*diag(C(1,:))*B.' A*diag(C(2,:))*B.']).^2));
    disp([' GRAM fitted data with a sum-squared error of ',num2str(fit)])
 end
 
 
 
 function [A,B,C,fit]=cdtld(X,DimX,F,SmallMode);

% DIRECT TRILINEAR DECOMPOSITION
% calculate the parameters of the three-
% way PARAFAC model directly. The model
% is not the least-squares but will be close
% to for precise data with little model-error
%
% This implementation works with an optimal
% compression using least-squares Tucker3 fitting
% to generate two pseudo-observation matrices that
% maximally span the variation of all samples. per
% default the mode of smallest dimension is compressed
% to two samples, while the remaining modes are 
% compressed to dimension F.
% 
% For large arrays it is fastest to have the smallest
% dimension in the first mode
%
% INPUT
% [A,B,C]=dtld(X,DimX,F);
% X is the I x J x K array unfolded to an I x JK matrix
% DimX = [I J K]
% F is the number of factors to fit
% An optional parameter may be given to enforce which
% mode is to be compressed to dimension two
%
% Copyright 1998
% Rasmus Bro, KVL
% rb@kvl.dk


DontShowOutput = 1;

%rearrange X so smallest dimension is in first mode
if nargin<4
  [a,SmallMode] = min(DimX);
  X = nshape(X,DimX,SmallMode);
  DimX = DimX([SmallMode 1:SmallMode-1 SmallMode+1:3]);
  Fac   = [2 F F];
else
  X = nshape(X,DimX,SmallMode);
  DimX = DimX([SmallMode 1:SmallMode-1 SmallMode+1:3]);
  Fac   = [2 F F];
end

if DimX(1) < 2
  error(' The smallest dimension must be > 1')
end
if any(DimX(2:3)-Fac(2:3)<0)
  error(' This algorithm requires that two modes are of dimension not less the number of components')
end

% Compress data into a 2 x F x F array. Only 10 iterations are used since exact SL fit is insignificant; only obtaining good truncated bases is important
[At,Bt,Ct,G]=tucker3(X,DimX,Fac,10);

% Fit GRAM to compressed data
[Bg,Cg,Ag]=cgram(reshape(G(1,:),F,F),reshape(G(2,:),F,F),F);

% De-compress data and find A
BB = Bt*Bg;
CC = Ct*Cg;
AA = X*pinv(ppp(BB,CC)).';

if SmallMode == 1
  A=AA;
  B=BB;
  C=CC;
elseif SmallMode == 2 
  A=BB;
  B=AA;
  C=CC;
elseif SmallMode == 3
  A=BB;
  B=CC;
  C=AA;
end


if ~DontShowOutput
  fit = sum(sum(abs(X - AA*ppp(BB,CC).').^2));
  disp([' DTLD fitted raw data with a sum-squared error of ',num2str(fit)])
end


function [A,B,C,G,fit,it]=tucker3(X,DimX,W,maxit);

% [A,B,C,G]=tucker3(X,R,W,maxit);
% This is an SVD-based algorithm for finding the
% parameters of the three-way Tucker3 model when
% the array as well as parameters are complex
% 
% Copyright 1998
% Rasmus Bro & Claus A. Andersson
% KVL, Denmark, rb@kvl.dk

UseNIPALS = 0; % Use NIPALS (svdf.m) instead of SVD. It's cheaper in terms of flops for large arrays. For small there's no big difference
if UseNIPALS == 1 & any(imag(X(:)))
   disp(' Apparently the loadings in this NIPALS are a little oblique. Check that (decrease convergence criterion)')
   error(' NIPALS has not yet been changed to handle complex numbers. Use SVD (set UseNIPALS to 0)')
end

DontShowOutput = 0;

if nargin<4
  maxit=100;
end

%Initialising counters and others
SSX=sum(sum(abs(X).^2));
it=0;
Oldfit=1e100;
Diff=1e100;
I=DimX(1);J=DimX(2);K=DimX(3);
B=orth(rand(J,W(2)));
C=orth(rand(K,W(3)));

while Diff>1e-6&it<maxit

  it=it+1;

  %Updating A
    TA1=C'*reshape(X,I*J,K).';
    TA2=B'*reshape(TA1,W(3)*I,J).';
    TA3=reshape(TA2,W(2)*W(3),I).';
    if UseNIPALS
      [TA4 TA5 TA6]=svdf(TA3,W(:,1));
    else
      [TA4 TA5 TA6]=svd(TA3,0);
    end
    A=TA4(:,1:W(1));

  %Updating C
    TC1=reshape(A'*X,W(1)*J,K).';
    TC2=B'*reshape(TC1,W(1)*K,J).';
    TC3=reshape(TC2,K*W(2),W(1)).';
    TC3=reshape(TC3,W(1)*W(2),K).';
    if UseNIPALS
      [TC4 TC5 TC6]=svdf(TC3,W(3));
    else
      [TC4 TC5 TC6]=svd(TC3,0);
    end
    C=TC4(:,1:W(3));

  %Updating B
    TB1=reshape(C'*TC1,W(1)*W(3),J).';
    if UseNIPALS
       [TB2 TB3 TB4]=svdf(TB1,W(2));
    else
       [TB2 TB3 TB4]=svd(TB1,0);
    end
    B=TB2(:,1:W(2));

  %Calculate core & fit
    G1=reshape(A'*X,W(1)*J,K).';
    G2=reshape(C'*G1,W(1)*W(3),J).';
    G=reshape(B'*G2,W(2)*W(3),W(1)).';
    fit=sum(sum(abs(G.^2)));
    fit=SSX-fit;

    Diff=abs(Oldfit-fit);
    Oldfit=fit;

end

function [u,s,v] = svdf(X,F);

% Rand-reduced SVD based on NIPALS (actually 
% the power-method the way it's implemented here)

maxit = 30; % Use 30
crit = 1e-4; % Use 1e-4
[I,J] = size(X);
u = zeros(I,F);
v = zeros(J,F);

if J > I
   x = X*X';
else
   x = X'*X;
end


for f = 1:F
   p = sum(x)';
   converged=0;
   it = 0;
   while ~converged
      it = it +1;
      pold = p;
      p = x*p;
      p = p/norm(p);
      if norm(p-pold)/norm(pold)<crit | it>maxit
         converged = 1;
      end
   end
  
   if J > I
     u(:,f) = p;
     v(:,f) = X'*p;
     s(f) = norm(v(:,f));     
     v(:,f) = v(:,f)/s(f);
   else
     v(:,f) = p;
     u(:,f) = X*p;
     s(f) = norm(u(:,f));
     u(:,f) = u(:,f)/norm(u(:,f));
   end
   x = x - s(f)^2*p*p';
end


function [X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17]=nshape(X,DimX,f);

% $ Version 1.03 $ Date 18. July 1999 $ Not compiled $
% $ Version 1.031 $ Date 18. July 1999 $ Error in help figure and now outputs new DimX $ Not compiled $
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
%
% [Xf,DimXf] = nshape(X,DimX,f);
%
% Refolds an N-way array so that Xf is X with index
% f as row-index, and the remaining in succesive order. For an 
% I x J x K x L four-way array this means X1 is I x JKL, X2 is
% J x ÌKL, X3 is K x IJL, and X4 is L x IJK
%
%
%    K  _______             
%      /      /|           1      J     2·J    J·K
%     /______/ |         1  _____________________
%    |      |  |           |      |      |      |
%    |      | /    -->     |      |      |      |        f = (Mode) 1 (same as original array)
% I  |______|/          I  |______|______|______|
%           J
%
%                          1      I     2·I    K·I
%                        1  _____________________
%                          |      |      |      |
%                  -->     |      |      |      |        f = (Mode) 2
%                        J |______|______|______|
%
%  
%                          1      I     2·I    I·J
%                        1  _____________________
%                          |      |      |      |
%                  -->     |      |      |      |        f = (Mode) 3
%                        K |______|______|______|
%
%
% If the last input is not given all rearrangements are given.
% For a fourway array this would read
% [X1,X2,X3,X4]=nshape(X,DimX);
%
%	Copyright
%	Rasmus Bro & Claus A. Andersson 1995
%	Denmark
%	E-mail rb@kvl.dk

ord=chkpfdim(X,DimX,NaN);
elemen=prod(DimX);
if nargin==2
  do_it=ones(1,ord);
else
  do_it=zeros(1,ord);
  do_it(f)=1;
end

if do_it(1)==1
  X1=X;
end


% _____Make X2_____

if do_it(2)==1
  X2=X(:,1:DimX(2)).';
  for R2=DimX(2)+1:DimX(2):elemen/DimX(1)
    X2=[X2 X(:,R2:R2+DimX(2)-1).'];
  end
end

if ord>3

% _____Make X3 - Xord-1_____

for comp=3:ord-1
  if do_it(comp)==1
    xx=[];  % Denne kan opbygges til flere
    for R2=1:prod(DimX(2:comp-1)):prod(DimX(2:comp))
      x=[];
      for R3=R2:prod(DimX(2:comp)):prod(DimX(2:ord))
        x=[x reshape(X(:,R3:R3+prod(DimX(2:comp-1))-1),1,prod(DimX(1:comp-1)))];
      end % for
      xx=[xx;x];
    end
    eval(['X',num2str(comp),'=xx;']);
  end % for comp
end
end % if ord>3


% _____Make Xord_____

if do_it(ord)==1
  xx=[];
  for R3=1:elemen/(DimX(1)*DimX(ord)):elemen/DimX(1)
    xx=[xx; reshape(X(:,R3:R3+elemen/(DimX(1)*DimX(ord))-1),1,elemen/DimX(ord))];
  end % for
  eval(['X',num2str(ord),'=xx;']);
end

if nargin==3
   eval(['X1=X',num2str(f),';']);
   X2 = [DimX(f) DimX([1:f-1 f+1:ord])];
end

function ord=chkpfdim(X,DimX,show);

% show == NaN => no text
%
% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $

if nargin<3
   show=0;
end


% Find the order, i.e., number of ways
ord=length(DimX);

% Check if DimX corresponds to size of X
if DimX(1)~=size(X,1)|prod(DimX)/DimX(1)~=size(X,2)
  disp(' ')
  disp(' Size of array does not correspond to dimensions given in DimX')
  error(['disp('' The matrix input must be of size ',num2str(DimX(1)),' x ',num2str(prod(DimX)/DimX(1)),' if DimX is correctly given'')'])
end

if ~isnan(show)
  txt=[];
  for i=1:ord-1
    txt=[txt num2str(DimX(i)) ' x '];
  end
  txt=[txt num2str(DimX(ord))];
  disp([' The array is a ',num2str(ord),'-way array with'])
  disp([' dimensions: ' txt])
end

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
% Copyright 1998
% Rasmus Bro
% KVL,DK
% rb@kvl.dk

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

function [NewA,NewB,NewC,DeltaMin]=linesrch(X,DimX,A,B,C,Ao,Bo,Co,Delta);

dbg=0;

if nargin<5
  Delta=5;
else
  Delta=max(2,Delta);
end

dA=A-Ao;
dB=B-Bo;
dC=C-Co;
Fit1=sum(sum(abs(X-A*ppp(B,C).').^2));
regx=[1 0 0 Fit1];
Fit2=sum(sum(abs(X-(A+Delta*dA)*ppp((B+Delta*dB),(C+Delta*dC)).').^2));
regx=[regx;1 Delta Delta.^2 Fit2];

while Fit2>Fit1
  if dbg
    disp('while Fit2>Fit1')
  end
  Delta=Delta*.6;
  Fit2=sum(sum(abs(X-(A+Delta*dA)*ppp((B+Delta*dB),(C+Delta*dC)).').^2));
  regx=[regx;1 Delta Delta.^2 Fit2];
end

Fit3=sum(sum(abs(X-(A+2*Delta*dA)*ppp((B+2*Delta*dB),(C+2*Delta*dC)).').^2));
regx=[regx;1 2*Delta (2*Delta).^2 Fit3];

while Fit3<Fit2
  if dbg
    disp('while Fit3<Fit2')
  end
  Delta=1.8*Delta;
  Fit2=Fit3;
  Fit3=sum(sum(abs(X-(A+2*Delta*dA)*ppp((B+2*Delta*dB),(C+2*Delta*dC)).').^2));
  regx=[regx;1 2*Delta (2*Delta).^2 Fit2];
end

% Add one point between the two smallest fits
[a,b]=sort(regx(:,4));
regx=regx(b,:);
Delta4=(regx(1,2)+regx(2,2))/2;
Fit4=sum(sum(abs(X-(A+Delta4*dA)*ppp((B+Delta4*dB),(C+Delta4*dC)).').^2));
regx=[regx;1 Delta4 Delta4.^2 Fit4];

%reg=pinv([1 0 0;1 Delta Delta^2;1 2*Delta (2*Delta)^2])*[Fit1;Fit2;Fit3]
reg=pinv(regx(:,1:3))*regx(:,4);
%DeltaMin=2*reg(3);

DeltaMin=-reg(2)/(2*reg(3));

%a*x2 + bx + c = fit
%2ax + b = 0
%x=-b/2a

NewA=A+DeltaMin*dA;
NewB=B+DeltaMin*dB;
NewC=C+DeltaMin*dC;
Fit=sum(sum(abs(X-NewA*ppp(NewB,NewC).').^2));

if dbg
  regx
  plot(regx(:,2),regx(:,4),'o'),
  hold on
  x=linspace(0,max(regx(:,2))*1.2);
  plot(x',[ones(100,1) x' x'.^2]*reg),
  hold off
  drawnow
  [DeltaMin Fit],pause
end

[minfit,number]=min(regx(:,4));
if Fit>minfit
  DeltaMin=regx(number,2);
  NewA=A+DeltaMin*dA;
  NewB=B+DeltaMin*dB;
  NewC=C+DeltaMin*dC;
end
