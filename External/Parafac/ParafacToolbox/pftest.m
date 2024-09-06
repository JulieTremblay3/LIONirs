function [ssX,Corco,Factors] = pftest(NumRep,X,Fac,Options,const,OldLoad,Weights);

%PFTEST find the number of PARAFAC components
%
% See also:
% 'corcondia'
%
%
% TEST HOW MANY COMPONENTS TO USE IN PARAFAC MODEL
% 
% Used for testing the appropriate number of components in a PARAFAC
% model. Input the appropriate input options (see PARAFAC for help) as
% well as the number of replicate fittings to use (NumRep). For one to Fac 
% component models are fitted each NumRep times.
% Three measures are output in three matrices, where the f'th row hold the values
% for the f-component model and each column correspond to a specific replicate run
%
% ssX:   THE SUM-SQUARED ERROR
%
%        What-To-Look-For:
%        look for sudden change as in a Scree-plot (often difficult)
%        and look for sudden increase in number of local minima (replicate
%        points for one component are not identical). This is often a good
%        indication that noise is being modeled.
%
% Corco: CORE CONSISTENCY DIAGNOSTIC (CORCONDIA)
%
%        What-To-Look-For:
%        CORCONDIA is a percentage below or equal to 100%. A value of 80-100% 
%        means that the model is valid, while a value below, say 40% means that
%        the model is not valid. A value between 40 and 80% means that the model
%        is probably valid but somehow difficult to estimate, e.g., due to 
%        slight misspecification or correlations. The Corcondia will mostly
%        decrease with number of components but very sharply where the correct
%        number of components are exceeded. Hence, the appropriate number of
%        components is the model with the highest number of components and a
%        valid CORCONDIA
%
% I/O FORMAT as PARAFAC with an additional parameter, NumRep, for the number
% of replicate runs
%
% [ssX,Corco] = pftest(NumRep,X,Fac,Options,const,OldLoad,Weights);
%
% Note that plots are only produced in Matlab ver. > 5.2

% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
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
% Rasmus Bro
% Chemometrics Group, Food Technology
% Department of Food and Dairy Science
% Royal Veterinary and Agricultutal University
% Rolighedsvej 30, DK-1958 Frederiksberg, Denmark
% Phone  +45 35283296
% Fax    +45 35283245
% E-mail rb@kvl.dk


if nargin<4
  Options=[];
end
if nargin<5
  const=[];
end
if nargin<6
  OldLoad=[];
end
if nargin<7
  FixMode=[];
end
if nargin<8
  Weights=[];
end

ssX=zeros(Fac,NumRep);
Corco=zeros(Fac,NumRep);


FactorsOut=[];
for f=1:Fac
  Options(2)=0; % DTLD init
  [Factors,it,err,corcondia]=parafac(X,f,Options,const,OldLoad,FixMode,Weights);
  ssX(f,1)=err; 
  BestErr=err;
  Corco(f,1)=corcondia;
  FactorsNow=Factors;
  
  Options(2)=1; % SVD init
  [Factors,it,err,corcondia]=parafac(X,f,Options,const,OldLoad,FixMode,Weights);
  ssX(f,2)=err;
  Corco(f,2)=corcondia;
  if err<BestErr
    FactorsNow=Factors;
    BestErr=err;
  end

  Options(2)=2; % random init
  for rep=3:NumRep
    [Factors,it,err,corcondia]=parafac(X,f,Options,const,FixMode,OldLoad,Weights);
    ssX(f,rep)=err;
    Corco(f,rep)=corcondia;
    if err<BestErr
      FactorsNow=Factors;
      BestErr=err;
    end
  end

  FactorsOut=[FactorsOut;FactorsNow];
end


% Only plot if MATLAB version >= 5
VER = version;
if VER(1)~='4'

delta=.25/NumRep; % shift points a little horizontally for appearance

subplot(2,1,1)
for r=1:NumRep
   plot(1+(r-1)*delta:1:Fac+(r-1)*delta,ssX(:,r), ...
      'MarkerEdgeColor','k','MarkerFaceColor','r', ...
      'LineWidth',2,'Marker','o','LineStyle','none', ...
      'MarkerSize',8)
   hold on
end
axis([1 Fac+1 0 max(ssX(:)) ])
set(gca,'XTick',[1:Fac])
ylabel('Residual sum of squares','FontWeight','bold')
title('PARAFAC TEST','FontWeight','bold')
hold off

subplot(2,1,2)
for r=1:NumRep
   plot(1+(r-1)*delta:1:Fac+(r-1)*delta,Corco(:,r), ...
      'MarkerEdgeColor','k','MarkerFaceColor','r', ...
      'LineWidth',2,'Marker','o','LineStyle','none', ...
      'MarkerSize',8)
   hold on
end
hold off
MinCo=min(Corco(:));
axis([1 Fac+1 min([MinCo 0]) 100 ]);
set(gca,'XTick',[1:Fac])

%title('Corcondia')
ylabel('Core consistency','FontWeight','bold')
xlabel('Number of components','FontWeight','bold')

end

if nargout>2
  disp(' ')
  whic = input([' For which number of components (1:',num2str(Fac),') do you want the parameters : ']);
  Factors=FactorsOut( sum(1:whic-1)*sum(DimX)+1: sum(1:whic)*sum(DimX) );
end