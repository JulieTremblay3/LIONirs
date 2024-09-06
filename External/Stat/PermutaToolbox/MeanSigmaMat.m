function [MX,SX] = MeanSigmaMat(X,v)


MX = mean(X(v,:,:),1,'omitnan');			%%Se suma a travez de la primera dimension
MX = squeeze(MX)';

SX = std(X(v,:,:),0,1,'omitnan').^2;    %%El 0 indica que se normaliza respecto a N-1
SX = squeeze(SX)';


      
         
