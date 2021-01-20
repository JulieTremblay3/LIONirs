function [FSupSup,FSupDeriv,FSupTime,FUniv,T0ij] = TestPermut2Grupos(ESTAD,INDEP,DatGrup1,DatGrup2,NPERM)

%global T0 
%
%  ESTAD - Statistic to use
%         1 - t-student with "original" values
%         2 - sum of difference with "original" values
%         4 - t-student with absolute values
%         5 - sum of difference with absolute values
% 
%  INDEP- Indice if groups are independent 
%         0 - DEPENDENT 
%         1 - INDEPENDENT
%
%  DatGrupo1 - Group1 data matrix : 3D --> (subjects x ROI/channel x time point/interval) 
%  DatGrupo2 - Group2 data matrix : 3D --> (subjects x ROI/channel x time point/interval) 
%
%  NPERM  - Number of permutations to conduct
%%OUTPUT
%   FSupSup : global significant probability that the mean is equal to zero 
%   FSupDeriv: vector (1:numROI) of significant probabilities that the mean is equal to zero 
%       for each ROI/channel with respect of all the time points (multivariate)
%   Fsuptime: vector (1:numtime) of significant probabilities that the mean is equal to zero 
%       for each time point/interval with respect of all ROI/channels (multivariate)
%   Funiv : matriz (numROI x numtime) of significant probabilities that the mean is equal 
%       to zero for each univariate variable (for each combination of ROIxtime)  
%


[nx1,ny1,nz1] = size(DatGrup1);
[nx2,ny2,nz2] = size(DatGrup2);
Dat = cat(1,DatGrup1,DatGrup2);



switch ESTAD
   case {1,4}
      if ESTAD ==1
         ValAbs = 0;
      else
         ValAbs = 1;
      end
      switch INDEP
      case 0   %t-student , caso dependiente
         tmpdat = DatGrup1 - DatGrup2;
         [T0ij,T0_deriv,T0_time,T0_global] = TStudentDep(tmpdat,ValAbs);
      case 1   %t-student , caso independiente
         v1 = (1:nx1)';
         v2 = nx1 + (1:nx2)';
         [T0ij,T0_deriv,T0_time,T0_global] = TStudentInDep(Dat,v1,v2,ValAbs);
         % FIN --------  t-student , caso independiente   ----------
      end
%%%%%%%%%%%%%%%%%%%%%%%%%  SUMA DE DIFERENCIAS  %%%%%%%%%%%%%%%%%%   
   case {2,5}
      if ESTAD ==2
         ValAbs = 0;
      else
         ValAbs = 1;
      end
      switch INDEP
      case {0,1}
         v1 = (1:nx1)';
         v2 = nx1 + (1:nx2)';
         [T0ij,T0_deriv,T0_time,T0_global] = SumaDiferencias(Dat,v1,v2,ValAbs);
      end
            
end
% t = -104:4:596;
% plot(t,T0ij)

% 
[FSupSup,FSupDeriv,FSupTime,FUniv] = Permutaciones(ESTAD,INDEP,Dat,nx1,nx2,T0_global,T0_deriv,T0_time,T0ij,NPERM,ValAbs);
FSupSup = 1 - FSupSup;
FSupDeriv = 1 - FSupDeriv ;
FSupTime  = 1 - FSupTime ;
FUniv= 1- FUniv;
% 
% 
