function [FSupSup,FSupDeriv,FSupTime,FUniv,T0ij] = TestPermut2Grupos(ESTAD,INDEP,DatGrup1,DatGrup2,NPERM)

%global T0 
%
%  ESTAD - Estadigrafo a utilizar
%         1 - t-student con Identidad
%         2 - Suma de diferencias con Identidad
%         4 - t-student con Valor Absoluto
%         5 - Suma de diferencias con Valor Absoluto
% 
%  INDEP- Indica si los grupos son independientes o no
%         0 - No son independientes
%         1 - Son independientes
%
%  DatGrupo1 -Datos del grupo 1 matrix de tres dimensiones (numero de
%  sujetos por numero de regiones por numero de tiempos
%  DatGrupo2 -Datos del grupo 2matrix de tres dimensiones (numero de
%  sujetos por numero de regiones por numero de tiempos
%
%  NPERM  - Numero de Permutaciones a realizar  
%%Salida
%FSupSup :probabilidad global de significacion si la media es igual a cero
%FSupDeriv: vector (1:numROI) de probabilidades de significacion si la media es cero en cada
%ROI con respecto a todos los times multivariadamente
%Fsuptime: vector (1:numtime) de probabilidades de significacion si la media es cero en cada
%time con respecto a todos las deriv multivariadamente
%Funiv : matriz (numROI x numtime) de probabilidades de significacion si la media es cero en cada
%variable univariadamente%  
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
