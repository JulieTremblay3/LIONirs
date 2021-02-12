function [FSupSup,FSupDeriv,FSupTime,FUniv] = Permutaciones(estadigrafo,INDEP,Dat,n1,n2,T0_global,T0_deriv,T0_time,T0ij,NPERM,ValAbs)

%global VG1 VG2 
%global TP

[nx1,ny1,nz1] = size(Dat);

mit = nx1/2;
v1 = (1:n1)';
v2 = n1 +(1:n2)';

FSupSup = 0;
dimderiv = length(T0_deriv);
dimtime  = length(T0_time);
FSupDeriv = zeros(dimderiv,1);
FSupTime = zeros(dimtime,1);

switch estadigrafo
case {1,2,4,5}
   FUniv = zeros(size(T0ij));
end

V1 = [];
%h = waitbar(0,'Computing Permutations...'); %aussi ligne 102 et 104
disp('Computing Permutations...')
for iperm = 1:NPERM
   
   %%   Calculo del vector de permutaciones
   switch INDEP
   case 0
      [vg1,vg2] = GenVecPermDep(mit);
   case 1
      [vg1,vg2] = GenVecPermInDep(n1,n2);
   end
   
   V1 = [V1,vg1];
   
   switch(estadigrafo)
   case {1,4}
      if estadigrafo ==1
         ValAbs = 0;
      else
         ValAbs = 1;
      end
      switch INDEP
      case 0
         tmpdat = Dat(vg1,:,:) - Dat(vg2,:,:);
         [TP,TP_deriv,TP_time,TP_global] = TStudentDep(tmpdat,ValAbs);
      case 1
         [TP,TP_deriv,TP_time,TP_global] = TStudentInDep(Dat,vg1,vg2,ValAbs);
      end
      [d,dim_deriv,dim_time] = size(Dat);
      d1 = dim_deriv; d2 = dim_time;
      %Distribucion del SUPSUP
      FSupSup = FSupSup + (TP_global < T0_global);
      %Distribucion del SUPj
      FSupDeriv = FSupDeriv + (TP_deriv < T0_deriv);
      %Distribucion del SUPi
      FSupTime = FSupTime + (TP_time < T0_time);
       FUniv = FUniv + ( TP < T0ij);
      
      %F1 = F1 + (TP_global < T0ij);
      %if((d1==1) | (d2==1))
       %  F2 = F2 + ( TP_deriv  < T0ij );
        % F3 = F3 + ( TP_time   < T0ij);
        %else
       %  F2 = F2 + ( (TP_deriv)*ones(1,d2)  < T0ij);
        % F3 = F3 + ( ones(d2,1)*(TP_time')  < T0ij);
        %end
      
      
   case {2,5}
      if estadigrafo ==2
         ValAbs = 0;
      else
         ValAbs = 1;
      end
      
      [TP,TP_deriv,TP_time,TP_global] = SumaDiferencias(Dat,vg1,vg2,ValAbs);
      
      [d,dim_deriv,dim_time] = size(Dat);
      d1 = dim_deriv; d2 = dim_time;
      
      %Distribucion del SUPSUP
      FSupSup = FSupSup + (TP_global < T0_global);
      %Distribucion del SUPj
      FSupDeriv = FSupDeriv + (TP_deriv < T0_deriv);
      %Distribucion del SUPi
      FSupTime = FSupTime + (TP_time < T0_time);
      FUniv = FUniv + ( TP < T0ij);
      
      %F1 = F1 + (TP_global < T0ij);
      %if((d1==1) | (d2==1))
        % F2 = F2 + ( TP_deriv < T0ij );
         %F3 = F3 + ( TP_time  < T0ij);
         %else
         %F2 = F2 + ( (TP_deriv)*ones(1,d2)  < T0ij);
         %F3 = F3 + ( ones(d2,1)*(TP_time')  < T0ij);
         %end
         
end

%waitbar(iperm/NPERM,h)      
end
%close(h)

%disp(V1)   %%%%   Para ver las permutaciones generadas
% bar(V1(6,:))

% FSupSup  = (FSupSup )/(NPERM );
% FSupDeriv  = (FSupDeriv  )/(NPERM );
% FSupTime  = (FSupTime )/(NPERM );
% FUniv  = (FUniv )/(NPERM );


%Correción propuesta por Pesarin   (puede estar con o sin corrección)
FSupSup  = (FSupSup + 0.5 )./(NPERM + 1);
FSupDeriv  = (FSupDeriv + 0.5 )./(NPERM + 1);
FSupTime  = (FSupTime + 0.5 )./(NPERM + 1);
FUniv  = (FUniv+ 0.5 )./(NPERM+ 1 );
