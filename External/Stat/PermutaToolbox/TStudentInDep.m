function [Tij,T_deriv,T_time,T_global] = TStudentInDep(Datos,v1,v2,ValAbs)

n1 = length(v1);
n2 = length(v2);
%[M1,S1] = meanstd2(Datos,v1);% Media y Sigma del Nuevo Grupo1
%[M2,S2] = meanstd2(Datos,v2);% Media y Sigma del Nuevo Grupo2
[M1,S1] = MeanSigmaMat(Datos,v1);% Media y Sigma del Nuevo Grupo1
[M2,S2] = MeanSigmaMat(Datos,v2);% Media y Sigma del Nuevo Grupo2
Tij = (M1 - M2)./((sqrt(((n1-1)*S1+(n2-1)*S2)/(n1+n2-2)))*(sqrt(1/n1+1/n2)));
%Tij = (M1 - M2)./((sqrt(((n1-1)*S1+(n2-1)*S2)/(n1+n2-2)))*2);
% figure;plot(Tij(:,1))
if (ValAbs == 1)
   Tij = abs(Tij);
end
Tij = Tij';

[d,dim_deriv,dim_time] = size(Datos);
if dim_deriv == 1
   T_deriv = Tij;
else
   T_deriv = max(Tij')';
end

if(dim_time == 1)   
   T_time = Tij';
else
   T_time = (max(Tij))';
end

T_global   = max([T_time ; T_deriv]);
