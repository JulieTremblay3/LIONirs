function [Tij,T_deriv,T_time,T_global] = TStudentDep(Datos,ValAbs)

[nx1,nx2,nx3] = size(Datos);
%[M,S] = meanstd2(Datos,(1:nx1));
[M,S] = MeanSigmaMat(Datos,(1:nx1));
Tij = M./sqrt(S/nx1);
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

   



