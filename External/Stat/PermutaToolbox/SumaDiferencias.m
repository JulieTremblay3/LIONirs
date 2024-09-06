function [Tij,T_deriv,T_time,T_global] = SumaDiferencias(Datos,v1,v2,ValAbs)

%S1 = sumave(Datos,(1:n1));
%S2 = sumave(Datos,(1:n2));
S1 = SumMat2(Datos,v1);
S2 = SumMat2(Datos,v2);
Tij = (S1 - S2);
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

%[d,dim_deriv,dim_time] = size(Datos);
%if dim_deriv == 1
%   T_deriv = Tij;
%else
%   T_deriv = max(Tij')';
%   T_deriv = max(Tij);
%end

%if(dim_time == 1)   
%   T_time = Tij;
%else
 %  T_time = (max(Tij))';
%   T_time = (max(Tij'))';
%end

%T_global   = max([T_time ; T_deriv]);
 