function [s] = triggerthreshold(adc, taille_min) % adc = canal d'adc provenant de ADCACC, taille_min = nombre de points des vecteurs t, d... pour savoir sur combien de points mettre s
%FUNCTION JEROME 
curveDigital = [];
threshold = 2800;

curveDigital = (adc>threshold);
curveDigital = double(curveDigital);

s1 = size(adc);

curveFinal = zeros(1,s1);

for i=31:s1
    point = curveDigital(i) == 1 && not(ismember(1,curveDigital(i-2:i-1)));
    curveFinal(i) = point;
end

curveFinal = curveFinal';

s = zeros(taille_min,1);
Trigg_uns = find(curveFinal == 1);
Trigg_adapt = Trigg_uns.*taille_min/s1(1);

for i = 1:size(Trigg_adapt)
    s(floor(Trigg_adapt(i)))=1;
end