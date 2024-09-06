function  Srs_n = SDpairs2Srs_n(SDPairs,is_srs);
% Description : permet d'associer un numéro de source ou détecteur du
% format SDPairs a une sources ou detecteur du format IOmtg
% avec un numéro de format homer permet de reconnaitre s'il s'agit d'une
% source ou d'un détecteur
% Convention suivante pour associer les sources :
%     1 à 8 correspond aux muxs a1b2, a3b4 à a15b16 
%     9 à 16 correspond aux sources b1a2, b3a4, b15a16
%     17 à 24 correspond aux sources c1d2...
%     25 à 32 correspond aux sources d1c2...
%     a1b2 = SD_pairs1
%     a3b4 = SD_pairs2
%     a5b6 = SD_pairs3
%     a7b8 = SD_pairs4
%     a9b10 = SD_pairs5
%     a11b12 = SD_pairs6
% Convention pour associer les détecteurs
% srs = 1000000 => SD_Pairs = 1
% is_srs = 0 si ce n'est ni une source ni un détecteur
%        = 1 si c'est une source
%        = 2 si c'est un détecteur
if is_srs == 2  % les détecteurs
      Srs_n = SDPairs*1000000;
elseif is_srs == 1  % les sources   
    temp = SDPairs * 2 - 1;
    if temp < 16
        Srs_n = (temp + 17)*1000 + temp; 
    elseif temp > 16 && temp < 32
        Srs_n = (temp - 15) + temp*1000;
    elseif temp > 32 && temp < 48
        Srs_n = (temp + 17)*1000 + temp;
    elseif temp > 48 && temp < 64
        Srs_n = (temp - 15) + temp*1000;
    elseif temp >=65 && temp <=80
        Srs_n = (temp + 17)*1000 + temp; 
    elseif temp >= 81 && temp <= 96
        Srs_n = (temp - 15) + temp*1000;
    elseif temp >= 97 && temp <=112
        Srs_n = (temp + 17)*1000 + temp;
    elseif temp >=113 && temp <=128
        Srs_n = (temp - 15) + temp*1000;       
    end
elseif is_srs == 0
    Srs_n = 0;
end
