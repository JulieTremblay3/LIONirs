function  [SDPairs,is_srs] = Srs_n2SDpairs(srs_n)
% Description : permet d'associer un numéro de source ou détecteur du format d'IOmtg
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

if srs_n > 999999     % détecteurs
    SDPairs = srs_n/1000000;
    is_srs = 2;
elseif srs_n < 999999 % sources
    srsA = mod(srs_n,1000);
    srsB = floor(srs_n/1000);
    %garder le numeros correspondant a la source impair
    if mod(srsA,2)== 1
        SDPairs = ceil(srsA/2);
    elseif mod(srsB,2)==1
        SDPairs = ceil(srsB/2);
    end
    is_srs = 1;
else
    is_srs = 0;
end