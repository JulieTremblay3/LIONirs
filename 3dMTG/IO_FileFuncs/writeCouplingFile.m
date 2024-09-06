function [  ] = writeCouplingFile( FileName, HelmetStruct )
Helmet = HelmetStruct.m_Helmet;

intensityVector = zeros(1,32);
listMatch = [18001 20003 22005 24007 26009 28011 30013 32015]; % correpond aux codes dans IOMtg pour chaque émetteur ou détecteur, ici les 8 premiers émetteurs a1b2...
listMatch = [listMatch 17002 19004 21006 23008 25010 27012 29014 31016]; % émetteurs 9 à 16 a2b1...
listMatch = [listMatch 50033 52035 54037 56039 58041 60043 62045 64047]; % émetteurs 17 à 24 c1d2...
listMatch = [listMatch 49034 51036 53038 55040 57042 59044 61046 63048]; % émetteurs 25 à 32 c2d1...
listMatch = [listMatch 1000000 2000000 3000000 4000000 5000000 6000000 7000000 8000000]; % détecteurs 1 à 8 (A à H)
listMatch = [listMatch 9000000 10000000 11000000 12000000 13000000 14000000 15000000 16000000]; % détecteurs 9 à 16 (I à P), les 16 autres détecteurs ne sont pas instantiables...

HelmetHoles = Helmet.v_Holes;
Src = Helmet.Mtg_Data.v_pSrc;
Det = Helmet.Mtg_Data.v_pDet;
HolesMtg = Helmet.Mtg_Data.v_HolesMtg;
[z sizeHoles] = size(HolesMtg);
listHolesED = zeros(1,64);

for i=1:sizeHoles
    if HolesMtg(i)~=0
        indice = find(listMatch == HolesMtg(i));
        listHolesED(indice) = i; % liste des numéros de trous pour chaque émetteur puis détecteurs (16 numéros, 0 si aucun rien)
    end
end

matDist = 1000.*ones(32,32);

for j1=1:32 % pour chaque émetteur
    if listHolesED(j1) ~= 0
        for j2=33:64 % pour chaque détecteur
            if listHolesED(j2) ~= 0
                e0 = HelmetHoles(listHolesED(j1)).Coord;
                d0 = HelmetHoles(listHolesED(j2)).Coord;
                matDist(j1,j2-32) = sqrt((e0.x-d0.x)^2 + (e0.y-d0.y)^2 +  (e0.z-d0.z)^2 );
            end
        end
    end
end

couplingMat = zeros(32,4);

for j3=1:32
    distVector = sort(matDist(j3,:));
    for j4=1:4
        listFind = find(matDist(j3,:) == distVector(j4));
        couplingMat(j3,j4) = listFind(1);
    end    
end

% for j5=1:32
%     couplingMat(j5,:) = sort(couplingMat(j5,:));
% end

fileID = fopen(FileName,'w');

for k=1:32
    fprintf(fileID, 'E%d D%d D%d D%d D%d WL1: %d WL2: %d\n', k, couplingMat(k,1), couplingMat(k,2), couplingMat(k,3), couplingMat(k,4), intensityVector(k), intensityVector(k));
end

fclose(fileID);

end

