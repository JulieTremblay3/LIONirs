function NIRS = nirs_boxy_associate_sources_detectors(NIRS,dist_max,dist_min,mat_Mtg,nb_det,...
                    mat_Ele,pEle,HolesEle,MaxSources)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create arrays and cells SD

%SD.Pairs: rows of 1xN cells, where rows: sources identified by their number and 
%in each cell are the detectors associated with the source that satisfy conditions 
%in other words: SD_pairs(1) = [1 5 6];  --> indicates Source position #1
%was seen by detectors 1,5,and 6'

%SD.PosLasers: matrix where rows are sources and 2 columns: column 1 is source
% for wavelength 830 nm and column 2 is source number for wavelength 690 nm 

%SD.SrcPos, SD.DetPos: source and detector positions 

%Task 1:

%       Specify which measurements to use:
%       SD_pairsSrc#=[Det# # #];
%       i.e. SD_pairs1=[1 5 6];  --> indicates Source position #1 was seen by detectors 1,5,and 6
%       List all detectors which pair with each source OPTODE position

%on utilise la convention suivante pour associer les sources :
% 1 à 8 correspond aux muxs a1b2, a3b4 à a15b16 
% 9 à 16 correspond aux sources b1a2, b3a4, b15a16
% 17 à 24 correspond aux sources c1d2...
% 25 à 32 correspond aux sources d1c2...
% 33 à 40 correspond aux sources e1f2
% 41 à 48 correspond aux sources f1e2
% 49 à 56 correspond aux sources f1h2
% 57 à 64 correspond aux sources h1f2
% a1b2 = SD_pairs(1)
% a3b4 = SD_pairs(2)
% a5b6 = SD_pairs(3)
% a7b8 = SD_pairs(4)
% a9b10 = SD_pairs(5)

%Task 2:
%       Define which lasers were in each source position
%       PosLasers(1) =[#of 830nm laser #of 690nm laser  ]
%       each sources is define by a code 1 to 16 egal a1 to a16 and
%       code 17 to 32 correspond at b1 to b16.

%Task 4: [x,y,z] coordinates for each optode (in cm)

%Task 3: Specify the probe geometry
%i.e. PosLasers(1,:)=[1 2 3 4];  --> laser 1 is 830nm, laser 2 is 755nm etc

nb_srs = 2*MaxSources-1; %63
%SD.Pairs{ceil(nb_srs/2)}=0; %need SD.pairs to be cells because the number 
%detectors associated to each source varies
for srs = 1:2:nb_srs %for each odd numbered source
    [row,col] = find_coor(mat_Mtg(:,:,4),srs);  
    if row ~= 0
        %Task 1
        [det, dist]= find_srs_det(mat_Mtg,srs,dist_max,dist_min,nb_det);
        NIRS.Pairs{ceil(srs/2),1} = det;
        NIRS.Dist{ceil(srs/2),1} = dist*100; %distance source detector
        %Task 2 for wavelengths illuminating in turn 690-830 nm
        %PROPRE AU SYSTEME MUX 32 de boxy
%         if NIRS.nb_Mux==32
            srs1 = mod(mat_Mtg(row,col,4),1000);
            srs2 = floor(mat_Mtg(row,col,4)/1000);
            lambda =  find_Lambda(srs1); 
            if lambda == 830 %odd lines in Boxy files
                NIRS.PosLasers{ceil(srs/2),1}=[srs1,srs2];
            elseif lambda == 690 %even lines in Boxy files  
                NIRS.PosLasers{ceil(srs/2),1}=[srs2,srs1];
            end 
%         elseif NIRS.nb_Mux==16
%             srs1 = mod(mat_Mtg(row,col,4),1000);
%             srs2 = floor(mat_Mtg(row,col,4)/1000);
%             lambda =  find_Lambda16(srs1);             
%             NIRS.PosLasers{ceil(srs/2),1}= lambda;            
%         end
        
        
        
        %Task 4: [x,y,z] coordinates for each source (in cm)
        NIRS.SrcPos(ceil(srs/2),:) = [mat_Mtg(row,col,1)*100 ...
           mat_Mtg(row,col,2)*100 mat_Mtg(row,col,3)*100];
    %else %perhaps incorrect, cell might need to be void rather than 0
    %    SD.Pairs{ceil(srs/2)} = 0; 
    %    SD.PosLasers{ceil(srs/2)}=[0 0];
    %    SD.SrcPos(ceil(srs/2),:) = [0 0 0];
    end
end

%Task 4: [x,y,z] coordinates for each detector (in cm) 
for iDet = 1:nb_det
    [row,col] = find_coor(mat_Mtg(:,:,4), iDet*1000000);
    if row ~= 0  
       NIRS.DetPos(iDet,:) = [mat_Mtg(row,col,1)*100 ...
           mat_Mtg(row,col,2)*100 mat_Mtg(row,col,3)*100];
    else
       NIRS.DetPos(iDet,:) = [0 0 0]; %OK, some detector might be unused
    end   
end

%Task 5: add the electrodes (positions and names)
for iEle = 1:size(mat_Ele,1)
    NIRS.ElePos(iEle,:) = [mat_Ele(iEle,1,1)*100 mat_Ele(iEle,1,2)*100 ...
        mat_Ele(iEle,1,3)*100]; %Electrodes Positions
    NIRS.EleNames{iEle,1} = HolesEle{pEle(iEle)}; %Electrode Names, not necessarily in the right order!
end

end

function [row,col] = find_coor(mat_Mtg, int_Source)
if int_Source > 999999
    [row,col]=find(mat_Mtg==int_Source);
else
    [row,col]=find( (mod(mat_Mtg,1000)==int_Source) + ...
        (floor(mat_Mtg/1000)==int_Source));
end  
end 

function [det, dist] = find_srs_det(mat_Mtg, srs, dist_max, dist_min,nb_det)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Description : Trouve les détecteurs valide près de la source indiquée
%                 en fonction de la distance min et max spécifiée
%   Input: mat_Mux
%          srs      -> code de la source recherche
%          dist_max -> la source est associe au detecteur si < dist_max 
%
%   Output: det -> detecteur valide qui est associe a la source
%
%   Author: Julie Tremblay
%   Date:          01/18/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[rS, cS] = find_coor(mat_Mtg(:,:,4),srs);
det = [];
dist = [];
if isscalar(rS)
    for int_det = 1:nb_det
        [rD ,cD] = find_coor(mat_Mtg(:,:,4),int_det*1000000);
        if isscalar(rD)
          distance = ((mat_Mtg(rS,cS,1)-mat_Mtg(rD,cD,1))^2+...
              (mat_Mtg(rS,cS,2)-mat_Mtg(rD,cD,2))^2+...
              (mat_Mtg(rS,cS,3)-mat_Mtg(rD,cD,3))^2)^0.5;
            if (distance < dist_max && distance > dist_min)
                    det = [det,int_det];
                    dist = [dist, distance];
            end
        end
    end
end
end

function [lambda] = find_Lambda(int_Src)
% Description : Retourne la longueur d'onde de la source utilisée
% Variable :numéro de la source mux32 
% Concepteur: Julie Tremblay 

if int_Src <= 16
    lambda = 830;
elseif int_Src >= 17 && int_Src <= 32
    lambda = 690;
elseif int_Src >= 33 && int_Src <= 48
    lambda = 830;
elseif int_Src >= 49 && int_Src <= 64
    lambda = 690;
elseif int_Src >=65 && int_Src <=80
    lambda = 830;
elseif int_Src >= 81 && int_Src <=96
    lambda = 690;
elseif int_Src >= 97 && int_Src <=112
    lambda = 830;
elseif int_Src >= 113 && int_Src <=128
    lambda = 690;    
end 
end
   
% function [lambda] = find_Lambda16(int_Src)
% % Description : Retourne directement la position de la source dans le fichier boxy 
% % utilisé en mux 16
% % Concepteur: Julie Tremblay 
% 
% if int_Src == 1 | int_Src == 33 %a1b2 %c1d2
%     lambda = [1,2];
% elseif int_Src == 3 | int_Src == 35  %a3b4 
%     lambda = [3,4];
% elseif int_Src == 5 | int_Src == 37  %a5b6
%     lambda = [5,6];
% elseif int_Src == 7 | int_Src == 39
%     lambda = [7,8];
% elseif int_Src == 9 | int_Src == 41
%     lambda = [9,10];
% elseif int_Src == 11 | int_Src == 43
%     lambda = [11,12];
% elseif int_Src == 13 | int_Src == 45
%     lambda = [13,14];
% elseif int_Src == 15 | int_Src == 47
%     lambda = [15,16];   
% elseif int_Src == 2 | int_Src == 34 %a2b1 %c2d1
%     lambda = [2,1];
% elseif int_Src == 4 | int_Src == 36  %a3b4 
%     lambda = [4,3];
% elseif int_Src == 6 | int_Src == 38  %a5b6
%     lambda = [6,5];
% elseif int_Src == 8 | int_Src == 40
%     lambda = [8,7];
% elseif int_Src == 10 | int_Src == 42
%     lambda = [10,9];
% elseif int_Src == 12 | int_Src == 44
%     lambda = [12,11];
% elseif int_Src == 14 | int_Src == 46
%     lambda = [14,13];
% elseif int_Src == 16 | int_Src == 48
%     lambda = [16,15];  
% end
% end