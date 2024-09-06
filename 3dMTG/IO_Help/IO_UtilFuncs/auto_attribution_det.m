function oHelmet = auto_attribution_det(oHelmet)

% Avoir une distribution automatique des sources et d�tecteurs
% Previous label gauche on essai de mettre de A
% droite on met C
       
       sMtg = get_Mtg(oHelmet);
       vHoles = get_vHoles(oHelmet);       
       nbsrc = numel(sMtg.v_pSrc);
       nbdet = numel(sMtg.v_pDet);   

%Voulez-vous placer les d�tecteurs
%On place les d�tecteurs
       %Determiner indice � Gauche et � Droite, on privil�gie le cot� en
       %fonction de la position du d�tecteur
     
           rightind = [];
           leftind =[];
           for ind_det = 1:nbdet
                labeldet= vHoles( sMtg.v_pDet(ind_det)).Label;
                if labeldet(1)=='D'
                    rightind = [rightind,sMtg.v_pDet(ind_det) ];
                elseif labeldet(1)== 'G' 
                    leftind =  [leftind, sMtg.v_pDet(ind_det)];
                elseif labeldet(1)== 'Z' 
                    leftind =  [leftind, sMtg.v_pDet(ind_det)];                
                end
           end

           %On place les d�tecteurs gauche colonne 1 et � droite
           %colonne 2
           listDET = [1000000,3000000;...
                      2000000,40000000;...
                      5000000,7000000;...
                      6000000,8000000;...
                      9000000,11000000;...
                      10000000,12000000;...
                      13000000,15000000;...
                      14000000,16000000];                       

           indlistdet = 1;
           for ind = 1:numel(rightind)
               sMtg.v_HolesMtg(rightind(ind)) = listDET(indlistdet);       
               indlistdet =indlistdet + 1;
           end

           for ind = 1:numel(leftind)
               sMtg.v_HolesMtg(leftind(ind)) = listDET(indlistdet);   
               indlistdet =indlistdet + 1;
           end
       %
       % fin du placement des d�tecteurs priorisations gauche droite
       % du � la position g�ographique des d�tecteurs sur le montage
         oHelmet= set_Mtg(oHelmet,sMtg);