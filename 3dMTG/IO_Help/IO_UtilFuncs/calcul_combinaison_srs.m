function calcul_combinaison_srs(oHelmet)
% Avoir une distribution automatique des sources et détecteurs
% Previous label gauche on essai de mettre de A
% droite on met C

% On defini la distance entre a1 c1 de contamination       
       answer = (inputdlg('Enter maximal distance','User define',1,{'6'}))
       distancemax = str2num(answer{1})/100;
       sMtg = get_Mtg(oHelmet);
       vHoles = get_vHoles(oHelmet);       
       nbsrc = numel(sMtg.v_pSrc);
       nbdet = numel(sMtg.v_pDet);      
       
   
       %Calcul des distance entre toute les sources et détecteurs
       distancemat = zeros(nbsrc,nbdet); %initialisation ligne srs, colonne det
       for ind_srs = 1:nbsrc
           psrs = sMtg.v_pSrc(ind_srs);
            for ind_det = 1:nbdet
                 pdet = sMtg.v_pDet(ind_det);
                 distancemat(ind_srs,ind_det) = find_distance(vHoles(pdet).Coord.x,vHoles(pdet).Coord.y,vHoles(pdet).Coord.z,...
                    vHoles(psrs).Coord.x,vHoles(psrs).Coord.y,vHoles(psrs).Coord.z);
            end
       end         
       
           
       near = distancemat(:,:) < distancemax ; %on determine quel source à proche de quel détector
    
       
       % 1 test pour chaque source on regarde quelle autre sources ne cause
       % pas de contamination on établie pour chaque source une liste de
       % sources sans danger de contamination         

       
       %liste des combinaisons sans contamination pour chaque source
       listsrs = cell(nbsrc,1); 
       totcombinaison2 = 0; 
       for ind_srs = 1:nbsrc
           list_srs_nocontamination = [];
           ind_near = find(near(ind_srs,:));
           if ~isempty(ind_near)
               for ind_srs_check = 1:nbsrc
                     if isempty(find(near(ind_srs_check,ind_near)))
                         list_srs_nocontamination = [list_srs_nocontamination ind_srs_check ];
                     end
               end
               totcombinaison2 = totcombinaison2 + numel(list_srs_nocontamination);
               listsrs{ind_srs} = list_srs_nocontamination;               
           end         
       end        
       hwaitbar = waitbar(0)
       
       %2 em test pour chaque combination source non contaminer (suite au test1) entre elle
       %on regarde si une troisieme source existe ne causant pas de
       listsrs2 = cell(totcombinaison2,1 ); %initialisation listsrs2 distance sans contamination maximum
       lst_srs2_nocontamination = [];       %liste des sources utilisés
       listsrs1bin = zeros(nbsrc,1);  
       listsrslevel2 = [];
       totcombinaison3 = 0;
       %contamination entre elles
       indcombination= 1;
       for ind_src = 1:nbsrc
           waitbar(ind_src/nbsrc*1/4,hwaitbar)
           listpossible = listsrs{ind_src};
           listsrs1bin(listpossible)=1;           
            for ind_srs_check = 1:numel(listpossible)
                srs_ok1 = listpossible(ind_srs_check );
                ind_srsok1 = listsrs{srs_ok1};
                listsrs2bin = zeros(nbsrc,1);
                listsrs2bin(ind_srsok1) = 1;               
                 % Ce qui ne présente pas de contabination srs1 & srs2   
                 if isempty(listsrslevel2) 
                  listsrs2{indcombination} =  find(listsrs1bin & listsrs2bin);
                         listsrslevel2 =  [listsrslevel2;ind_src srs_ok1];
                         indcombination = indcombination + 1;
                         totcombinaison3 = totcombinaison3 + numel(find(listsrs1bin & listsrs2bin));  
                 end                 
                if ~isempty(listsrslevel2) %eviter les doublons %verifier si la combinaison inverse n'est pas déjà là 1 15 ou 15 1
                    testsrsok1=(listsrslevel2(:,1)==srs_ok1|listsrslevel2(:,2)==srs_ok1);
                    testsrsok2=(listsrslevel2(:,1)==ind_src|listsrslevel2(:,2)==ind_src );
                    if find(testsrsok1&testsrsok2) ;       
                      
                    else
                    	 listsrs2{indcombination} =  find(listsrs1bin & listsrs2bin);
                         listsrslevel2 =  [listsrslevel2;ind_src srs_ok1];
                         indcombination = indcombination + 1;
                         totcombinaison3 = totcombinaison3 + numel(find(listsrs1bin & listsrs2bin));  
                    end
                end
                                        
            end
            listsrs1bin = zeros(nbsrc,1);
       end
       indproblem = []
       for i=1:size(listsrslevel2,1)
          tot = near(listsrslevel2(i,1),:)+...
           near(listsrslevel2(i,2),:);       
           if ~isempty(find(tot>1))
            indproblem = [indproblem,i];
           end       
       end
       listsrslevel2(indproblem,:)=[];
       
       %On ajuste la taille de listsrs2 sans les doublons       
       %3 em test pour chaque 3 sources non contaminer existe-t-il une
       %quatrieme source ok             
       listsrs3 = cell(totcombinaison3,1 ); 
       listsrs1bin = zeros(nbsrc,1);
       totcombinaison4 = 0;
       indcombination= 1;
       listsrslevel3 =[];
       for ind_src = 1:size(listsrslevel2,1)
           waitbar(ind_src/size(listsrslevel2,1)*0.25+ 0.25,hwaitbar)
           listpossible = listsrs2{ind_src};
           listsrs1bin(listpossible)=1;
           listsrs2bin = zeros(nbsrc,1);
           for ind_srs_check = 1:numel(listpossible) %on regarde les choix voir ceux qui fonctionne
               srs_ok2 = listpossible(ind_srs_check);
               ind_srsok2 = listsrs{srs_ok2};
               listsrs2bin = zeros(nbsrc,1);
               listsrs2bin(ind_srsok2) = 1;
               % Ce qui ne présente pas de contabination srs1 & srs2    
               if isempty(listsrslevel3)
                    indcombination = indcombination + 1;
                    totcombinaison4 = totcombinaison4 + numel(find(listsrs1bin & listsrs2bin));
                    listsrs3{indcombination} =  find(listsrs1bin & listsrs2bin);              
                    listsrslevel3 =  [listsrslevel3;listsrslevel2(ind_src,:) srs_ok2];
               end
               srs_ok1 = listsrslevel2(ind_src,1);
               srs_ok3 = listsrslevel2(ind_src,2);
                if ~isempty(listsrslevel3) %eviter les doublons %verifier si la combinaison inverse n'est pas déjà là 1 15 ou 15 1
                    testind_srs1 = listsrslevel3(:,1)==srs_ok2|listsrslevel3(:,2)==srs_ok2|listsrslevel3(:,3)==srs_ok2;
                    testind_srs2 = listsrslevel3(:,1)==srs_ok1|listsrslevel3(:,2)==srs_ok1|listsrslevel3(:,3)==srs_ok1;
                    testind_srs3 = listsrslevel3(:,1)==srs_ok3|listsrslevel3(:,2)==srs_ok3|listsrslevel3(:,3)==srs_ok3;
                        if find(testind_srs2 &  testind_srs1 & testind_srs3 );   
                            
                        else
                            indcombination = indcombination + 1;
                            totcombinaison4 = totcombinaison4 + numel(find(listsrs1bin & listsrs2bin));
                            listsrs3{indcombination} =  find(listsrs1bin & listsrs2bin);              
                            listsrslevel3 =  [listsrslevel3;listsrslevel2(ind_src,:) srs_ok2];
                      end
                 end                
           end
            listsrs1bin = zeros(nbsrc,1);
       end
       %Double check distance
       indproblem = []
       for i=1:size(listsrslevel3,1)
          tot = near(listsrslevel3(i,1),:)+...
           near(listsrslevel3(i,2),:)+...
           near(listsrslevel3(i,3),:);          
           if ~isempty(find(tot>1))
            indproblem = [indproblem,i];
           end       
       end
       listsrslevel3(indproblem,:)=[];
       
       %4 em test pour définir les combinaisons de 4 sources sans
       %contamination         
           listsrs4 = cell(totcombinaison4,1 );
           listsrs1bin = zeros(nbsrc,1);
           indcombination= 1;
           listsrslevel4 =[];
           for ind_src = 1:size(listsrslevel3,1)
               waitbar(ind_src/size(listsrslevel3,1)*0.5+ 0.50,hwaitbar)
               listpossible = listsrs3{ind_src};
               listsrs1bin(listpossible)=1;
               listsrs2bin = zeros(nbsrc,1);
               if ~isempty(listpossible)
                   for ind_srs_check = 1:numel(listpossible) %on regarde les choix voir ceux qui fonctionne
                       srs_ok2 = listpossible(ind_srs_check);
                       ind_srsok2 = listsrs{srs_ok2};
                       listsrs2bin = zeros(nbsrc,1);
                       listsrs2bin(ind_srsok2) = 1;
                      % listsrs4{indcombination} =  find(listsrs1bin & listsrs2bin); 
                       indcombination = indcombination + 1;
                       listsrslevel4 =  [listsrslevel4;listsrslevel3(ind_src,:) srs_ok2];
                   end
               end
               listsrs1bin = zeros(nbsrc,1);
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
       indproblem = []
       for i=1:size(listsrslevel4,1)
          tot = near(listsrslevel4(i,1),:)+...
           near(listsrslevel4(i,2),:)+...
           near(listsrslevel4(i,3),:)+...
           near(listsrslevel4(i,4),:);
           if ~isempty(find(tot>1))
            indproblem = [indproblem,i];
           end       
       end
       listsrslevel4(indproblem,:)=[];
       
       %Fin calcul des possibilité
       [name,path] = uiputfile([num2str(distancemax*100),'.com'],'savepossibility matrix');
       save([path, name],'listsrslevel4','listsrslevel3','listsrslevel2');
    
       close(hwaitbar)