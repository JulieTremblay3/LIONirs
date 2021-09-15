%Description : Utilise les valeurs de concentration a un temps donnee et
%produit un colormap sur le cortex ou la peau. Projection radial, moyenne les canaux
%superposés.
function PrjStruct = display_MRIcolor(PrjStruct,PMI,d1,type)
%type 0  = skin
%type 1  = cortex

global currentsub
cf = PMI{currentsub}.currentFile;
oMRI = get_MRI_Data(PrjStruct);
if type == 0 %Peau
    [VertexBuffer, IndexBuffer] = get_SkinMesh( oMRI );
    if numel(d1)~=(size(PMI{currentsub}.data(cf).MeasList,1)/2) %detect vcolor display
        vcolor =  d1(1:size(VertexBuffer,1));
        oMRI = set_SkinVcolor(oMRI,vcolor);
        PrjStruct = set_MRI_Data(PrjStruct, oMRI);
        return
    end
end
if type == 1 %Cortex
    %NOUVEAU DISPLAY SUR LA PEAU ET ENSUITE PROJECTION SUR CORTEX
    [VertexBuffer, IndexBuffer] = get_CortexMeshLowRes( oMRI );
    if numel(d1)~=(size(PMI{currentsub}.data(cf).MeasList,1)/2) %detect vcolor display
        vcolor =  d1(1:size(VertexBuffer,1));
        oMRI = set_CortexLowResVcolor(oMRI,vcolor);
        PrjStruct = set_MRI_Data(PrjStruct, oMRI);
        return
    end
    [VertexBuffer, IndexBuffer] = get_SkinMesh( oMRI ); %On prend la peau pour commencer

end


meanchannel = 1;
VertexBuffer(:,4) = 1;
T = get_matTransform(oMRI);   % Transfert de plan pour le mapping
Vertex_tmp = VertexBuffer*T;
[m n] = size(VertexBuffer);
Vcolor = zeros(m,1);
DispHelm = get_Helmet(PrjStruct);
sMtg = get_Mtg( DispHelm );
vHoles = get_vHoles( DispHelm );
cf = PMI{currentsub}.currentFile;
ml = PMI{currentsub}.data(cf).MeasList;
pzone =[];
tic
iscase0 = 1;
iscase1 = 1;
iscase2 = 1;
iscase3 = 1;
%adaptation pour plan 2d spherical angle work well on 3d sphere as head but
%not on 2d plane,use only x y step projection. 
if iscase0
for ch = 1:size(ml,1)/2 %Passage 1 coordonné 
    % Coordonnées du canal à tracer
    if PMI{currentsub}.data(cf).MeasListAct(ch)& ml(ch,1)~= 0;
        Srs_n = SDpairs2Srs_n(ml(ch,1),1);
        Det_n = SDpairs2Srs_n(ml(ch,2),2);
        pSrs = find( Srs_n == sMtg.v_HolesMtg);
        pDet = find( Det_n == sMtg.v_HolesMtg);
        if ~isempty(pSrs)&~isempty(pDet)
            if type == 0 | type == 1  %Ne donne pas de bon resultats avec les normals sur le cortex remplacé par la projection de
            %par une projection sur la peau de l'activation et ensuite
            %transfert sur le cortex.
                x_srs = vHoles(pSrs).Coord.x-vHoles(pSrs).Normal.x*vHoles(pSrs).SkinDepth;
                y_srs = vHoles(pSrs).Coord.y-vHoles(pSrs).Normal.y*vHoles(pSrs).SkinDepth;
                z_srs = vHoles(pSrs).Coord.z-vHoles(pSrs).Normal.z*vHoles(pSrs).SkinDepth;
                x_det = vHoles(pDet).Coord.x-vHoles(pDet).Normal.x*vHoles(pDet).SkinDepth;
                y_det = vHoles(pDet).Coord.y-vHoles(pDet).Normal.y*vHoles(pDet).SkinDepth;
                z_det = vHoles(pDet).Coord.z-vHoles(pDet).Normal.z*vHoles(pDet).SkinDepth;
           
          if  z_srs == 0 & z_det ==0
             delta_x = (x_det - x_srs)/7;
             delta_x =x_srs :((x_det - x_srs)/7):x_det;
             delta_y =y_srs :((y_det - y_srs)/7):y_det;
          else
              break
          end
                       pzone = [];
            for i = 1:(numel(delta_x)-1)
                x_min = min([delta_x(i),delta_x(i+1)])-.002;;
                x_max = max([delta_x(i),delta_x(i+1)])+.002;;
                y_min = min([delta_y(i),delta_y(i+1)])-.002;;
                y_max = max([delta_y(i),delta_y(i+1)])+.002;;
                p_in = find(Vertex_tmp(:,1) > x_min & Vertex_tmp(:,1) < x_max & Vertex_tmp(:,2) > y_min & Vertex_tmp(:,2) < y_max);
                pzone = [pzone; p_in];
                if ~isempty(p_in)
                    for p_vertex = 1 : numel(p_in)
                        if ~(Vcolor(p_in((p_vertex)),1) == 0)
                            if meanchannel == 0
                                Vcolor(p_in(p_vertex),1) = (Vcolor(p_in(p_vertex),1) + d1(ch))/2; %moyennage des canaux superposés
                            elseif 1
                                val = [Vcolor(p_in(p_vertex),1),d1(ch)];
                                [maxval,id]=max(abs(val));
                                Vcolor(p_in(p_vertex),1) = val(id);            %la plus grande valeur est selectionné minimal ou maximal
                            else
                                Vcolor(p_in(p_vertex),1) = (Vcolor(p_in(p_vertex),1) + d1(ch)); %Somme des canaux superposés
                            end
                        else
                            Vcolor(p_in(p_vertex),1) = d1(ch);
                        end
                    end
                end
            end
            end
            %no 3 d
            iscase1 = 0;
            iscase2 = 0;
            iscase3 = 0;
        end
    end
end
end
%L'arrière est moins bien couvert changement de plan pour eviter l'effet
%des coordonnées sphérique ou phi est moins bien représenté dans l'axe
%z on effectue un changement de plan pour éviter ce défaut
%CAS1 x,y,z standard
[Vertex_sph(:,1),Vertex_sph(:,2),Vertex_sph(:,3)] = cart2sph(Vertex_tmp(:,1),Vertex_tmp(:,2),Vertex_tmp(:,3));
if iscase1
for ch = 1:size(ml,1)/2 %Passage 1 coordonné 
    % Coordonnées du canal à tracer
    if PMI{currentsub}.data(cf).MeasListAct(ch)& ml(ch,1)~= 0;
        Srs_n = SDpairs2Srs_n(ml(ch,1),1);
        Det_n = SDpairs2Srs_n(ml(ch,2),2);
        pSrs = find( Srs_n == sMtg.v_HolesMtg);
        pDet = find( Det_n == sMtg.v_HolesMtg);
        if ~isempty(pSrs)&~isempty(pDet)
            if type == 0 | type == 1  %Ne donne pas de bon resultats avec les normals sur le cortex remplacé par la projection de
            %par une projection sur la peau de l'activation et ensuite
            %transfert sur le cortex.
                x_srs = vHoles(pSrs).Coord.x-vHoles(pSrs).Normal.x*vHoles(pSrs).SkinDepth;
                y_srs = vHoles(pSrs).Coord.y-vHoles(pSrs).Normal.y*vHoles(pSrs).SkinDepth;
                z_srs = vHoles(pSrs).Coord.z-vHoles(pSrs).Normal.z*vHoles(pSrs).SkinDepth;
                x_det = vHoles(pDet).Coord.x-vHoles(pDet).Normal.x*vHoles(pDet).SkinDepth;
                y_det = vHoles(pDet).Coord.y-vHoles(pDet).Normal.y*vHoles(pDet).SkinDepth;
                z_det = vHoles(pDet).Coord.z-vHoles(pDet).Normal.z*vHoles(pDet).SkinDepth;
            end
            [theta_srs,phi_srs,r_srs]= cart2sph(x_srs,y_srs,z_srs);
            [theta_det,phi_det,r_det]= cart2sph(x_det,y_det,z_det);
            if theta_srs * theta_det < 0
                if theta_srs < -pi/2
                    theta_srs = theta_srs + 2*pi;
                elseif theta_det < -pi/2
                    theta_det = theta_det + 2*pi;
                end
            end
            if phi_srs * phi_det < 0
                if  phi_srs < -pi/2
                    phi_srs = phi_srs + 2*pi;
                elseif phi_det < -pi/2
                    phi_det = phi_det + 2*pi;
                end
            end
            if phi_srs * phi_det == 0
                if  phi_srs < -pi/2
                    phi_srs = phi_srs + 2*pi;
                elseif phi_det < -pi/2
                    phi_det = phi_det + 2*pi;
                end
            end
            
            delta_theta = theta_det-theta_srs;
            delta_phi = phi_det-phi_srs;
            %Définition de la zone à tracer en coordonnée polaire
            d_theta = theta_srs:(delta_theta/7):theta_det ;
            d_phi =  [phi_srs:delta_phi/7:phi_det];
            %Attribution de la couleur au voxel dans la zone du canal
            pzone = [];
            for i = 1:(numel(d_theta)-1)
                x_min = min([d_theta(i),d_theta(i+1)])-.03;
                x_max = max([d_theta(i),d_theta(i+1)])+.03;
                y_min = min([d_phi(i),d_phi(i+1)])-.03;
                y_max = max([d_phi(i),d_phi(i+1)])+.03;
                p_in = find(Vertex_sph(:,1) > x_min & Vertex_sph(:,1) < x_max & Vertex_sph(:,2) > y_min & Vertex_sph(:,2) < y_max);
                pzone = [pzone; p_in];
                if ~isempty(p_in)
                    for p_vertex = 1 : numel(p_in)
                        if ~(Vcolor(p_in((p_vertex)),1) == 0)
                            if meanchannel == 0
                                Vcolor(p_in(p_vertex),1) = (Vcolor(p_in(p_vertex),1) + d1(ch))/2; %moyennage des canaux superposés
                            elseif 1
                                val = [Vcolor(p_in(p_vertex),1),d1(ch)];
                                [maxval,id]=max(abs(val));
                                Vcolor(p_in(p_vertex),1) = val(id);            %la plus grande valeur est selectionné minimal ou maximal
                            else
                                Vcolor(p_in(p_vertex),1) = (Vcolor(p_in(p_vertex),1) + d1(ch)); %Somme des canaux superposés
                            end
                        else
                            Vcolor(p_in(p_vertex),1) = d1(ch);
                        end
                    end
                end
            end
        end
    end
end
end
%CAS2 x,z,y standard
if iscase2
[Vertex_sph(:,1),Vertex_sph(:,2),Vertex_sph(:,3)] = cart2sph(Vertex_tmp(:,1),Vertex_tmp(:,3),Vertex_tmp(:,2));
for ch = 1:size(ml,1)/2 %Passage 1 coordonné 
    % Coordonnées du canal à tracer
    if PMI{currentsub}.data(cf).MeasListAct(ch)& ml(ch,1)~= 0;
        Srs_n = SDpairs2Srs_n(ml(ch,1),1);
        Det_n = SDpairs2Srs_n(ml(ch,2),2);
        if ch==57
            1
        end
        pSrs = find( Srs_n == sMtg.v_HolesMtg);
        pDet = find( Det_n == sMtg.v_HolesMtg);
        if ~isempty(pSrs)&~isempty(pDet)
            if type == 0 | type == 1  %Ne donne pas de bon resultats avec les normals sur le cortex remplacé par la projection de
            %par une projection sur la peau de l'activation et ensuite
            %transfert sur le cortex.
                x_srs = vHoles(pSrs).Coord.x-vHoles(pSrs).Normal.x*vHoles(pSrs).SkinDepth;
                z_srs = vHoles(pSrs).Coord.y-vHoles(pSrs).Normal.y*vHoles(pSrs).SkinDepth;
                y_srs = vHoles(pSrs).Coord.z-vHoles(pSrs).Normal.z*vHoles(pSrs).SkinDepth;
                x_det = vHoles(pDet).Coord.x-vHoles(pDet).Normal.x*vHoles(pDet).SkinDepth;
                z_det = vHoles(pDet).Coord.y-vHoles(pDet).Normal.y*vHoles(pDet).SkinDepth;
                y_det = vHoles(pDet).Coord.z-vHoles(pDet).Normal.z*vHoles(pDet).SkinDepth;
            end
            [theta_srs,phi_srs,r_srs]= cart2sph(x_srs,y_srs,z_srs);
            [theta_det,phi_det,r_det]= cart2sph(x_det,y_det,z_det);
            if theta_srs * theta_det < 0
                if theta_srs < -pi/2
                    theta_srs = theta_srs + 2*pi;
                elseif theta_det < -pi/2
                    theta_det = theta_det + 2*pi;
                end
            end
            if phi_srs * phi_det < 0
                if  phi_srs < -pi/2
                    phi_srs = phi_srs + 2*pi;
                elseif phi_det < -pi/2
                    phi_det = phi_det + 2*pi;
                end
            end
            delta_theta = theta_det-theta_srs;
            if delta_theta <2
              
            delta_phi = phi_det-phi_srs;
            %Définition de la zone à tracer en coordonnée polaire
            d_theta = theta_srs:(delta_theta/7):theta_det ;
            d_phi =  [phi_srs:delta_phi/7:phi_det];
            %Attribution de la couleur au voxel dans la zone du canal
            pzone = [];
            for i = 1:(numel(d_theta)-1)
                x_min = min([d_theta(i),d_theta(i+1)])-.03;
                x_max = max([d_theta(i),d_theta(i+1)])+.03;
                y_min = min([d_phi(i),d_phi(i+1)])-.03;
                y_max = max([d_phi(i),d_phi(i+1)])+.03;
                p_in = find(Vertex_sph(:,1) > x_min & Vertex_sph(:,1) < x_max & Vertex_sph(:,2) > y_min & Vertex_sph(:,2) < y_max);
                pzone = [pzone; p_in];
                if ~isempty(p_in)
                    for p_vertex = 1 : numel(p_in)
                        if ~(Vcolor(p_in((p_vertex)),1) == 0)
                            if meanchannel == 0
                                Vcolor(p_in(p_vertex),1) = (Vcolor(p_in(p_vertex),1) + d1(ch))/2; %moyennage des canaux superposés
                            elseif 1
                                val = [Vcolor(p_in(p_vertex),1),d1(ch)];
                                [maxval,id]=max(abs(val));
                                Vcolor(p_in(p_vertex),1) = val(id);            %la plus grande valeur est selectionné minimal ou maximal
                            else
                                Vcolor(p_in(p_vertex),1) = (Vcolor(p_in(p_vertex),1) + d1(ch)); %Somme des canaux superposés
                            end
                        else
                            Vcolor(p_in(p_vertex),1) = d1(ch);
                        end
                    end
                end
            end
            else
                1
            end
        end
    end
end
end
%CAS3 z,y,x standard
if iscase3 == 1
[Vertex_sph(:,1),Vertex_sph(:,2),Vertex_sph(:,3)] = cart2sph(Vertex_tmp(:,3),Vertex_tmp(:,2),Vertex_tmp(:,1));
for ch = 1:size(ml,1)/2 %Passage 1 coordonné 
    % Coordonnées du canal à tracer
    if PMI{currentsub}.data(cf).MeasListAct(ch)& ml(ch,1)~= 0;
        Srs_n = SDpairs2Srs_n(ml(ch,1),1);
        Det_n = SDpairs2Srs_n(ml(ch,2),2);
        pSrs = find( Srs_n == sMtg.v_HolesMtg);
        pDet = find( Det_n == sMtg.v_HolesMtg);
        if ~isempty(pSrs)&~isempty(pDet)
            if type == 0 | type == 1  %Ne donne pas de bon resultats avec les normals sur le cortex remplacé par la projection de
            %par une projection sur la peau de l'activation et ensuite
            %transfert sur le cortex.
                z_srs = vHoles(pSrs).Coord.x-vHoles(pSrs).Normal.x*vHoles(pSrs).SkinDepth;
                y_srs = vHoles(pSrs).Coord.y-vHoles(pSrs).Normal.y*vHoles(pSrs).SkinDepth;
                x_srs = vHoles(pSrs).Coord.z-vHoles(pSrs).Normal.z*vHoles(pSrs).SkinDepth;
                z_det = vHoles(pDet).Coord.x-vHoles(pDet).Normal.x*vHoles(pDet).SkinDepth;
                y_det = vHoles(pDet).Coord.y-vHoles(pDet).Normal.y*vHoles(pDet).SkinDepth;
                x_det = vHoles(pDet).Coord.z-vHoles(pDet).Normal.z*vHoles(pDet).SkinDepth;
            end
            [theta_srs,phi_srs,r_srs]= cart2sph(x_srs,y_srs,z_srs);
            [theta_det,phi_det,r_det]= cart2sph(x_det,y_det,z_det);
            if theta_srs * theta_det < 0
                if theta_srs < -pi/2
                    theta_srs = theta_srs + 2*pi;
                elseif theta_det < -pi/2
                    theta_det = theta_det + 2*pi;
                end
            end
            if phi_srs * phi_det < 0
                if  phi_srs < -pi/2
                    phi_srs = phi_srs + 2*pi;
                elseif phi_det < -pi/2
                    phi_det = phi_det + 2*pi;
                end
            end
            delta_theta = theta_det-theta_srs;
            delta_phi = phi_det-phi_srs;
            %Définition de la zone à tracer en coordonnée polaire
            d_theta = theta_srs:(delta_theta/7):theta_det ;
            d_phi =  [phi_srs:delta_phi/7:phi_det];
            %Attribution de la couleur au voxel dans la zone du canal
            pzone = [];
            for i = 1:(numel(d_theta)-1)
                x_min = min([d_theta(i),d_theta(i+1)])-.03;
                x_max = max([d_theta(i),d_theta(i+1)])+.03;
                y_min = min([d_phi(i),d_phi(i+1)])-.03;
                y_max = max([d_phi(i),d_phi(i+1)])+.03;
                p_in = find(Vertex_sph(:,1) > x_min & Vertex_sph(:,1) < x_max & Vertex_sph(:,2) > y_min & Vertex_sph(:,2) < y_max);
                pzone = [pzone; p_in];
                if ~isempty(p_in)
                    for p_vertex = 1 : numel(p_in)
                        if ~(Vcolor(p_in(p_vertex),1) == 0)
                            if meanchannel == 0
                                Vcolor(p_in(p_vertex),1) = (Vcolor(p_in(p_vertex),1) + d1(ch))/2; %moyennage des canaux superposés
                            elseif 1
                                val = [Vcolor(p_in(p_vertex),1),d1(ch)];
                                [maxval,id]=max(abs(val));
                                Vcolor(p_in(p_vertex),1) = val(id);            %la plus grande valeur est selectionné minimal ou maximal
                            else
                                Vcolor(p_in(p_vertex),1) = (Vcolor(p_in(p_vertex),1) + d1(ch)); %Somme des canaux superposés
                            end
                        else
                            Vcolor(p_in(p_vertex),1) = d1(ch);
                        end
                    end
                end
            end
        end
    end
end
end
toc
if type == 0
    oMRI = set_SkinVcolor(oMRI,Vcolor);
elseif type == 1
    [VertexBuffer_s, IndexBuffer_s] = get_SkinMesh( oMRI );
    [VertexBuffer_c, IndexBuffer_c] = get_CortexMeshLowRes( oMRI );
    %pour l'utilisation des coordonnées sphériques on déplace le centre du
    %cerveau au centre de la tête, ainsi on peut utiliser la variation d'angle
    %de sphérique pour passer du référentiel de la peau au référentiel du
    %cerveau.
    VertexBuffer_sc(:,1) = VertexBuffer_s(:,1)-max(VertexBuffer_s(:,1))/2;
    VertexBuffer_sc(:,2) = VertexBuffer_s(:,2)-max(VertexBuffer_s(:,2))/2+30;
    VertexBuffer_sc(:,3) = VertexBuffer_s(:,3)-max(VertexBuffer_s(:,3))/2;
    %Même transformation vertex du cortex Replacer au centre avec les mêmes translation que la peau
    VertexBuffer_cc(:,1) = VertexBuffer_c(:,1)-max(VertexBuffer_s(:,1))/2;     %gauche droite
    VertexBuffer_cc(:,2) = VertexBuffer_c(:,2)-max(VertexBuffer_s(:,2))/2+30;  %haut bas
    VertexBuffer_cc(:,3) = VertexBuffer_c(:,3)-max(VertexBuffer_s(:,3))/2;     %avant arrière
    clear VertexBuffer_c;
    clear VertexBuffer_s;
    [Vertex_sph_s(:,1),Vertex_sph_s(:,2),Vertex_sph_s(:,3)] = cart2sph(VertexBuffer_sc(:,1),VertexBuffer_sc(:,2),VertexBuffer_sc(:,3));
    [Vertex_sph_c(:,1),Vertex_sph_c(:,2),Vertex_sph_c(:,3)] = cart2sph(VertexBuffer_cc(:,1),VertexBuffer_cc(:,2),VertexBuffer_cc(:,3));

    %L'arrière est moins bien couvert changement de plan pour eviter l'effet
    %des coordonnées sphérique ou phi est moins bien représenté dans l'axe z
    [Vertex_sph_s2(:,1),Vertex_sph_s2(:,2),Vertex_sph_s2(:,3)] = cart2sph(VertexBuffer_sc(:,2),VertexBuffer_sc(:,3),VertexBuffer_sc(:,1));
    [Vertex_sph_c2(:,1),Vertex_sph_c2(:,2),Vertex_sph_c2(:,3)] = cart2sph(VertexBuffer_cc(:,2),VertexBuffer_cc(:,3),VertexBuffer_cc(:,1));

    [Vertex_sph_s3(:,1),Vertex_sph_s3(:,2),Vertex_sph_s3(:,3)] = cart2sph(VertexBuffer_sc(:,3),VertexBuffer_sc(:,1),VertexBuffer_sc(:,2));
    [Vertex_sph_c3(:,1),Vertex_sph_c3(:,2),Vertex_sph_c3(:,3)] = cart2sph(VertexBuffer_cc(:,3),VertexBuffer_cc(:,1),VertexBuffer_cc(:,2));
    clear VertexBuffer_sc;
    clear VertexBuffer_cc;
    vColor_c = zeros(size(Vertex_sph_c,1),1);
    %Projection
    ind = find(Vcolor);
    [val,indordre] = sort(Vcolor(ind));
    deltaangle = 0.04;
    for i = 1:numel(ind)
        i/numel(ind);
        theta = Vertex_sph_s(ind(indordre(i)),1); %centre de l'activation topo
        phi = Vertex_sph_s(ind(indordre(i)),2);
        theta2 = Vertex_sph_s2(ind(indordre(i)),1);
        phi2 = Vertex_sph_s2(ind(indordre(i)),2);
        theta3 = Vertex_sph_s3(ind(indordre(i)),1);
        phi3 = Vertex_sph_s3(ind(indordre(i)),2);
        indcolor =  find(theta+deltaangle > Vertex_sph_c(:,1)&...
            theta-deltaangle < Vertex_sph_c(:,1)&...
            phi+deltaangle > Vertex_sph_c(:,2)&...
            phi-deltaangle < Vertex_sph_c(:,2));
        indcolor2 =  find(theta2+deltaangle > Vertex_sph_c2(:,1)&...
            theta2-deltaangle < Vertex_sph_c2(:,1)&...
            phi2+deltaangle > Vertex_sph_c2(:,2)&...
            phi2-deltaangle < Vertex_sph_c2(:,2));
        indcolor3 =  find(theta3+deltaangle > Vertex_sph_c3(:,1)&...
            theta3-deltaangle < Vertex_sph_c3(:,1)&...
            phi3+deltaangle > Vertex_sph_c3(:,2)&...
            phi3-deltaangle < Vertex_sph_c3(:,2));
        both=[indcolor;indcolor2;indcolor3];
        vColor_c(both) = Vcolor(ind(indordre(i)));
    end
    oMRI = set_CortexLowResVcolor(oMRI,vColor_c);
end
PrjStruct = set_MRI_Data(PrjStruct, oMRI);