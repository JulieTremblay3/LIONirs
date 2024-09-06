%Description : Utilise les valeurs de concentration a un temps donnee et
%produit un colormap sur le cortex ou la peau. Moyenne les canaux
%superposés. 
function PrjStruct = display_MRIcolor(PrjStruct,PMI,d1,type)
%type 0  = skin
%type 1  = cortex

    global currentsub
    cf = PMI{currentsub}.currentFile;
    oMRI = get_MRI_Data(PrjStruct);
    if type == 0
        [VertexBuffer, IndexBuffer] = get_SkinMesh( oMRI );
        if numel(d1)==size(VertexBuffer,1)
             oMRI = set_SkinVcolor(oMRI,d1);
              PrjStruct = set_MRI_Data(PrjStruct, oMRI);
              return
        end    
    end
    if type == 1
        [VertexBuffer, IndexBuffer] = get_CortexMeshLowRes( oMRI );
         if numel(d1)==size(VertexBuffer,1)
             oMRI = set_CortexLowResVcolor(oMRI,d1);
              PrjStruct = set_MRI_Data(PrjStruct, oMRI);
              return
        end    
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
    [Vertex_sph(:,1),Vertex_sph(:,2),Vertex_sph(:,3)] = cart2sph(Vertex_tmp(:,1),Vertex_tmp(:,2),Vertex_tmp(:,3));
      pzone =[];
    for ch = 1:size(ml,1)/2
        % Coordonnées du canal à tracer
        if PMI{currentsub}.data(cf).MeasListAct(ch)& ml(ch,1)~= 0;
        Srs_n = SDpairs2Srs_n(ml(ch,1),1); 
        Det_n = SDpairs2Srs_n(ml(ch,2),2);
        pSrs = find( Srs_n == sMtg.v_HolesMtg);
        pDet = find( Det_n == sMtg.v_HolesMtg);
        if ~isempty(pSrs)&~isempty(pDet)
%         if isempty(pDet)
%             1
%         end
        if type == 0         
            x_srs = vHoles(pSrs).Coord.x-vHoles(pSrs).Normal.x*vHoles(pSrs).SkinDepth;
            y_srs = vHoles(pSrs).Coord.y-vHoles(pSrs).Normal.y*vHoles(pSrs).SkinDepth;
            z_srs = vHoles(pSrs).Coord.z-vHoles(pSrs).Normal.z*vHoles(pSrs).SkinDepth;
            x_det = vHoles(pDet).Coord.x-vHoles(pDet).Normal.x*vHoles(pDet).SkinDepth;
            y_det = vHoles(pDet).Coord.y-vHoles(pDet).Normal.y*vHoles(pDet).SkinDepth;
            z_det = vHoles(pDet).Coord.z-vHoles(pDet).Normal.z*vHoles(pDet).SkinDepth;
        end
        if type == 1
%             load('D:\src_det_NIRS_SPM_2_HSJ.mat');
%             x_srs = nos_sources(1,ml(ch,1))/100;
%             y_srs = nos_sources(2,ml(ch,1))/100;
%             z_srs = nos_sources(3,ml(ch,1))/100;
%             x_det = nos_detecteurs(1,ml(ch,2))/100;
%             y_det = nos_detecteurs(2,ml(ch,2))/100;
%             z_det = nos_detecteurs(3,ml(ch,2))/100;          
            
            x_srs = vHoles(pSrs).Coord.x-vHoles(pSrs).Normal.x*vHoles(pSrs).CortexDepth;
            y_srs = vHoles(pSrs).Coord.y-vHoles(pSrs).Normal.y*vHoles(pSrs).CortexDepth;
            z_srs = vHoles(pSrs).Coord.z-vHoles(pSrs).Normal.z*vHoles(pSrs).CortexDepth;
            x_det = vHoles(pDet).Coord.x-vHoles(pDet).Normal.x*vHoles(pDet).CortexDepth;
            y_det = vHoles(pDet).Coord.y-vHoles(pDet).Normal.y*vHoles(pDet).CortexDepth;
            z_det = vHoles(pDet).Coord.z-vHoles(pDet).Normal.z*vHoles(pDet).CortexDepth;
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
%         d_theta = theta_srs:(delta_theta/7):theta_det ;
%         d_phi =  [phi_srs:delta_phi/7:phi_det];
        d_theta = theta_srs:(delta_theta/4):theta_det ;
        d_phi =  [phi_srs:delta_phi/4:phi_det];
        %Attribution de la couleur au voxel dans la zone du canal
           for i = 1:(numel(d_theta)-1)
            x_min = min([d_theta(i),d_theta(i+1)])-.06;
            x_max = max([d_theta(i),d_theta(i+1)])+.06;
            y_min = min([d_phi(i),d_phi(i+1)])-.06;
            y_max = max([d_phi(i),d_phi(i+1)])+.06;
            p_in = find(Vertex_sph(:,1) > x_min & Vertex_sph(:,1) < x_max & Vertex_sph(:,2) > y_min & Vertex_sph(:,2) < y_max);
            pzone = [pzone; p_in]; 
            if ~isempty(p_in)
                      for p_vertex = 1 : numel(p_in)
                          if ~(Vcolor(p_in,1) == 0)
                              if meanchannel == 0
                                    Vcolor(p_in(p_vertex),1) = (Vcolor(p_in(p_vertex),1) + d1(ch))/2; %moyennage des canaux superposés
                              elseif 1
                                  val = [Vcolor(p_in(p_vertex),1),d1(ch)];
                                  [maxval,id]=max(abs(val));
                                    Vcolor(p_in(p_vertex),1) = val(id);            %la plus grande valeur est selectionné
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
     if type == 0
        oMRI = set_SkinVcolor(oMRI,Vcolor);
     else type == 1
       oMRI = set_CortexLowResVcolor(oMRI,Vcolor);  
     end
    PrjStruct = set_MRI_Data(PrjStruct, oMRI);