function IO_DisplayHelmet(PrjStruct,PMI,DispParameter)
%**************************************************************************
% function displayHelmet(helmet)
% Fonction permettant gerant l'affichage de tous les elements du casque
% Entr�e : 
%          Prjstrut      structure de donnee IOmtg
%          PMI structure structure de donn�e d'homer 
%          DispParameter parametre d'affichage     
%               DispParameter.axes1 = handles.axes_projection;
%               DispParameter.viewhelmet = 1;
%               DispParameter.reset = 1;  
%               DispParameter.channel = 0;
%               DispParameter.viewfront = 0;
%               DispParameter.dist = 0;
%               DispParameter.idnb = 0;
%               DispParameter.fiducies = 0;
%               DispParameter.viewselected = 1;
%**************************************************************************
   DispHelm = get_Helmet(PrjStruct);
   oMRI = get_MRI_Data(PrjStruct);
   if ~isfield(DispParameter,'viewhelmet')
        DispParameter.viewhelmet = 1;
   end
   if ~isfield(DispParameter,'reset')
       DispParameter.reset = 1;
   end
   if ~isfield(DispParameter,'channel')
       DispParameter.channel = 0;
   end
   if ~isfield(DispParameter,'viewfront')
       DispParameter.viewfront = 0;
   end
   if ~isfield(DispParameter,'viewskin')
       DispParameter.viewskin = 0;
   end
   if ~isfield(DispParameter,'dist')
       DispParameter.dist = 0;
   end
   if ~isfield(DispParameter,'idnb')
    DispParameter.idnb = 0;
   end
   if ~isfield(DispParameter,'fiducies')
       DispParameter.fiducies = 0;
   end
   if ~isfield(DispParameter, 'viewselected')
       DispParameter.viewselected = 0;
   end
   if ~isfield(DispParameter, 'viewskin')
       DispParameter.viewskin = 0;
   end
   if ~isfield(DispParameter, 'viewcortex')
       DispParameter.viewcortex = 0;
   end
   if ~isfield(DispParameter,'MRIviewlight')
       DispParameter.MRIviewlight = 0;
   end
   if ~isfield(DispParameter,'MRIviewtransparent')
       DispParameter.MRIviewtransparent = 0;
   end
   if ~isfield(DispParameter, 'viewMeasListAct')
       DispParameter.viewMeasListAct = 1;
   end
   if ~isfield(DispParameter, 'viewcortexandatlas')
       DispParameter.viewcortexandatlas = 0;
   end
   if ~isfield(DispParameter,'viewd1chcolor') 
        DispParameter.viewd1chcolor = 0;
   end
   if ~isfield(DispParameter,'nbavg');
       DispParameter.nbavg = 0;
   end 
   axes(DispParameter.axes1)
   if DispParameter.reset
        cla
   end
   if ~isfield(DispParameter,'viewscale')
        DispParameter.viewscale = 0;
   end
   if ~isfield(DispParameter,'HideNotCover')
        DispParameter.HideNotCover = 0;
   end
   if ~isfield(DispParameter,'mapoption') %Jet or hot color map
        DispParameter.mapoption = 1; %1 JET, 2 HOT
   end
   if ~isfield(DispParameter,'backgroundcolor')
        DispParameter.backgroundcolor = 'white';
   end
   set(DispParameter.axes1,'xcolor','white');   
   set(DispParameter.axes1,'ycolor','white');   
   set(DispParameter.axes1,'zcolor','white');
   set(DispParameter.axes1,'DataAspectRatio',[1 1 1]);
   hold on
    %Affichage des trous et capteurs
    if DispParameter.viewhelmet & DispParameter.reset
        CameraPosition = get( DispParameter.axes1, 'CameraPosition' );
        display_holes(DispHelm, CameraPosition,DispParameter);
    end
    if DispParameter.viewcortexintersection & DispParameter.reset
        CameraPosition = get( DispParameter.axes1, 'CameraPosition' );
        display_holescortex(DispHelm, CameraPosition,DispParameter);
    end
    %Affichage des canaux s�lectionn�es
    if  DispParameter.channel
        displaychannel(DispHelm,PMI,DispParameter);
    end
%     %temp LOAD 10-20 ref
%     load('C:\data\Data_NIRS\Validation_R_FV_H\AnalyserNEUTRAL\ZoneCreation\DoubleBananaBrainstorm.mat')
%     for ieeg =1:numel(Channel)
%         plot3(Channel(1,ieeg).Loc(1),Channel(1,ieeg).Loc(2),Channel(1,ieeg).Loc(3),'x','displayname',Channel(1,ieeg).Name)
%         text(Channel(1,ieeg).Loc(1),Channel(1,ieeg).Loc(2),Channel(1,ieeg).Loc(3),Channel(1,ieeg).Name)
%     end

    
    %Affichage des canaux s�lectionn�es
    if 0%DispParameter.viewcortexintersection & DispParameter.reset
        displaychannelcortex(DispHelm,PMI,DispParameter);
    end
    %Affichage des points de fiducie
    if DispParameter.fiducies & DispParameter.reset
        displayfiducies(DispHelm,DispParameter);
        displayreference();
    end 
    if DispParameter.viewskin & DispParameter.reset | DispParameter.viewcortex & DispParameter.reset | DispParameter.viewatlas & DispParameter.reset
        displayMRI(oMRI,DispParameter);
    end

%     up = get(DispParameter.axes1,'cameraupvector')
%     Ry = makehgtform('yrotate',pi/2);
%     up(4) = 1
%     up2 = Ry*up'
%     set(DispParameter.axes1,'cameraupvector',up2(1:3))
end

function display_holes(DispHelm, CameraPosition,DispParameter)
%**************************************************************************
% function display_holes(helmet)
% Fonction permettant d'afficher les trous du casque 
%**************************************************************************
 % Affichage des trous du montage         
    ind_s = 1;  
    ind_d = 1;  
    sMtg = get_Mtg( DispHelm );
    vHoles = get_vHoles( DispHelm );
DispParameter.viewd1chcolor = 1;
 if ~isfield(sMtg,'v_HolesMtg_Selected') 
    sMtg.v_HolesMtg_Selected = ones(1,numel(vHoles));
 elseif isempty(sMtg.v_HolesMtg_Selected)
    sMtg.v_HolesMtg_Selected = ones(1,numel(vHoles));
 end
  for( p = 1:numel(sMtg.v_HolesMtg_Selected));
      if ~(DispParameter.viewselected == 1  & ~sMtg.v_HolesMtg_Selected(p) == 1) 
            %affichage des trous de face  
            if find(sMtg.v_pDet==p)  
                rayon = 0.0057; %Changer taille des detectors
            else
                rayon = 0.0027; %Changer taille des sources
            end
            
            t = 0:pi/6:2*pi;
            vn = [ vHoles(p).Normal.x; vHoles(p).Normal.y; vHoles(p).Normal.z ];
            vz = [ 0; 0; 1; ];
            % Produit vectoriel de n et z
            A = vn(2)*vz(3)-vn(3)*vz(2);
            B = vn(3)*vz(1)-vn(1)*vz(3);
            C = vn(1)*vz(2)-vn(2)*vz(1);                  
            tetha = -2*atan2( ((vn-vz)'*(vn-vz))^0.5, ((vn+vz)'*(vn+vz))^0.5 );
            
            %Matrice de rotation autour du vecteur v=(A,B,C)
            matRot = [ (1-A^2)*cos(tetha)+A^2,          A*B*(1-cos(tetha))-C*sin(tetha), A*C*(1-cos(tetha))+B*sin(tetha); ...
                       A*B*(1-cos(tetha))+C*sin(tetha), (1-B^2)*cos(tetha)+B^2,          B*C*(1-cos(tetha))-A*sin(tetha); ...
                       A*C*(1-cos(tetha))-B*sin(tetha), B*C*(1-cos(tetha))+A*sin(tetha), (1-C^2)*cos(tetha)+C^2 ];
        
                vCercle = (matRot*[ rayon*cos(t); rayon*sin(t); 0*t ])';
                
                 if ( CameraPosition*vn > 0  || ~DispParameter.viewfront )%affiche seulement le plan de face
                        if  DispParameter.viewd1chcolor==2
                            posx = vHoles(p).Coord.x;
                            posy = vHoles(p).Coord.y;
                            posz = vHoles(p).Coord.z;
                        else
                            posx = vHoles(p).Coord.x;
                            posy = vHoles(p).Coord.y;
                            posz = vHoles(p).Coord.z;
%                             if posx > 0
%                                 posx = posx+0.003;
%                             else
%                                 posx = posx-0.003;
%                             end
%                              if posy > 0
%                                 posy = posy+0.003;
%                             else
%                                 posy = posy-0.003;
%                             end
%                             if posz > 0
%                                 posz = posz-0.003;
%                             else
%                                 posz = posz+0.003;
%                             end                            
                        plot3( posx+vCercle(:,1), ...
                                              posy+vCercle(:,2), ...
                                             posz+vCercle(:,3));
                        end
                 end
            %Affichage des sources 
            %if  p == sMtg.v_pSrc(ind_s);
             if find(sMtg.v_pSrc==p)
                if ind_s < length(sMtg.v_pSrc)
                    ind_s = ind_s + 1;
                end  
                 [VertexBuffer, IndexBuffer] = CreateSphereMesh( 0.005, 10 );
                
                                         
                if ( CameraPosition*vn > 0  || ~DispParameter.viewfront )%affiche seulement le plan de face
                   % fill3( posx+vCercle(:,1),  posy+vCercle(:,2),  posz+vCercle(:,3), [204/255, 204/255 204/255]);   
                    
                        v_hPrevDispItems(p) = patch( 'Vertices', [VertexBuffer(:,1)+posx,...
                                                          VertexBuffer(:,2)+ posy, ...
                                                          VertexBuffer(:,3)+ posz], ...
                                             'Faces', IndexBuffer, ...
                                             'FaceColor', [204/255 204/255 204/255], ...
                                             'EdgeColor', 'none', ...
                                             'FaceLighting', 'gouraud',... 
                                             'FaceAlpha',0.5);
                    strFiber = get_HoleFiberID( DispHelm , p );
                    if ~isfield(sMtg.Gen_Params,'AcqSystem')
                        sMtg.Gen_Params.AcqSystem = 'ISS';
                    end
                    if strcmp(sMtg.Gen_Params.AcqSystem,'ISS')
                        strFiber = strFiber(1:ceil(end/2));
                    end
                    if DispParameter.viewd1chcolor==2 
                    
                    else
                    text(   posx+vHoles(p).Normal.x*0.0025, ...
                                                       posy+vHoles(p).Normal.y*0.0025, ...
                                                        posz+vHoles(p).Normal.z*0.0025, ...
                                                       strFiber, 'FontSize',14, 'Color','black')%,   ...'VerticalAlignment', 'top');
                    end
                end
             end
            %Affichage des d�tecteurs     
            %if  p == sMtg.v_pDet(ind_d);
            if find(sMtg.v_pDet==p)  
                if ind_d < length(sMtg.v_pDet)
                    ind_d = ind_d + 1;
                end
                  [VertexBuffer, IndexBuffer] = CreateSphereMesh( 0.005, 5 );
                if ( CameraPosition*vn > 0  || ~DispParameter.viewfront )%affiche seulement le plan de face
                   % fill3( posx+vCercle(:,1),posy+vCercle(:,2),posz+vCercle(:,3), [255/255,128/255,255/255]);      
                     v_hPrevDispItems(p) = patch( 'Vertices', [VertexBuffer(:,1)+posx,...
                                                          VertexBuffer(:,2)+ posy, ...
                                                          VertexBuffer(:,3)+ posz], ...
                                             'Faces', IndexBuffer, ...
                                             'FaceColor', [204/255 204/255 204/255], ...
                                             'EdgeColor', 'none', ...
                                             'FaceLighting', 'gouraud',... 
                                             'FaceAlpha',0.5);
                    strFiber = get_HoleFiberID( DispHelm , p );
                    if DispParameter.viewd1chcolor==2 
                        
                    else
                        text(  posx+vHoles(p).Normal.x*0.0025, ...
                                                       posy+vHoles(p).Normal.y*0.0025, ...
                                                       posz+vHoles(p).Normal.z*0.0025, ...
                                                       strFiber, 'FontSize',14, 'Color','black');%,   ...'VerticalAlignment', 'top');    
                    end
                end
            end   
            if find(sMtg.v_pEle==p)  
                if ind_d < length(sMtg.v_pDet)
                    ind_d = ind_d + 1;
                end
                 [VertexBuffer, IndexBuffer] = CreateSphereMesh( 0.005, 5 );
                if ( CameraPosition*vn > 0  || ~DispParameter.viewfront )%affiche seulement le plan de face
                    fill3( vHoles(p).Coord.x+vCercle(:,1), vHoles(p).Coord.y+vCercle(:,2), vHoles(p).Coord.z+vCercle(:,3), [255/255,128/255,255/255]);      
                     v_hPrevDispItems(p) = patch( 'Vertices', [VertexBuffer(:,1)+posx,...
                                                          VertexBuffer(:,2)+ posy, ...
                                                          VertexBuffer(:,3)+ posz], ...
                                             'Faces', IndexBuffer, ...
                                             'FaceColor', [204/255 204/255 204/255], ...
                                             'EdgeColor', 'none', ...
                                             'FaceLighting', 'gouraud',... 
                                             'FaceAlpha',0.5);
                    strFiber = sMtg.v_HolesEle(p);                    
                    text(  posx+vHoles(p).Normal.x*0.0025, ...
                                                       posy+vHoles(p).Normal.y*0.0025, ...
                                                       posz+vHoles(p).Normal.z*0.0025, ...
                                                       strFiber, 'FontSize',14, 'Color','r');%,   ...'VerticalAlignment', 'top');     
                end
            end   
        end 
  end   

end 

function display_holescortex(DispHelm, CameraPosition,DispParameter)
%**************************************************************************
% function display_holes(helmet)
% Fonction permettant d'afficher les trous du casque 
%**************************************************************************
 % Affichage des trous du montage         
    ind_s = 1;  
    ind_d = 1;  
    sMtg = get_Mtg( DispHelm );
    vHoles = get_vHoles( DispHelm );

 if ~isfield(sMtg,'v_HolesMtg_Selected') 
    sMtg.v_HolesMtg_Selected = ones(1,numel(vHoles));
 elseif isempty(sMtg.v_HolesMtg_Selected)
    sMtg.v_HolesMtg_Selected = ones(1,numel(vHoles));
 end
depth_r = 1;

[VertexBuffer, IndexBuffer] = CreateSphereMesh( 0.05/10, 16 );
Mtg_Selected = [sMtg.v_pDet, sMtg.v_pSrc];
  for( i = 1:numel(Mtg_Selected));
      p = Mtg_Selected(i);
      if find(p==sMtg.v_pDet)
          vColor = [0,1,0]; %Green
          vColor = [0.3,0.3,0.3]; %Color holecortex Olivia
          [VertexBuffer, IndexBuffer] = CreateSphereMesh( 0.05/5, 16 );
      else
         % vColor = [1,1,0];%Yellow
         vColor = [1,1,1];
         [VertexBuffer, IndexBuffer] = CreateSphereMesh( 0.05/10, 16 );
      end
      
      if ~(DispParameter.viewselected == 1  & ~sMtg.v_HolesMtg_Selected(p) == 1) 
          
%        load('D:\src_det_NIRS_SPM_2_HSJ.mat');
%         if i <= numel(sMtg.v_pDet)
%             x = nos_detecteurs(1,i)/100;
%             y = nos_detecteurs(2,i)/100;
%             z = nos_detecteurs(3,i)/100;    
%         else 
%             d = i-numel(sMtg.v_pDet)
%             x = nos_sources(1,d)/100;
%             y = nos_sources(2,d)/100;
%             z = nos_sources(3,d)/100;
%         end

       
       x= vHoles(p).Coord.x-vHoles(p).Normal.x*vHoles(p).CortexDepth*depth_r;
       y= vHoles(p).Coord.y-vHoles(p).Normal.y*vHoles(p).CortexDepth*depth_r;
       z= vHoles(p).Coord.z-vHoles(p).Normal.z*vHoles(p).CortexDepth*depth_r;
       x= vHoles(p).Coord.x; %-vHoles(p).Normal.x*vHoles(p).CortexDepth*depth_r;
       y= vHoles(p).Coord.y; %-vHoles(p).Normal.y*vHoles(p).CortexDepth*depth_r;
       z= vHoles(p).Coord.z; %-vHoles(p).Normal.z*vHoles(p).CortexDepth*depth_r;
              v_hPrevDispItems(p) = patch( 'Vertices',  [VertexBuffer(:,1)+x,...
                                                         VertexBuffer(:,2)+y, ...
                                                         VertexBuffer(:,3)+z], ...
                                             'Faces', IndexBuffer, ...
                                             'FaceColor', vColor, ...
                                             'EdgeColor', 'none', ...
                                             'FaceLighting', 'gouraud' );
%             %affichage des trous de face  
%             CoordsH(p).x = vHoles(p).Coord.x-vHoles(p).Normal.x*vHoles(p).CortexDepth*depth_r;
%             CoordsH(p).y = vHoles(p).Coord.y-vHoles(p).Normal.y*vHoles(p).CortexDepth*depth_r;
%             CoordsH(p).z = vHoles(p).Coord.z-vHoles(p).Normal.z*vHoles(p).CortexDepth*depth_r;
%             if vHoles(p).CortexDepth ~= 0                
%                 rayon = 0.0027;
%                 t = 0:pi/6:2*pi;
%                 vn = [ vHoles(p).Normal.x; vHoles(p).Normal.y; vHoles(p).Normal.z ];
%                 vz = [ 0; 0; 1; ];
%                 % Produit vectoriel de n et z
%                 A = vn(2)*vz(3)-vn(3)*vz(2);
%                 B = vn(3)*vz(1)-vn(1)*vz(3);
%                 C = vn(1)*vz(2)-vn(2)*vz(1);                  
%                 tetha = -2*atan2( ((vn-vz)'*(vn-vz))^0.5, ((vn+vz)'*(vn+vz))^0.5 );
% 
%                 %Matrice de rotation autour du vecteur v=(A,B,C)
%                 matRot = [ (1-A^2)*cos(tetha)+A^2,          A*B*(1-cos(tetha))-C*sin(tetha), A*C*(1-cos(tetha))+B*sin(tetha); ...
%                            A*B*(1-cos(tetha))+C*sin(tetha), (1-B^2)*cos(tetha)+B^2,          B*C*(1-cos(tetha))-A*sin(tetha); ...
%                            A*C*(1-cos(tetha))-B*sin(tetha), B*C*(1-cos(tetha))+A*sin(tetha), (1-C^2)*cos(tetha)+C^2 ];
% 
%                     vCercle = (matRot*[ rayon*cos(t); rayon*sin(t); 0*t ])';
%                      if ( CameraPosition*vn > 0  || ~DispParameter.viewfront )%affiche seulement le plan de face
%                              plot3(  CoordsH(p).x+vCercle(:,1), ...
%                                                  CoordsH(p).y+vCercle(:,2), ...
%                                                  CoordsH(p).z+vCercle(:,3));
%                      end
%                 %Affichage des sources 
%                 %if  p == sMtg.v_pSrc(ind_s);
%                  if find(sMtg.v_pSrc==p)
%                     if ind_s < length(sMtg.v_pSrc)
%                         ind_s = ind_s + 1;
%                      end  
%                     if ( CameraPosition*vn > 0  || ~DispParameter.viewfront )%affiche seulement le plan de face
%                         fill3( CoordsH(p).x+vCercle(:,1), CoordsH(p).y+vCercle(:,2), CoordsH(p).z+vCercle(:,3), 'r')      
%                         strFiber = get_HoleFiberID( DispHelm , p );
%                         text(   CoordsH(p).x + vHoles(p).Normal.x*0.0025, ...
%                                                             CoordsH(p).y+vHoles(p).Normal.y*0.0025, ...
%                                                             CoordsH(p).z+vHoles(p).Normal.z*0.0025, ...
%                                                            strFiber, 'FontSize',8, 'Color','black')%,   ...'VerticalAlignment', 'top');
%                     end
%                  end
%                 %Affichage des d�tecteurs     
%                 %if  p == sMtg.v_pDet(ind_d);
%                 if find(sMtg.v_pDet==p)  
%                     if ind_d < length(sMtg.v_pDet)
%                         ind_d = ind_d + 1;
%                     end
%                     if ( CameraPosition*vn > 0  || ~DispParameter.viewfront )%affiche seulement le plan de face
%                         fill3( CoordsH(p).x+vCercle(:,1), CoordsH(p).y+vCercle(:,2), CoordsH(p).z+vCercle(:,3), 'g');      
%                         strFiber = get_HoleFiberID( DispHelm , p );
%                         text(  CoordsH(p).x+vHoles(p).Normal.x*0.0025, ...
%                                                            CoordsH(p).y+vHoles(p).Normal.y*0.0025, ...
%                                                            CoordsH(p).z+vHoles(p).Normal.z*0.0025, ...
%                                                            strFiber, 'FontSize',12, 'Color','black');%,   ...'VerticalAlignment', 'top');     
%                     end
%                 end   
%             end
%         end 
      end
  end   
end 

function displaychannelcortex(DispHelm,PMI,DispParameter)
%**************************************************************************
% function displaychannel
% Fonction permettant d'afficher les canaux sur le casque
%**************************************************************************
  global currentsub;
    sMtg = get_Mtg( DispHelm );
    vHoles = get_vHoles( DispHelm );
    cf = PMI{currentsub}.currentFile;
  for ind = 1:length(PMI{currentsub}.plotLst);
   
    srs = PMI{currentsub}.plot(ind,1);
    det = PMI{currentsub}.plot(ind,2);     
    if (srs ~= 0 && det ~= 0 )     
        Srs_n = SDpairs2Srs_n(srs,1); 
        Det_n = SDpairs2Srs_n(det,2);
        p_srs = find(sMtg.v_HolesMtg == Srs_n);
        flag_display = 1;
        depth_r =1;% 0.80
        
        if ~isempty(p_srs)
            srsx = vHoles(p_srs).Coord.x-vHoles(p_srs).Normal.x*vHoles(p_srs).CortexDepth*depth_r;
            srsy = vHoles(p_srs).Coord.y-vHoles(p_srs).Normal.y*vHoles(p_srs).CortexDepth*depth_r;
            srsz =vHoles(p_srs).Coord.z-vHoles(p_srs).Normal.z*vHoles(p_srs).CortexDepth*depth_r;
        else 
            flag_display = 0 ;
        end
         
         p_det = find(sMtg.v_HolesMtg == Det_n);
         if ~isempty(p_det)
            detx =  vHoles(p_det).Coord.x-vHoles(p_det).Normal.x*vHoles(p_det).CortexDepth*depth_r;
            dety =  vHoles(p_det).Coord.y-vHoles(p_det).Normal.y*vHoles(p_det).CortexDepth*depth_r;
            detz = vHoles(p_det).Coord.z-vHoles(p_det).Normal.z*vHoles(p_det).CortexDepth*depth_r;
         else 
            flag_display = 0 ;
         end
      if flag_display
        d = sqrt((detx-srsx)^2+(dety-srsy)^2+ (detz-srsz)^2)*100;
        d_str = sprintf('%1.1f',d); 
        if PMI{currentsub}.data(cf).MeasListAct(PMI{currentsub}.plotLst(ind))         
                idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList,PMI{currentsub}.plotLst(ind),numel(PMI{currentsub}.color)/3);
                plot3( [srsx,detx], [srsy,dety], [srsz,detz],'color',PMI{currentsub}.color(idxc,:),'LineWidth',2);        
        else
             if DispParameter.viewMeasListAct
               idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList,PMI{currentsub}.plotLst(ind),numel(PMI{currentsub}.color)/3);
               plot3( [srsx,detx], [srsy,dety], [srsz,detz],'color',PMI{currentsub}.color(idxc,:),'LineWidth',2,'linestyle','--');
            end
        end
        
%           DispParameter.D1label = 1;
%     DispParameter.dist = 0;
%     DispParameter.idnb = 0;
%     DispParameter.d1 = handles.d1;
    
        if  DispParameter.D1label
            x = ((srsx+detx)/2)+0.005;
            y = (srsy+dety)/2+0.005;
            z = (srsz+detz)/2+0.005;
            text(x,y,z,d_str,'color','black','FontSize',14);
            DispParameter.D1label
        elseif DispParameter.dist
            x = ((srsx+detx)/2)+0.005;
            y = (srsy+dety)/2+0.005;
            z = (srsz+detz)/2+0.005;
            text(x,y,z,d_str,'color','black','FontSize',14);
        elseif DispParameter.idnb
            x = ((srsx+detx)/2)+0.005;
            y = (srsy+dety)/2+0.005;
            z = (srsz+detz)/2+0.005;
            text(x,y,z,num2str(PMI{currentsub}.plotLst(ind)),'color','black','FontSize',14);      
        elseif DispParameter.navg
            x = ((srsx2+detx2)/2)+0.005;
            y = (srsy2+dety2)/2+0.005;
            z = (srsz2+detz2)/2+0.005;                      
            text(x,y,z,num2str(PMI{currentsub}.data(cf).HRF.navg(PMI{currentsub}.plotLst(ind))),'color','black','FontSize',14,'BackgroundColor',[1,1,1]);
        end
    end
    end
  end

end

function displaychannel(DispHelm,PMI,DispParameter)
%**************************************************************************
% function displaychannel
% Fonction permettant d'afficher les canaux sur le casque
%**************************************************************************
  global currentsub;
    sMtg = get_Mtg( DispHelm );
    vHoles = get_vHoles( DispHelm );
    cf = PMI{currentsub}.currentFile;
    colorjet = colormap(jet(101));

    
  for ind = 1:(numel(PMI{currentsub}.plot)/2);  
      if ~size(PMI{currentsub}.plot,2)
          msgbox('Size wrong')
      end
      try
    srs = PMI{currentsub}.plot(ind,1);
    det = PMI{currentsub}.plot(ind,2);
      catch
          1
      end
    if (srs ~= 0 && det ~= 0 )     
        Srs_n = SDpairs2Srs_n(srs,1); 
        Det_n = SDpairs2Srs_n(det,2);
        p_srs = find(sMtg.v_HolesMtg == Srs_n);
        flag_display = 1;
        if ~isempty(p_srs)
            srsx = vHoles(p_srs).Coord.x;
            srsy = vHoles(p_srs).Coord.y;
            srsz = vHoles(p_srs).Coord.z;
        else 
            flag_display = 0 ;
        end
         
         p_det = find(sMtg.v_HolesMtg == Det_n);
         if ~isempty(p_det)
            detx = vHoles(p_det(1)).Coord.x;
            dety = vHoles(p_det(1)).Coord.y;
            detz = vHoles(p_det(1)).Coord.z;
         else 
            flag_display = 0 ;
         end
        if DispParameter.viewd1chcolor==1 & isfield(PMI{currentsub},'coloramp')         
%          ind = find(isnan(PMI{currentsub}.coloramp));
%          PMI{currentsub}.coloramp(ind) = 0.001;      
        cmin = DispParameter.cmin;
        cmax = DispParameter.cmax;
        d1=round(0+(PMI{currentsub}.coloramp-cmin)./(cmax-cmin)*100);
        idtop = find(d1>100);
        idless = find(d1<=0);
        if ~isempty(idtop);d1(idtop)=100; end
        if ~isempty(idless);d1(idless)=1;; end
        colorch =  colorjet(d1,:);      
        elseif   DispParameter.viewd1chcolor==2        
          colorch = zeros(numel(PMI{currentsub}.color),3); %+1
        else
          colorch = PMI{currentsub}.color;
        end
        
      if flag_display
        d = sqrt((detx-srsx)^2+(dety-srsy)^2+ (detz-srsz)^2)*100;
        d_str = sprintf('%1.1f',d); 
        deltadisplay = 0; %0.003
        if srsx<0 | detx<0
            srsx2 = srsx-deltadisplay;
            detx2 = detx-deltadisplay;       
        else
            srsx2 = srsx+deltadisplay;
            detx2 = detx+deltadisplay;         
        end
        if srsy<0 | dety<0
            srsy2 = srsy-deltadisplay;
            dety2 = dety-deltadisplay;       
        else
            srsy2 = srsy+deltadisplay;
            dety2 = dety+deltadisplay;         
        end
        if srsz<0 | detz<0
            srsz2 = srsz-deltadisplay;
            detz2 = detz-deltadisplay;       
        else
            srsz2 = srsz+deltadisplay;
            detz2 = detz+deltadisplay;         
        end
        
        if PMI{currentsub}.data(cf).MeasListAct(PMI{currentsub}.plotLst(ind))           
                idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList,...
                    PMI{currentsub}.plotLst(ind),numel(colorch)/3);                
                plot3( [srsx2,detx2], [srsy2,dety2], [srsz2,detz2],'color',colorch(idxc,:),'LineWidth',DispParameter.LineWidth  );
                
        else
             if DispParameter.viewMeasListAct
                idxc = find_idx_color(PMI{currentsub}.data(cf).MeasList,...
                    PMI{currentsub}.plotLst(ind),numel(colorch)/3);
                plot3( [srsx2,detx2], [srsy2,dety2], [srsz2,detz2],'color',colorch(idxc,:),'LineWidth',DispParameter.LineWidth ,'linestyle','--');
            end
        end
            x = ((srsx+detx)/2)+0.005*((srsx+detx)/2)/abs(((srsx+detx)/2));
            y = (srsy+dety)/2+0.005*((srsy+dety)/2)/abs(((srsy+dety)/2));
            z = (srsz+detz)/2+0.005*((srsz+detz)/2)/abs(((srsz+detz)/2));
            if isnan(z);%support 2d label position
                 vHoles(p_srs).Coord.x;
                 vHoles(p_det).Coord.x;
                 x = srsx + (detx-srsx)/2;
                 y = srsy + (dety-srsy)/2;
                 z = 0.01;  end 
        if  DispParameter.D1label & isfield(DispParameter,'d1')         
            if PMI{currentsub}.plotLst(ind)<=numel(DispParameter.d1)
                val = DispParameter.d1(PMI{currentsub}.plotLst(ind));  
            else 
                val = '';
            end
            text(x,y,z,sprintf('%1.2f',val),'color','black','FontSize',14);
        elseif DispParameter.dist
            text(x,y,z,d_str,'color','black','FontSize',14);
        elseif DispParameter.idnb
            text(x,y,z,num2str(PMI{currentsub}.plotLst(ind)),'color','black','FontSize',14);         
        elseif DispParameter.nbavg&isfield(PMI{currentsub}.data(cf).HRF,'navg')                 
            text(x,y,z,num2str(PMI{currentsub}.data(cf).HRF.navg(PMI{currentsub}.plotLst(ind))),'color','black','FontSize',14,'BackgroundColor',[1,1,1]);
        end
        
    
    end
    end
  end
 
% ListToWrite = [];
% for i=1:(size(PMI{currentsub}.data(1).MeasList,1)/2)
%               srs = PMI{currentsub}.data(1).MeasList(i,1);
%               det = PMI{currentsub}.data(1).MeasList(i,2);        
%               Srs_n = SDpairs2Srs_n(srs,1); 
%               Det_n = SDpairs2Srs_n(det,2);
%               p_det = find(sMtg.v_HolesMtg == Det_n);
%               p_srs = find(sMtg.v_HolesMtg == Srs_n);            
%             srsx(i) = vHoles(p_srs).Coord.x;
%             srsy(i) = vHoles(p_srs).Coord.y;
%             srsz(i) = vHoles(p_srs).Coord.z;
%             detx(i) = vHoles(p_det).Coord.x;
%             dety(i) = vHoles(p_det).Coord.y;
%             detz(i) = vHoles(p_det).Coord.z;
%             middlepoint(i,1) = srsx(i) + (detx(i)-srsx(i))/2;
%             middlepoint(i,2) = srsy(i) + (dety(i)-srsy(i))/2;
%             middlepoint(i,3) = srsz(i) + (detz(i)-srsz(i))/2; 
%             middlepoint(i,4) = sqrt((detx(i)-srsx(i))^2+(dety(i)-srsy(i))^2+ (detz(i)-srsz(i))^2); 
%             plot3(middlepoint(i,1), middlepoint(i,2), middlepoint(i,3),'or')
%         end
% for ind = 1:size(PMI{currentsub}.plotLst,1)
%        srs = PMI{currentsub}.plot(ind,1);
%        det = PMI{currentsub}.plot(ind,2);
%        ListToWrite(ind,1)=PMI{currentsub}.plotLst(ind);
%        indch = PMI{currentsub}.plotLst(ind);
%        thispoint = middlepoint(indch,:);
%        for iclose = 1:size(middlepoint,1)     
%            if iclose == ind
%                 dist(iclose,:)=1;
%            else
%             if middlepoint(iclose,1)>0
%                 dist(iclose,1) = middlepoint(iclose,1) - thispoint(1);
%             else
%                 dist(iclose,1) = middlepoint(iclose,1) + thispoint(1);
%             end
%             dist(iclose,2) = middlepoint(iclose,2) + thispoint(2);
%             if middlepoint(iclose,3)>0
%                  dist(iclose,3) = middlepoint(iclose,3) - thispoint(3);
%             else
%                  dist(iclose,3) = middlepoint(iclose,3) + thispoint(3);
%             end            
%             dist(iclose,4) = middlepoint(iclose,4) - thispoint(4);
%            end
%                
%        end
%        sume = sum(abs(dist'))';
%        [val,indop]=min(sum(abs(dist')));
%        ListToWrite(ind,2)=indop;    
% end
%    save('filecoordonne.txt','ListToWrite','-ascii')
end

 
  
function displayfiducies(DispHelm,DispParameter)
%**************************************************************************
% function displayfiducies
% Fonction permettant d'afficher les points de fiducies
%**************************************************************************
             sMtg = get_Mtg( DispHelm );
            %Points de fiducie (X+; Y+; Y-)  
            matFids = sMtg.matFiducials;      
            plot3( matFids(1,1), matFids(1,2), matFids(1,3), ...
                                     '.', 'Color','b');
            text(  matFids(1,1), matFids(1,2), matFids(1,3), ...
                                     ' NAS','FontSize',18, 'Color','black');
            plot3( matFids(2,1), matFids(2,2), matFids(2,3), ...
                                     '.', 'Color', 'b' );
            text(  matFids(2,1), matFids(2,2), matFids(2,3), ...
                                     ' LPA', 'FontSize',18, 'Color','black');
            plot3( matFids(3,1), matFids(3,2), matFids(3,3), ...
                                     '.', 'Color','b' );
            text(  matFids(3,1), matFids(3,2), matFids(3,3), ...
                                     ' RPA','FontSize',18, 'Color','black');                            
end

function displayreference()
%**************************************************************************
% function displayreference()
% Fonction permettant le systeme d'axe x,y,z
%**************************************************************************
        plot3( [0;0.04], [0;0], [0;0], 'Color', 'black', 'linewidth',2);
        plot3( [0;0], [0;0.04], [0;0], 'Color', 'black', 'linewidth',2 );
        plot3( [0;0], [0;0], [0;0.04], 'Color', 'black','linewidth',2 );
        text( [0.04], [0], [0], 'x', 'Color', 'black', 'FontSize', 10 );
        text( [0], [0.04], [0], 'y', 'Color', 'black', 'FontSize', 10 );
        text( [0], [0], [0.04], 'z', 'Color', 'black','FontSize', 10 );
end

function displayMRI(oMRI,DispParameter)
%Description: affiche le cortex ou la peau du mri

    axes(DispParameter.axes1);
    DispParameter.mapoption = 1; %2 Couleur hot as fmri
    if DispParameter.mapoption == 1  
        map = jet(300);
    elseif DispParameter.mapoption == 2     
        map = hot(300);        % map = autumn(300)
    end
     cmin=-1;
     cmax=1;
  if (DispParameter.viewskin | DispParameter.viewcortex) &  DispParameter.cthresh    
    cmin = DispParameter.cmin;
    cmax = DispParameter.cmax;
    if cmin>=cmax
        msgbox('Please check the color scale, imposible value')
        return
    end
    nb_color = size(map,1);
    delta_map = (cmax-cmin)/nb_color;
    cthresh =  abs(DispParameter.cthresh);
    nb_map_min = floor(abs(cmin+cthresh)/delta_map);  
    nb_map_max = ceil(abs(cmin-cthresh)/delta_map);
    nb_map_list = nb_map_min:nb_map_max;
    nb_map = size(nb_map_list,2);
   scalelist = cmin : delta_map:cmax;
   sum(scalelist<0);
  
   if sum(scalelist<0)==0
        zeroid = 1;
   else
        zeroid =[ sum(scalelist<0),sum(scalelist<0)+1,sum(scalelist<0)+2];
   end
    for i_map = 1:nb_map;
        if DispParameter.viewskin
            map( nb_map_list(i_map),:) = [255/256 255/256 179/256];
            map(zeroid,:)=repmat([255/256 255/256 179/256],numel(zeroid),1) ;
        elseif  DispParameter.viewcortex
            map( nb_map_list(i_map),:) =  [192/256 192/256 192/256];%[88/256 88/256 88/256];%
             map(zeroid,:)=repmat([192/256 192/256 192/256],numel(zeroid),1); 
        end
    end
    if DispParameter.mapoption == 2
       mapend= hot(300-nb_map_list(end));
        map(nb_map_list(end)+1:end,:)= mapend;
       for i = 1:nb_map_list(1)
        if DispParameter.viewskin
            map( i,:) = [255/256 255/256 179/256];
        elseif  DispParameter.viewcortex
            map( i,:) =  [192/256 192/256 192/256];%[192/256 192/256 192/256];%
        end
       end
      if strcmp(DispParameter.backgroundcolor,'white')
          set(gca,'color',[1 1 1]) 
      else
        set(gca,'color',[0 0 0]) %Couleur Background JT
      end
        %set(gca,'color',[1 1 1]) 
    else
      if strcmp(DispParameter.backgroundcolor,'white')
          set(gca,'color',[1 1 1]) 
      else
        set(gca,'color',[0 0 0]) %Couleur Background JT
      end
    end
      if DispParameter.HideNotCover; %Black uncovered region
          if DispParameter.viewskin
              
              map(zeroid,:)=repmat([0.5,0.5,0.5],numel(zeroid),1);
          else
              map(zeroid,:)=repmat([192/256 192/256 192/256],numel(zeroid),1);
          end
%           
%           if DispParameter.viewskin
%               i_map = ceil(nb_map/2)+1;
%               map( nb_map_list(i_map),:) = [0.5,0.5,0.5];
%               i_map = ceil(nb_map/2);
%               map( nb_map_list(i_map),:) = [0.5,0.5,0.5];
%           else
%          	i_map = ceil(nb_map/2)+1;
%             map( nb_map_list(i_map),:) = [192/256 192/256 192/256];%
%             i_map = ceil(nb_map/2);
%             map( nb_map_list(i_map),:) = [192/256 192/256 192/256];%
%           end
      end
end
    if DispParameter.viewskin
    [VertexBuffer, IndexBuffer] = get_SkinMesh( oMRI );
      if isempty(VertexBuffer)
         % msgbox('Sorry skin is''nt available on this project')
          return
      end
    VertexBuffer(:,4) = 1;
    vColor = get_SkinVcolor(oMRI);      
    if ~(size(vColor,1) == size(IndexBuffer,1)|size(vColor,1) == size(VertexBuffer,1))
        vColor = zeros( size(VertexBuffer,1)); 
    end
    colormap(map);
    elseif DispParameter.viewcortex
      [VertexBuffer, IndexBuffer] = get_CortexMeshLowRes( oMRI );  
     if isempty(VertexBuffer)
          %msgbox('Sorry cortex is''nt available on this project')
          return
     end
      VertexBuffer(:,4) = 1;
      vColor = get_CortexLowResVcolor(oMRI);    
      if ~(size(vColor,1) == size(IndexBuffer,1)|size(vColor,1) == size(VertexBuffer,1))
        vColor = zeros( size(VertexBuffer,1));
      end
      colormap(map);
        
    if DispParameter.viewcortexandatlas       
     [VertexBufferatlas, IndexBufferatlas] = get_CortexMeshHiRes( oMRI );  
          if isempty(VertexBufferatlas)
              %msgbox('Sorry cortex is''nt available on this project')
              return
          end 
      VertexBufferatlas(:,4) = 1;
      
      vColoratlas = get_CortexHiResVcolor(oMRI);   
      if 1
          vColorline = find(vColoratlas(IndexBuffer(:,1))~= vColoratlas(IndexBuffer(:,2))|...
          vColoratlas(IndexBuffer(:,1))~=vColoratlas(IndexBuffer(:,3))|...
          vColoratlas(IndexBuffer(:,2))~=vColoratlas(IndexBuffer(:,3)));%==vColor(IndexBuffer(:,3))) 
          vColor(IndexBufferatlas(vColorline,1))=-100;
          vColor(IndexBufferatlas(vColorline,2))=-100;
          vColor(IndexBufferatlas(vColorline,3))=-100;
      end
    end
    elseif DispParameter.viewatlas
          [VertexBuffer, IndexBuffer] = get_CortexMeshHiRes( oMRI );  
          if isempty(VertexBuffer)
             % msgbox('Sorry cortex is''nt available on this project')
              return
          end 
      VertexBuffer(:,4) = 1;
      vColor = get_CortexHiResVcolor(oMRI);   
      nb_atlas = max(vColor);
%       vColorline = find()
%       vColorline = find(vColor(IndexBuffer(:,1))~= vColor(IndexBuffer(:,2))|...
%                 vColor(IndexBuffer(:,1))~=vColor(IndexBuffer(:,3))|...
%                 vColor(IndexBuffer(:,2))~=vColor(IndexBuffer(:,3)));%==vColor(IndexBuffer(:,3))) 
%       vColor(IndexBuffer(vColorline,1))=1000;
%       vColor(IndexBuffer(vColorline,2))=1000;
%       vColor(IndexBuffer(vColorline,3))=1000;
      

      if size(vColor,1) ~= size(VertexBuffer,1)
        vColor = zeros( size(VertexBuffer,1),1);
      end
      nb_color = nb_atlas;
      mapcolorcube = colorcube(300);
      map = mapcolorcube(1:nb_color,:);       
      if ~isempty(DispParameter.enumaltas)
          map = ones(nb_color,3)*(159/256); 
          nb_enumatlas = numel(DispParameter.enumaltas);
           mixcolor = lines(nb_enumatlas); 
%           mixcolor = [255/255,128/255,0;
%                       255/255,128/255,0;
%                       255/255,255/255,0;
%                       255/255,255/255,0;
%                       0,255/255,0;
%                       0,255/255,0;]
%           mixcolor = [0,128/255,0;
%                      0,128/255,0;
%                      0,128/255,0];
          for ind_map = 1:numel(DispParameter.enumaltas)
              map(DispParameter.enumaltas(ind_map),:) = mixcolor(ind_map,:);
          end
          colormap (map);      
      else
          colormap (map); 
      end
      cmin = 1;
      cmax = nb_atlas;
    end
    if DispParameter.viewscale %timecolor
         map = jet(300);
         nb_color = size(map,1);
        delta_map = (cmax-cmin)/nb_color;
        cthresh =  DispParameter.cthresh;
        nb_map_min = nb_color - floor((cmax-cthresh)/delta_map);  
        scalelist = cmin : delta_map:cmax;
        sum(scalelist<0);
         if sum(scalelist<0)==0
            zeroid = 1;
         else
            zeroid =[ sum(scalelist<0),sum(scalelist<0)+1,sum(scalelist<0)+2];
        end
%     nb_map_list = nb_map_min:nb_map_max;
%     nb_map = size(nb_map_list,2);
         if DispParameter.viewskin
            map( 1:nb_map_min,:) = repmat([255/256 255/256 179/256],nb_map_min,1) ;
        elseif  DispParameter.viewcortex
            map(1:nb_map_min,:) = repmat([192/256 192/256 192/256],nb_map_min,1);%[192/256 192/256 192/256];%
         end
        
        if DispParameter.HideNotCover; %Black uncovered region
          if DispParameter.viewskin
%               i_map = ceil(nb_map/2)+1;
%               map( nb_map_list(i_map),:) = [0.5,0.5,0.5];
%               i_map = ceil(nb_map/2);
%               map( nb_map_list(i_map),:) = [0.5,0.5,0.5];
               map(zeroid,:)=repmat([0.5,0.5,0.5],numel(zeroid),1); 
          else
%          	i_map = ceil(nb_map/2)+1;
%             map( nb_map_list(i_map),:) = [192/256 192/256 192/256];%
%             i_map = ceil(nb_map/2);
%             map( nb_map_list(i_map),:) = [192/256 192/256 192/256];%
             map(zeroid,:)=repmat([192/256 192/256 192/256],numel(zeroid),1); 
          end
      end
         
         colormap (map);
    end
    T = get_matTransform(oMRI);
    Vertex_tmp = VertexBuffer*T;
    
   caxis([cmin,cmax])
   c = colorbar;
    set(c,'fontsize',14)
   
% %   %Julie pour dima
%   %vColor atlas
% vColoratlas = vColor
%  idBrL = find((vColor==44 |vColor==45) .*Vertex_tmp(:,2)>0 );
%   idBrR = find((vColor==44 |vColor==45) .* Vertex_tmp(:,2)<0 );
% 
%   idWerL = find((vColor==22 ).*Vertex_tmp(:,2)>0  );
%   idWerR = find((vColor==22) .* Vertex_tmp(:,2)<0 );
%   
%   idFusL = find((vColor==37) .*Vertex_tmp(:,2)>0  );
%   idFusR = find((vColor==37) .* Vertex_tmp(:,2)<0 );
%   
%   idOccL = find((vColor==17 |vColor==18) .*Vertex_tmp(:,2)>0  );
%   idOccR = find((vColor==17 |vColor==18) .* Vertex_tmp(:,2)<0 );
%   
% 
%   vColor(:) = 0
% %Mot
%   val = [0.742995746	0.681963494	0.905367301	0.82797232	0.698903031	0.129505308	0.409583534	0.368797563]
% % %Non mot
%    val = [ 0.873246097	0.679594567	0.935210603	0.984508999	1.414158603	0.140421993	0.721521617	0.347314877]
%  %mot
%    val = [3.79892426948574,1.96488584540338,3.51790492856442,2.51174417989247,3.08049398453432,-0.126780845033441,1.13068755890819,1.01536433647187;]
%  %non mot 
%  val  =[2.95032042479363,4.14087727607235,1.51161718148363,2.28820205564092,2.78626658810554,0.500551368127353,1.27667900132864,0.575250117343282;]
%    vColor(idBrL) = val(1)
%   vColor(idBrR) =  val(2)
%   vColor(idWerL) = val(3)
%   vColor(idWerR) = val(4)
%   vColor(idFusL) = val(5)
%   vColor(idFusR) =  val(6)
%   vColor(idOccL) = val(7)
%   vColor(idOccR) = val(8)
%   %Brocaleft = 
%    map = jet(150);
%    map(1,:) =[88/256 88/256 88/256];
%   colormap(map);
% Broca = [44,45]
% Wernicke = 22
% Occipital = 17,18
% Fusiforme = 39
%	0,681963494	0,905367301	0,82797232	0,698903031	0,129505308	0,409583534	0,368797563

%caxis([0,5])
%finDima

    %IndexBuffer
%     fv.faces = IndexBuffer;
%     fv.vertices = Vertex_tmp(:,1:3);
%     fv.facevertexcdata = vColor;
%     nfv = reducepatch(fv,0.5);
%     IndexBuffer = nfv.faces;
%     vColor = vColor(1:5122,:
%     clear Vertex_tmp
%     Vertex_tmp = nfv.vertices;
    
    
    if DispParameter.MRIviewtransparent        
        v = patch('Vertices',Vertex_tmp(:,1:3),...
                                     'Faces',IndexBuffer, ...
                                     'EdgeColor', 'none', ...
                                     'FaceVertexCData',  vColor, ...
                                     'Facecolor','flat',...
                                     'FaceLighting', 'gouraud', ...
                                     'SpecularColorReflectance', 0.5, ...
                                     'DiffuseStrength', 0.5, ...
                                     'SpecularStrength', 0.5, ...
                                     'AmbientStrength', 0.5,...
                                     'FaceAlpha', 0.5);
    else
        if 1 
         v = patch('Vertices',Vertex_tmp(:,1:3),...
                                     'Faces',IndexBuffer, ...
                                     'EdgeColor', 'none', ...
                                     'FaceVertexCData',  vColor, ...
                                     'Facecolor','interp',...
                                     'FaceLighting', 'gouraud', ...
                                     'SpecularColorReflectance', 0.5, ...
                                     'DiffuseStrength', 0.5, ...
                                     'SpecularStrength', 0.5, ...
                                     'AmbientStrength', 0.5);  

       else           
        v = patch('Vertices',Vertex_tmp(:,1:3),...
                                     'Faces',IndexBuffer, ...
                                     'EdgeColor', 'none', ...
                                     'FaceVertexCData',  vColor, ...
                                     'Facecolor','flat',...
                                     'FaceLighting', 'gouraud', ...
                                     'SpecularColorReflectance', 0.5, ...
                                     'DiffuseStrength', 0.5, ...
                                     'SpecularStrength', 0.5, ...
                                     'AmbientStrength', 0.5);     
        end
    end

    if DispParameter.MRIviewlight
        v_hPrevDispItems(1) = light( 'Position', [500,500,500] );
        v_hPrevDispItems(2) = light( 'Position', [500,-500,500] );
        v_hPrevDispItems(3) = light( 'Position', [-500,0,0] );
    end
   
end