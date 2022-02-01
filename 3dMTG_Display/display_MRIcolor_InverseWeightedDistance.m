%Description : Utilise les valeurs de concentration a un temps donnee et
%produit un colormap sur le cortex ou la peau. Projection pondéré par
%l'inverse de la distance du centre des canaux. 
function PrjStruct = display_MRIcolor_InverseWeightedDistance(PrjStruct,PMI,d1,type)
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
%L'arrière est moins bien couvert changement de plan pour eviter l'effet
%des coordonnées sphérique ou phi est moins bien représenté dans l'axe
%z on effectue un changement de plan pour éviter ce défaut
%CAS1 x,y,z standard
%[Vertex_sph(:,1),Vertex_sph(:,2),Vertex_sph(:,3)] = cart2sph(Vertex_tmp(:,1),Vertex_tmp(:,2),Vertex_tmp(:,3));
for ch = 1:size(ml,1)/2 %Passage 1 coordonné 
    % Coordonnées du canal à tracer
    if PMI{currentsub}.data(cf).MeasListAct(ch)& ml(ch,1)~= 0;
        Srs_n = SDpairs2Srs_n(ml(ch,1),1);
        Det_n = SDpairs2Srs_n(ml(ch,2),2);
        pSrs = find( Srs_n == sMtg.v_HolesMtg);
        pDet = find( Det_n == sMtg.v_HolesMtg); 
        if ~isempty(pSrs)&~isempty(pDet)
            if type == 0 
            %transfert on the skin.
                x_srs = vHoles(pSrs).Coord.x;%-vHoles(pSrs).Normal.x*vHoles(pSrs).SkinDepth;
                y_srs = vHoles(pSrs).Coord.y;%-vHoles(pSrs).Normal.y*vHoles(pSrs).SkinDepth;
                z_srs = vHoles(pSrs).Coord.z;%-vHoles(pSrs).Normal.z*vHoles(pSrs).SkinDepth;
                x_det = vHoles(pDet(1)).Coord.x;%-vHoles(pDet).Normal.x*vHoles(pDet).SkinDepth;
                y_det = vHoles(pDet(1)).Coord.y;%-vHoles(pDet).Normal.y*vHoles(pDet).SkinDepth;
                z_det = vHoles(pDet(1)).Coord.z;%-vHoles(pDet).Normal.z*vHoles(pDet).SkinDepth;
            elseif type == 1 
            %transfert on the cortex
                x_srs = vHoles(pSrs).Coord.x;%-vHoles(pSrs).Normal.x*vHoles(pSrs).CortexDepth;
                y_srs = vHoles(pSrs).Coord.y;%-vHoles(pSrs).Normal.y*vHoles(pSrs).CortexDepth;
                z_srs = vHoles(pSrs).Coord.z;%-vHoles(pSrs).Normal.z*vHoles(pSrs).CortexDepth;
                x_det = vHoles(pDet(1)).Coord.x;%-vHoles(pDet).Normal.x*vHoles(pDet).CortexDepth;
                y_det = vHoles(pDet(1)).Coord.y;%-vHoles(pDet).Normal.y*vHoles(pDet).CortexDepth;
                z_det = vHoles(pDet(1)).Coord.z;%-vHoles(pDet).Normal.z*vHoles(pDet).CortexDepth;   
            else
                x_srs = vHoles(pSrs).Coord.x;%-vHoles(pSrs).Normal.x*vHoles(pSrs).SkinDepth;
                y_srs = vHoles(pSrs).Coord.y;%-vHoles(pSrs).Normal.y*vHoles(pSrs).SkinDepth;
                z_srs = vHoles(pSrs).Coord.z;%-vHoles(pSrs).Normal.z*vHoles(pSrs).SkinDepth;
                x_det = vHoles(pDet(1)).Coord.x;%-vHoles(pDet).Normal.x*vHoles(pDet).SkinDepth;
                y_det = vHoles(pDet(1)).Coord.y;%-vHoles(pDet).Normal.y*vHoles(pDet).SkinDepth;
                z_det = vHoles(pDet(1)).Coord.z;%-vHoles(pDet).Normal.z*vHoles(pDet).SkinDepth;
            end
              x_c = x_srs-(x_srs-x_det)/2;
              y_c = y_srs-(y_srs-y_det)/2;
              z_c = z_srs-(z_srs-z_det)/2;
            for ivertice = 1:size(Vertex_tmp,1)
                 x = Vertex_tmp(ivertice,1);
                 y = Vertex_tmp(ivertice,2);
                 z = Vertex_tmp(ivertice,3);
                 d =sqrt( (x-x_c) * (x-x_c) +  (y-y_c) * (y-y_c) + (z-z_c) * (z-z_c));
                 if d<0.03
                    weight(ivertice,ch) = 1/(d);      
                 else
                    weight(ivertice,ch) = 0;    
                 end
            end
        end
    end
end

sum(isnan(weight));
sum(isnan(d1));

weightN = weight./(ones(size(weight,1),1)*max(weight));
idremove = find(max(weight)==0);
weightN(:,idremove)=0; %set weight =0 if channel rejected
if size(d1,1)==1
    Vcolor = weightN*(d1');
elseif size(d1,2)==1
    Vcolor = weightN*(d1);
end
%Vcolor = weight(:,:);

toc
if type == 0
    oMRI = set_SkinVcolor(oMRI,Vcolor);
elseif type == 1
    oMRI = set_CortexLowResVcolor(oMRI,Vcolor);
end
PrjStruct = set_MRI_Data(PrjStruct, oMRI);