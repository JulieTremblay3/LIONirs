% Fonction permettant de calculer la profondeur des trous du casque qui
% fait en sorte que les optodes arrivent sur la peau (intersection de la
% segmentation et de la peau).
%
%
function oHelmet = Fit_Helmet_OnMRISegmentations_manual( oHelmet, oMRI_Data, bProcessSkin, bProcessCortex, bDisplay )

    vHoles = get_vHoles(oHelmet); 
    [VertexBuffer, IndexBuffer] = get_SkinMesh( oMRI_Data ); 
    VertexBuffer(:,4) = 1;
    T = get_matTransform(oMRI_Data);
    Vertex_tmp = VertexBuffer*T; 

    vHoles = get_vHoles( oHelmet );
    matCoords = zeros(numel(vHoles),3);
    matNormals = zeros(numel(vHoles),3);
      for p = 1:numel(vHoles)
        matCoords(p,:) = [vHoles(p).Coord.x,vHoles(p).Coord.y,vHoles(p).Coord.z ];  
        matNormals(p,:) = [vHoles(p).Normal.x,vHoles(p).Normal.y,vHoles(p).Normal.z];
      end
    
      1
 h = waitbar(0,'Please wait...');
SkinDepthSteps = -0.03:0.002:0.03;%step 5 cm avant et après la normal calculé à interval de 1m pour trouver l'intersection
CortexDepthSteps = -0.03:0.002:0.03;%step 5 cm avant et après la normal calculé à interval de 1m pour trouver l'intersection

 for p = 1:numel(vHoles)
      tic 
    HCoordsStepsX = matCoords(p,1)-matNormals(p,1)*SkinDepthSteps;
    HCoordsStepsY = matCoords(p,2)-matNormals(p,2)*SkinDepthSteps;
    HCoordsStepsZ = matCoords(p,3)-matNormals(p,3)*SkinDepthSteps;     
        
    errorx = abs(Vertex_tmp(:,1)*ones(1,numel(SkinDepthSteps))-ones(size(Vertex_tmp,1),1)*HCoordsStepsX); 
    errory = abs(Vertex_tmp(:,2)*ones(1,numel(SkinDepthSteps))-ones(size(Vertex_tmp,1),1)*HCoordsStepsY);
    errorz = abs(Vertex_tmp(:,3)*ones(1,numel(SkinDepthSteps))-ones(size(Vertex_tmp,1),1)*HCoordsStepsZ);
     error =errorx+errory+errorz;
    [val,ivertex] = (min(sum(error,2)));        
    [val,Idepth]=min(error(ivertex,:));
    vHoles(p).SkinDepth = SkinDepthSteps(Idepth);
    
    waitbar(p / numel(vHoles))

 % vHoles(pElem).CortexDepth = CortexDepth
 
 end
 
 
  
   
 
 close(h) 

%figure;plot(errorx(300,:))



%     for( p = 1:size(HCoords,1) )
%     %vHoles actual skin position
%         Hskin(p,1) = (vHoles(p).Coord.x-vHoles(p).Normal.x*vHoles(p).SkinDepth)
%         Hskin(p,2) = (vHoles(p).Coord.y-vHoles(p).Normal.y*vHoles(p).SkinDepth)
%         Hskin(p,3) = (vHoles(p).Coord.z-vHoles(p).Normal.z*vHoles(p).SkinDepth)
%         Hskin(p,4) = 1
%     end
   %Deplacement matrix new skin position. 
   
%      
%      
%     %Fonction de transfert des coordonnees de trous en voxel-space du MRI
%     [HCoords, HNorms ] = CoRegisterHelmetToMRI( oHelmet, oMRI_Data,1);
% 
%     %Memorisation des matrices de segmentations
%     matSkinSegmentation = get_SkinSegmentation( oMRI_Data );
%     A = zeros(size(matSkinSegmentation,1)+50,size(matSkinSegmentation,2)+50,size(matSkinSegmentation,3)+50) ;
%     A(26:25+size(matSkinSegmentation,1),26:25+size(matSkinSegmentation,2),26:25+size(matSkinSegmentation,3) ) =  matSkinSegmentation;
%     matCortexSegmentation = get_CortexSegmentation( oMRI_Data );
%     [nVoxSkinSegX,nVoxSkinSegY,nVoxSkinSegZ] = size( matSkinSegmentation );
%     [nVoxCortSegX,nVoxCortSegY,nVoxCortSegZ] = size( matSkinSegmentation );
%     
%     
%     
%     vPhysicalVoxDim = get_PhysicalVoxDim( oMRI_Data );
%     if( numel(find( vPhysicalVoxDim == vPhysicalVoxDim(1) ) ) < 3 )
%         warndlg( 'Voxel size is not cubic: possible errors in distances calculated on skin segmentation (ERR204)');
%     end
%     %MaxDepth_m = 0.05; %in meters (5 cm)
%     MaxDepth_m = 0.10;
%     if vPhysicalVoxDim == 0
%         MaxDepth_m = 0.10;
%         MaxDepth_vox = 1;
%     else
%         MaxDepth_vox = ceil(MaxDepth_m/mean(vPhysicalVoxDim));
%     end
%     %Structure de memorisation des profondeurs de trous.
%     %vHoles_SkinDepth = zeros( 1,size(HCoords,1) );
%     if( exist( 'bDisplay' ) && bDisplay )
%         hold on;
%         [VertexBuffer, IndexBuffer] = get_CortexMeshHiRes( oMRI_Data );
%         patch('Vertices',VertexBuffer,'Faces',IndexBuffer, ...
%               'EdgeColor', 'none', 'FaceColor', [0.5,0.5,0.5], ...
%               'FaceLighting', 'gouraud', ...
%               'SpecularColorReflectance', 0.5, ...
%               'DiffuseStrength', 0.5, ...
%               'SpecularStrength', 1);
%         light( 'Position', [500,500,500] );
%         light( 'Position', [500,-500,500] );
%         light( 'Position', [-500,0,0] );
%     end
%         %TESTJULIE
% %   
%         matTransformManual = get_matTransformManual(oMRI_Data)
%         New= kspace3d(A, matTransformManual );        
%         h = gcf
%         figure
%         subplot(3,2,1)
%         imagesc(squeeze(A(:,:,100)))
%         subplot(3,2,2)
%         imagesc(squeeze(New(:,:,100)))     
%         subplot(3,2,3)
%         imagesc(squeeze(A(:,100,:)))
%         subplot(3,2,4)
%         imagesc(squeeze(New(:,100,:)))    
%         subplot(3,2,5)
%         imagesc(squeeze(A(100,:,:)))
%         subplot(3,2,6)
%         imagesc(squeeze(New(100,:,:)))           
%         matSkinSegmentation=New;         
%         figure(h)
%         HCoords(:,1:3) = HCoords(:,1:3) + 25
%         nVoxSkinSegX = nVoxSkinSegX+50;
%         nVoxSkinSegY = nVoxSkinSegY+50;
%         nVoxSkinSegZ = nVoxSkinSegZ+50;
% 
%     %Pour chaque position de trou ce casque (en voxel-space)
%     for( pElem = 1:size(HCoords,1) )
%         if pElem==6
%             1
%         end
%        % DepthSteps = 0:1:MaxDepth_vox;  intial modifier pour sortir de la
%        % tête si digitaliser à l'intérieur de la peau. 
%          DepthSteps = -(MaxDepth_vox/2):1:MaxDepth_vox;
%         HCoordsStepsX = int16(HCoords(pElem,1)-HNorms(pElem,1)*DepthSteps);
%         HCoordsStepsY = int16(HCoords(pElem,2)-HNorms(pElem,2)*DepthSteps);
%         HCoordsStepsZ = int16(HCoords(pElem,3)-HNorms(pElem,3)*DepthSteps);
%         
%         %hold on;
%         if( exist( 'bDisplay' ) && bDisplay )
%             plot3( HCoordsStepsX, HCoordsStepsY,HCoordsStepsZ, 'r' );
%         end
%     
%         
%         %----------------------------------------------------------------------------------------------------------
%         %-- SKIN  -------------------------------------------------------------------------------------------------
%         %----------------------------------------------------------------------------------------------------------
%          if( bProcessSkin && ~isempty(matSkinSegmentation) ) %IN VOXEL JT 
%             SkinReached = false;
%             SkinDepth = 0;
%             %Pour au moins chaque voxel dans la trajectoire de la normale
%             for( iDepth=1:numel(DepthSteps) )
% 
%                 %Si a l'interieur de la limite de voxels (-1 et +1 pour cube)
%                 if(   ( HCoordsStepsX(iDepth) > 1 && HCoordsStepsX(iDepth) < nVoxSkinSegX ) ...
%                     &&( HCoordsStepsY(iDepth) > 1 && HCoordsStepsY(iDepth) < nVoxSkinSegY ) ...
%                     &&( HCoordsStepsZ(iDepth) > 1 && HCoordsStepsZ(iDepth) < nVoxSkinSegZ ) )
% 
%                     %Verifier que le nombre de voxels du cube centre en
%                     %HCoordsStepsX/Y/Z faisant partie de la segmentation de la
%                     %peau est d'au moins 13 (13/27 des voxels du cube)
%                     if( numel( find( matSkinSegmentation(HCoordsStepsX(iDepth)-1:HCoordsStepsX(iDepth)+1, ...
%                                                          HCoordsStepsY(iDepth)-1:HCoordsStepsY(iDepth)+1, ...
%                                                          HCoordsStepsZ(iDepth)-1:HCoordsStepsZ(iDepth)+1 )>20 ) ) >= 13 )
%                         %Conversion de profondeur vox->m (distance euclidienne
%                         %respectee)
%                         %SkinDepth = ((iDepth-1)/MaxDepth_vox)*MaxDepth_m;
%                         SkinDepth = ((DepthSteps(iDepth)-1)/MaxDepth_vox)*MaxDepth_m;
%                         SkinReached = true;
%                         if( exist( 'bDisplay' ) && bDisplay )
%                             plot3( [HCoordsStepsX(iDepth),HCoordsStepsX(iDepth)-1], ...
%                                    [HCoordsStepsY(iDepth),HCoordsStepsY(iDepth)-1], ...
%                                    [HCoordsStepsZ(iDepth),HCoordsStepsZ(iDepth)-1], ...
%                                    'color', [1,1,1], 'linewidth', 3 );
%                         end
%                         break;
%                     end
%                 end
%             end
%             vHoles(pElem).SkinDepth = SkinDepth;
%             if pElem == size(HCoords,1)
%                 vHoles(pElem).SkinDepth
%                 vHoles(pElem).Normal.x
%                 vHoles(pElem).Normal.y
%                 vHoles(pElem).Normal.z
%               
%             end
%             Skin(pElem,1) = SkinDepth
%          end     
%         %----------------------------------------------------------------------------------------------------------
%         %-- CORTEX  -------------------------------------------------------------------------------------------------
%         %----------------------------------------------------------------------------------------------------------
%         if( bProcessCortex && ~isempty(matCortexSegmentation) )
%             CortexReached = false;
%             CortexDepth = 0;
%             %Pour au moins chaque voxel dans la trajectoire de la normale
%             for( iDepth=1:numel(DepthSteps) )
% 
%                 %Si a l'interieur de la limite de voxels (-1 et +1 pour cube)
%                 if(   ( HCoordsStepsX(iDepth) > 1 && HCoordsStepsX(iDepth) < nVoxCortSegX ) ...
%                     &&( HCoordsStepsY(iDepth) > 1 && HCoordsStepsY(iDepth) < nVoxCortSegY ) ...
%                     &&( HCoordsStepsZ(iDepth) > 1 && HCoordsStepsZ(iDepth) < nVoxCortSegZ ) )
% 
%                     %Verifier que le nombre de voxels du cube centre en
%                     %HCoordsStepsX/Y/Z faisant partie de la segmentation de la
%                     %peau est d'au moins 13 (13/27 des voxels du cube)
%                     if( numel( find( matCortexSegmentation(HCoordsStepsX(iDepth)-1:HCoordsStepsX(iDepth)+1, ...
%                                                            HCoordsStepsY(iDepth)-1:HCoordsStepsY(iDepth)+1, ...
%                                                            HCoordsStepsZ(iDepth)-1:HCoordsStepsZ(iDepth)+1 ) ) ) >= 13 )
%                         %Conversion de profondeur vox->m (distance euclidienne respectee)
%                         CortexDepth = ((iDepth-1)/MaxDepth_vox)*MaxDepth_m;
%                         CortexReached = true;
%                         if( exist( 'bDisplay' ) && bDisplay )
%                             plot3( [HCoordsStepsX(iDepth),HCoordsStepsX(iDepth)-1], ...
%                                    [HCoordsStepsY(iDepth),HCoordsStepsY(iDepth)-1], ...
%                                    [HCoordsStepsZ(iDepth),HCoordsStepsZ(iDepth)-1], ...
%                                    'color', [1,1,1], 'linewidth', 3 );
%                         end
%                         break;
%                     end
%                 end
%             end
%             vHoles(pElem).CortexDepth = CortexDepth;
%         end
%         
%     end
%     1
    oHelmet = set_vHoles( oHelmet, vHoles );
end

function v = kspace3d(v,M)
% 3D rigid body transformation performed as shears in 1D Fourier space.
% FORMAT v1 = kspace3d(v,M)
% Inputs:
% v - the image stored as a 3D array.
% M - the rigid body transformation matrix.
% Output:
% v - the transformed image.
%
% The routine is based on the excellent papers:
% R. W. Cox and A. Jesmanowicz (1999)
% Real-Time 3D Image Registration for Functional MRI
% Submitted to MRM (April 1999) and avaliable from:
% http://varda.biophysics.mcw.edu/~cox/index.html.
% and:
% W. F. Eddy, M. Fitzgerald and D. C. Noll (1996)
% Improved Image Registration by Using Fourier Interpolation
% Magnetic Resonance in Medicine 36(6):923-931
%_______________________________________________________________________

[S0,S1,S2,S3] = shear_decomp(M);

d  = [size(v) 1 1 1];
g = 2.^ceil(log2(d));
if any(g~=d),
    tmp = v;
    v   = zeros(g);
    v(1:d(1),1:d(2),1:d(3)) = tmp;
    clear tmp;
end;

% XY-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(3)-1)/2) 0 (-g(3)/2+1):-1])/g(3);
for j=1:g(2),
    t        = reshape( exp((j*S3(3,2) + S3(3,1)*(1:g(1)) + S3(3,4)).'*tmp1) ,[g(1) 1 g(3)]);
    v(:,j,:) = real(ifft(fft(v(:,j,:),[],3).*t,[],3));
end;

% XZ-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(2)-1)/2) 0 (-g(2)/2+1):-1])/g(2);
for k=1:g(3),
    t        = exp( (k*S2(2,3) + S2(2,1)*(1:g(1)) + S2(2,4)).'*tmp1);
    v(:,:,k) = real(ifft(fft(v(:,:,k),[],2).*t,[],2));
end;

% YZ-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(1)-1)/2) 0 (-g(1)/2+1):-1])/g(1);
for k=1:g(3),
    t        = exp( tmp1.'*(k*S1(1,3) + S1(1,2)*(1:g(2)) + S1(1,4)));
    v(:,:,k) = real(ifft(fft(v(:,:,k),[],1).*t,[],1));
end;

% XY-shear
tmp1 = -sqrt(-1)*2*pi*([0:((g(3)-1)/2) 0 (-g(3)/2+1):-1])/g(3);
for j=1:g(2),
    t        = reshape( exp( (j*S0(3,2) + S0(3,1)*(1:g(1)) + S0(3,4)).'*tmp1) ,[g(1) 1 g(3)]);
    v(:,j,:) = real(ifft(fft(v(:,j,:),[],3).*t,[],3));
end;

if any(g~=d), v = v(1:d(1),1:d(2),1:d(3)); end;
end


function [S0,S1,S2,S3] = shear_decomp(A)
% Decompose rotation and translation matrix A into shears S0, S1, S2 and
% S3, such that A = S0*S1*S2*S3.  The original procedure is documented
% in:
% R. W. Cox and A. Jesmanowicz (1999)
% Real-Time 3D Image Registration for Functional MRI

A0 = A(1:3,1:3);
if any(abs(svd(A0)-1)>1e-7), error('Can''t decompose matrix'); end;


t  = A0(2,3); if t==0, t=eps; end;
a0 = pinv(A0([1 2],[2 3])')*[(A0(3,2)-(A0(2,2)-1)/t) (A0(3,3)-1)]';
S0 = [1 0 0; 0 1 0; a0(1) a0(2) 1];
A1 = S0\A0;  a1 = pinv(A1([2 3],[2 3])')*A1(1,[2 3])';  S1 = [1 a1(1) a1(2); 0 1 0; 0 0 1];
A2 = S1\A1;  a2 = pinv(A2([1 3],[1 3])')*A2(2,[1 3])';  S2 = [1 0 0; a2(1) 1 a2(2); 0 0 1];
A3 = S2\A2;  a3 = pinv(A3([1 2],[1 2])')*A3(3,[1 2])';  S3 = [1 0 0; 0 1 0; a3(1) a3(2) 1];

s3 = A(3,4)-a0(1)*A(1,4)-a0(2)*A(2,4);
s1 = A(1,4)-a1(1)*A(2,4);
s2 = A(2,4);
S0 = [[S0 [0  0 s3]'];[0 0 0 1]];
S1 = [[S1 [s1 0  0]'];[0 0 0 1]];
S2 = [[S2 [0 s2  0]'];[0 0 0 1]];
S3 = [[S3 [0  0  0]'];[0 0 0 1]];
end

