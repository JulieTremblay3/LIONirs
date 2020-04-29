% fonction AlignementVecteurs
%
% Intrants: vStrt   Vecteur de depart (orientation de depart)
%           vFnsh   Vecteur d'arrivee (orientation souhaitee)
%
% Extrants: matrix4_4   Matrice de rotation 4x4
% 
% Commentaires: Fonction permettant de construire la matrice de rotation 
%               servant a aligner l'orientation d'un vecteur 3D sur un
%               deuxieme vecteur 3D.
%
% Exemple: v2 = [ x, y, z, 1] * matRot;
%
% Le vecteur utilisant la matrice de transformation doit etre un
% vecteur-colonne (4x1). La matrice doit etre a droite dans la
% multiplication.
function matrix4_4 = AlignementVecteurs( vStrt, vFnsh )


    %Ne peut fonctionner avec des vecteurs nuls. Retourner la matrice
    %identitee
    if( isempty(find(vStrt)) || isempty(find(vFnsh)) )
        matrix4_4 = [ 1 0 0 0; ...
                      0 1 0 0; ...
                      0 0 1 0; ...
                      0 0 0 1 ];
        return;
    end
    
    % Normalisation
    vFnsh = vFnsh./(vFnsh'*vFnsh)^0.5;
    vStrt = vStrt./(vStrt'*vStrt)^0.5;

    % Produit vectoriel (axe de rotation)
    vABC = cross( vFnsh,vStrt );

    % Normalisation de l'axe de rotation
    vABC = vABC./(vABC'*vABC)^0.5;

    % Determination de l'angle de rotation
    tetha = 2*atan2( ((vFnsh-vStrt)'*(vFnsh-vStrt))^0.5, ((vFnsh+vStrt)'*(vFnsh+vStrt))^0.5 );

    % Creation du Quaternion de rotation: Q=(X,Y,Z,W);
    sin_a = sin( tetha / 2 );
    cos_a = cos( tetha / 2 );
    quat_rotation = [ vABC(1) * sin_a, vABC(2) * sin_a, ...
                      vABC(3) * sin_a, cos_a ];

    %Normalisation du quaternion
    quat_rotation = quat_rotation ./ (quat_rotation*quat_rotation')^0.5;

    XX = quat_rotation(1)*quat_rotation(1);
    XY = quat_rotation(1)*quat_rotation(2);
    XZ = quat_rotation(1)*quat_rotation(3);
    XW = quat_rotation(1)*quat_rotation(4);
    YY = quat_rotation(2)*quat_rotation(2);
    YZ = quat_rotation(2)*quat_rotation(3);
    YW = quat_rotation(2)*quat_rotation(4);
    ZZ = quat_rotation(3)*quat_rotation(3);
    ZW = quat_rotation(3)*quat_rotation(4);

    %Matrice de rotation autours de l'axe de rotation
    matrix4_4 = [ 1 - 2 * ( YY + ZZ ), 2 * ( XY - ZW ),     2 * ( XZ + YW ),      0; ...
                  2 * ( XY + ZW ),     1 - 2 * ( XX + ZZ ), 2 * ( YZ - XW ),      0; ...
                  2 * ( XZ - YW ),     2 * ( YZ + XW ),     1 - 2 * ( XX + YY ),  0; ...
                  0,            0,                   0,     1  ];