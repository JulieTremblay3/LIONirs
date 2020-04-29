%**************************************************************************
% CLASSE : Inter_Dlg_Comm (Constructeur)
%
% Utilité : Objet servant a memoriser les information concernant l'etat de
%           la fenetre principale. Contient le vecteur de sélections
%           courantes, le vecteur de booléen de visibilité de positions  
%           (la visibilité des élements dépend de la position de la caméra),
%           les informations de souris, l'indice du curseur
%           d'information ainsi que les indices des éléments mesurés avec
%           la règle.
%           
%**************************************************************************
function Obj = Inter_Dlg_Comm(args)

    %Constructeur par defaut
    if( nargin == 0 )
        
        %*****************************************************************
        %Structure contenant les infos du graphique et de ses éléments
        Obj.HelmetAxeDisp.IsAxeRotating = false; %Flag de camorbit en cours
        Obj.HelmetAxeDisp.axe_handle = 0; % Handle du graphique de casque
        
        %Indices lineaires dans 'Helmet::vHoles(p)'
        Obj.HelmetAxeDisp.pIdentification = 0; %Pour l'affichage des infos
        Obj.HelmetAxeDisp.pRulerFirst = 0;  %1er pt du calcul de distance
        Obj.HelmetAxeDisp.pRulerLast = 0;   %2è  pt du calcul de distance
                
        %Vecteur contenant les indices lineaires de selections pour 'Helmet::vHoles(p)'
        Obj.HelmetAxeDisp.v_pSelection = [];
        
        %Vecteur de visibilite des trous du casque. 
        %(Indices linéaires à correspondance directe pour Helmet::vHoles)
        Obj.HelmetAxeDisp.v_bVisible = [];
        
        %Vecteur de visibilite des labels du casque. 
        %(Indices linéaires à correspondance directe pour Helmet::vHoles)
        Obj.HelmetAxeDisp.v_bDisplayLabels = [];
        
        %*****************************************************************
        %Structure contenant les infos de souris
        Obj.MouseData.IsLBtnDown = 0;   %Flag de btn gauche enfoncé
        Obj.MouseData.IsRBtnDown = 0;   %Flag de btn droit enfoncé
        Obj.MouseData.LastPos = [0,0]; %Derniere position lorsque BtnDown
        Obj.MouseData.Mode = 01; %Doit etre une valeur parmi "Obj.MouseMode.Subscripts"
        Obj.MouseMode.Subscripts.CamOrbit = 01;
        Obj.MouseMode.Subscripts.CamZoom  = 02;
        Obj.MouseMode.Subscripts.Ruler    = 03;
        Obj.MouseMode.Subscripts.AddSrc   = 04;
        Obj.MouseMode.Subscripts.AddDet   = 05;
        Obj.MouseMode.Subscripts.AddEle   = 06;
        Obj.MouseMode.Subscripts.GetInf   = 07;
        
        
        %*****************************************************************
        %Structure contenant les infos de clavier
        Obj.KeyboardData.LastKey = 0; %Derniere touche enfoncee.
        Obj.KeyboardData.LKClock = [0,0,0,0,0,0]; %Instant d'enfoncement de la touche.
        
        Obj = class( Obj, 'Inter_Dlg_Comm' );
        
	%Constructeur copieur
    elseif( isa(args, 'Inter_Dlg_Comm') )
        Obj.HelmetAxeDisp    = args.HelmetAxeDisp;
        Obj.MouseData        = args.MouseData;
        Obj.KeyboardData     = args.KeyboardData;
        Obj = class( Obj, 'Inter_Dlg_Comm' );
        
    %Constructeur parametrique
    else
        disp( 'Constructeur parametrique inexistant' );
    end
    
    