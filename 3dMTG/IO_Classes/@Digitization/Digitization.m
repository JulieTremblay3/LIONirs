% Classe de gestion des options d'affichage
function oDig = Digitization(args)

    %Constructeur par defaut
    if( nargin == 0 )
        %                       [x,y,z,1] (last column for 4x4 matrices)
        oDig.vProbes(1).Coord = [0,0,0,0]; 
        oDig.vProbes(1).Type = 0;
        oDig.vProbes(1).Label = '';
        oDig.matFiducials = zeros(3,3); %[ - NA  - ; - LPA - ; - RPA - ]
        
        oDig = class( oDig, 'Digitization' );
        
	%Constructeur copieur
    elseif( isa(args, 'IO_DisplayOptions') )
        oDig.vProbes = args.vProbes;
        oDig.matFiducials = args.matFiducials;
        
        oDig = class( oDig, 'Digitization' );

    %Constructeur parametrique
    else
        disp( 'Digitization class: Parametric constructor does not exists' );
    end
    
    