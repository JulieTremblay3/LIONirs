% Permet d'activer ou de desactiver un element graphique
%
% Intrants: oDispOpt  : Objet IO_DisplayOptions
%           strField  : string de l'option d'affichage.
%           bIsActive : 0 = deasactive
%                       1 = active
%
% Extrants: oDispOpt : Objet modifie
%
%**************************************************************************
function oDispOpt = set_DispOptChecked( oDispOpt, strField, bIsActive )

    ItemNo = getfield( oDispOpt.Subscripts, strField );
    oDispOpt.v_bItemChk( ItemNo ) = bIsActive;
    