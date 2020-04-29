%**************************************************************************
% FONCTION : HELMET/set_Mtg 
%
% INTRANTS : obj          -> objet Helmet
%            new_Mtg      -> nouvelle structure d'information de montage

%
% EXTRANTS : Objet Helmet modifié.
%
% UTILITE : Fonction permettant de modifier les paramètres d'un montage.
%**************************************************************************
function obj = set_Mtg(obj, new_Mtg)

    obj.Mtg_Data = new_Mtg;
    
    %Classement des Sources, par leur indice lineaire. (0 en dernier)
    matTemp = sort( obj.Mtg_Data.v_pSrc );
    obj.Mtg_Data.v_pSrc = matTemp(find(matTemp));
    
    %Classement des Detecteurs, par leur indice lineaire. (0 en dernier)
    matTemp = sort( obj.Mtg_Data.v_pDet );
    obj.Mtg_Data.v_pDet = matTemp(find(matTemp));
    
    %Classement des Electrodes, par leur indice lineaire. (0 en dernier)
    matTemp = sort( obj.Mtg_Data.v_pEle );
    obj.Mtg_Data.v_pEle = matTemp(find(matTemp));
    
    %Regrouppement des elements (Sources, Detecteurs, Electrodes)
    v_pValid_Elems = [obj.Mtg_Data.v_pSrc,obj.Mtg_Data.v_pDet,obj.Mtg_Data.v_pEle];
    
    %Ajustement de la taille de la matrice de trous
    if( numel(obj.Mtg_Data.v_HolesMtg) < max(v_pValid_Elems) )
        obj.Mtg_Data.v_HolesMtg(max(v_pValid_Elems)) = 0;
    end
    
    %Retrait des elements du montage qui ne sont pas dans la liste de
    %src/det/ele
    tmp_v_HolesMtg = obj.Mtg_Data.v_HolesMtg;
    obj.Mtg_Data.v_HolesMtg = zeros(size(tmp_v_HolesMtg));
    %disp( sprintf( 'Numel(v_HolesMtg):%d', numel(obj.Mtg_Data.v_HolesMtg) ) );
    %v_pValid_Elems
    obj.Mtg_Data.v_HolesMtg(v_pValid_Elems) = tmp_v_HolesMtg(v_pValid_Elems);
  