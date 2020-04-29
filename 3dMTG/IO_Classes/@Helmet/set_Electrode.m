function [Helm, bOk] = set_Electrode(Helm,pEle)


sMtg = get_Mtg( Helm );  
if strcmp(sMtg.Gen_Params.ElectrodeType,'10-10') 
    list = liste_electrode(1);
elseif strcmp(sMtg.Gen_Params.ElectrodeType,'E128') 
    list = liste_electrode(2);
elseif strcmp(sMtg.Gen_Params.ElectrodeType,'Custom') 
    list =liste_electrode(3);
    if ~iscell(list)
     list ={list};
    end
end
%Enlever une electrode
if  ~isempty(sMtg.v_HolesEle{pEle})
    sMtg.v_HolesEle{pEle}=[];
    bOk = true;
    Helm = set_Mtg(Helm, sMtg);
    return
end   
if isempty(find(sMtg.v_pEle==pEle))
    bOk = true;
    Helm = set_Mtg(Helm, sMtg);
    return
end
%Ajout d'une électrode
% On utilise le premier element de la liste qui n'a pas été utilisé
isuse = 0;
for i = 1:numel(list) 
   for ind = 1:numel(sMtg.v_HolesEle)
        if iscellstr(sMtg.v_HolesEle(ind))
            if strmatch(list{i}, sMtg.v_HolesEle(ind)) 
               isuse = 1 ;          
            end
        end
    end
    if ~isuse
        label = list{i};
        break
    end
    isuse = 0;
end  
sMtg.v_HolesEle{pEle} = label;
bOk = true;
Helm = set_Mtg(Helm, sMtg);


%  %Verifier la possibilite du montage
%     EleFiberSelected = get( handles.popup_Electrode, 'Value' );
%     liste_ele =liste_electrode
%     label_ele = liste_ele(EleFiberSelected )
