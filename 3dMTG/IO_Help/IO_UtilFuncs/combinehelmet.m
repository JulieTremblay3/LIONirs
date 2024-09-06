function helmet1 = combinehelmet(helmet1,helmet2)
%helmet1 : contient un casque sans montage avec les vHoles du fichier digitalises
%helmet2 : contient un casque avec montage avec les vHoles du fichier d'un
%autre projet
    %Le montage !!! 
    sMtg2 = get_Mtg( helmet2);
    vHoles2 = get_vHoles( helmet2);
    
    %Le casque
    sMtg1 = get_Mtg( helmet1 );
    vHoles1 = get_vHoles( helmet1);
    bad = [];
    ifbad = 0;
    %LABEL que l'on recherche 
    %Les sources
    for i = 1:numel(sMtg2.v_pSrc)
        p_2 = sMtg2.v_pSrc(i);
        ifbad = 0;
        for p_1 = 1:numel(vHoles1)
            if strcmp(vHoles1(p_1).Label,vHoles2(p_2).Label)
                 sMtg1.v_HolesMtg(p_1)= sMtg2.v_HolesMtg(p_2);
                 sMtg1.v_pSrc = [sMtg1.v_pSrc, p_1];
                 ifbad = 1; %la source est placé
                 break
            end
        end
        if ifbad == 0 %La source n'est pas placé
             bad = [bad,' ', vHoles2(p_2).Label];
        end
    end
    %Les détecteurs
    for i = 1:numel(sMtg2.v_pDet)
        p_2 =  sMtg2.v_pDet(i);
        ifbad = 0;
        for p_1 = 1:numel(vHoles1)
            if strcmp(vHoles1(p_1).Label,vHoles2(p_2).Label)
                 sMtg1.v_HolesMtg(p_1)= sMtg2.v_HolesMtg(p_2);
                 sMtg1.v_pDet = [sMtg1.v_pDet, p_1];
                 ifbad = 1; %la source est placé
                 break               
            end              
        end
        if ifbad == 0 %La source n'est pas placé
             bad = [bad,' ', vHoles2(p_2).Label];
        end
    end
    %Les électrodes
    sMtg1.v_HolesEle = cell(1,numel(sMtg1.v_HolesMtg));
    for i = 1:numel(sMtg2.v_pEle)
        p_2 = sMtg2.v_pEle(i);
        ifbad = 0;
        for p_1 = 1:numel(vHoles1)
            if strcmp(vHoles1(p_1).Label,vHoles2(p_2).Label)
                 sMtg1.v_HolesEle{p_1}= sMtg2.v_HolesEle{p_2};
                 sMtg1.v_pEle = [sMtg1.v_pEle, p_1];
                 ifbad = 1; %la source est placé
                 break
            end              
        end
        if ifbad == 0 %La source n'est pas placé
             bad = [bad,' ', vHoles2(p_2).Label];
        end
    end
    if isempty(bad)   
        
    else
        h= msgbox(['Be carefull the hole '  bad ' are not found'])
        uiwait(h)
    end
    helmet1 = set_Mtg(helmet1, sMtg1);
    helmet1 = set_vHoles(helmet1,vHoles1);
    
%     %on regarde pour chaque source du montage helmet2
%     for p_1 = 1:numel(vHoles1)
%         for p_2 = 1:numel(sMtg2.v_HolesMtg)        
%             if strcmp(vHoles2(p_2).Label,vHoles1(p_1).Label)
%                 sMtg1.v_HolesMtg(p_1)= sMtg2.v_HolesMtg(p_2);
%                 if sMtg2.v_HolesMtg(p_2)> 999999
%                     sMtg1.v_pDet = [sMtg1.v_pDet, p_1];
%                 elseif sMtg2.v_HolesMtg(p_2)> 0 & sMtg2.v_HolesMtg(p_2) < 999999
%                     sMtg1.v_pSrc = [sMtg1.v_pSrc, p_1];
%                 end    
%                 break
%             end
%         end
%     end
%      sMtg1.v_HolesEle = cell(1,numel(vHoles2));
%     for p_ind_sMtg2 = 1:numel(sMtg2.v_pEle) 
%         p_2 = sMtg2.v_pEle(p_ind_sMtg2);       
%         for p_1 = 1:numel(vHoles1)
%             if strcmp(vHoles2(p_2).Label,vHoles1(p_1).Label)
%                 sMtg1.v_HolesEle{p_1}= sMtg2.v_HolesEle{p_2};                
%                 sMtg1.v_pEle = [sMtg1.v_pEle, p_1];                
%                 break
%             end
%         end
%     end         
