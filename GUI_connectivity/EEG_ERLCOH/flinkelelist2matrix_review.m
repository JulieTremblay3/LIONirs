function [option,idhalf2,idoption2,idhalf1] = flinkelelist2matrix_review(listelectrode)
% Input : Label Electrode listelectrode
% Link of the matrix are organised as  
% 1,2;1,3;1,4....1,n; 2,3.... to save space
% without the diagonal and the symmetry sush as nlink = (nele x nele - nele)/2
% matrix 20 x 20 = (20*20-20)/2 = 190
% option identify the link option{1).ele1 = ele1name, option{1}.ele2 =
% ele2name
% to reconstruct the matrix use 
% idhalf to reconstruct link in a matrix
% matgr1 = zeros(nele,nele);
% matgr1(idhalf)=A ;
% matgr1 = matgr1 +flipud(rot90(matgr1))
% imagesc(matgr1)
% idoption identify the indice in the matrix to option{id} location 
% matgr1(3,4) = option. 
% io = idoption(isig(id));
% option{io}.ele1,  option{io}.ele2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id = 1;
nele = numel(listelectrode)
matid = zeros(nele,nele);
matid2 = zeros(nele,nele);
    for ielex=2:nele
        ielex;
        ieley = 1;
        while ieley < ielex | ieley==1
            option{id}.ele1 =  listelectrode{ielex};
            option{id}.ele2 =listelectrode{ieley};
            option{id}.matposition = [ielex,ieley];     
            matid2(ieley,ielex)=id;
            ieley = ieley + 1;
            id = id + 1;
        end
    end
%idhalf1 = find(matid); %index dans la matrice 99x99
idhalf2 = find(matid2);
idoption2 = matid(idhalf2); 
%idoption = matid(idhalf); %index dans la matrice option 

