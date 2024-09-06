%Permet de changer un mux a1 par un autre a2, tout les autres sources
%du même mux x(sont automatiquement changer
function [oHelmet] = inverse_srs(oHelmet)

sMtg = get_Mtg(oHelmet);
mux = inputdlg('Enter mux to change on each line example a1 and c1 write as a1:c1','Personnal mux',1,{'a1:c1'})
srs = mux{1};
[muxa,muxb]= strtok(srs,':');
muxa = muxa;
muxb = muxb(2:end);
mtg = sMtg.v_HolesMtg(:);
if numel(muxa)==1
    nbsrsa = (double(muxa)-64) * 1000000;
    pa = find(nbsrsa==mtg); 
    aisdet = 1;
else
    banka = muxa(1:1);
    nba= str2num(muxa(2:end));
    if banka=='a'|banka=='b'
        nbsrsa = nba;
    elseif banka=='c'|banka=='d'
        nbsrsa = nba+32;
    elseif banka=='e'|banka=='f'
        nbsrsa = nba+64;
    elseif banka=='g'|banka=='h'
        nbsrsa = nba+96;
    end
    pa = find(nbsrsa==mod(mtg,1000)|nbsrsa==floor(mtg/1000));
    aisdet = 0;  
end

if numel(muxb)==1
    nbsrsb = (double(muxb)-64) * 1000000;
    pb = find(nbsrsb==mtg);
    bisdet = 1;
else
    bankb = muxb(1:1);
    nbb= str2num(muxb(2:end));
    if bankb=='a'|bankb=='b'
        nbsrsb = nbb;
    elseif bankb=='c'|bankb=='d'
        nbsrsb = nba+32;
    elseif bankb=='e'|bankb=='f'
        nbsrsb = nbb+64;
    elseif bankb=='g'|bankb=='h'
        nbsrsb = nbb+96;
    end
    pb = find(nbsrsb==mod(mtg,1000)|nbsrsb==floor(mtg/1000));
    bisdet = 0;
end

if bisdet == aisdet
    %pas de changement de place pour une source et un détecteur
    
else % on inverse
    if aisdet
        i = find(pa == sMtg.v_pDet);
        sMtg.v_pDet(i)=[]; 
        sMtg.v_pSrc = [sMtg.v_pSrc,pa];
        i = find(pb == sMtg.v_pSrc);
        sMtg.v_pSrc(i)=[];
        sMtg.v_pDet = [sMtg.v_pDet,pb];
    else
        i = find(pb == sMtg.v_pDet);
        sMtg.v_pDet(i)=[]; 
        sMtg.v_pSrc = [sMtg.v_pSrc,pb];
        i = find(pa == sMtg.v_pSrc);
        sMtg.v_pSrc(i)=[];
        sMtg.v_pDet = [sMtg.v_pDet,pa];
    end
end

if isempty(pa)| isempty(pb)
    msgbox('Sorry srs not found, no change')
    return
end

   
sMtg.v_HolesMtg(pb)=mtg(pa);
sMtg.v_HolesMtg(pa)=mtg(pb);
oHelmet= set_Mtg(oHelmet,sMtg);

% %listsrsname 
% listsrs = mod(listsrsnum,1000);
% mtg1 = mod(mtg,1000);
% p1 = zeros(1,4);
% p2 = zeros(1,4);
% listsrsA = zeros(1,4);
% listsrsB = zeros(1,4);
% for i = 1:numel(lista)
%     p = find(mtg1==lista(i));  
%     if ~isempty(p)
%     p1(i) = p;
%     end
%     p = find(listsrs ==lista(i))
%     listsrsA(i) = listsrsnum(p)
%     p = find(mtg1==listb(i));        
%     if ~isempty(p)
%     p2(i) = p;
%     end 
%     p = find(listsrs ==listb(i))
%     listsrsB(i) = listsrsnum(p)
% end
% mtgtemp = mtg;
% 
% for i = 1:4
%     if p2(i)~= 0
%         sMtg.v_HolesMtg(p2(i))= listsrsA(i); 
%     end
%     if p1(i)~= 0
%         sMtg.v_HolesMtg(p1(i))= listsrsB(i);       
%     end              
% end
