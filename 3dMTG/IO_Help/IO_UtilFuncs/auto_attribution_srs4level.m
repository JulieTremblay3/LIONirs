function [not_all,msg,oHelmet] = auto_attribution_srs(oHelmet)

% Avoir une distribution automatique des sources et détecteurs
% Previous label gauche on essai de mettre de A
% droite on met C

sMtg = get_Mtg(oHelmet);
vHoles = get_vHoles(oHelmet);
nbsrc = numel(sMtg.v_pSrc);
nbdet = numel(sMtg.v_pDet);

[path, name]= uigetfile('.com','open possibility');
if name==0
    return
end
load([name,path],'-mat');
answer = (inputdlg('Number of trial to optimise solution','User define',1,{'100'}));
max_essai =  str2num(answer{1});
not_all = 0;

%Priorisation des trous  a1 à droite à l'avant, e1 droite arrière, b1
%droite à l'avant, g1 droite arrière
% maxnb_by_letter_D = zeros(26,1);
% maxnb_by_letter_G = zeros(26,1);
% for p = 1:numel(vHoles)
%     label = vHoles(p).Label;
%     if label(1) == 'D'
%         indletter = double(label(2))-64;
%         nb = str2num(label(3:end));
%         if indletter >= 1 & indletter <= 26
%             if maxnb_by_letter_D(indletter,1) < nb
%                 maxnb_by_letter_D(indletter,1) = nb;
%             end
%         end
%     else
%         indletter = double(label(2))-64;
%         nb = str2num(label(3:end));
%         if indletter >= 1 & indletter <= 26
%             if maxnb_by_letter_G(indletter,1) < nb
%                 maxnb_by_letter_G(indletter,1) = nb;
%             end
%         end
%     end
% end
typemuxsrs = zeros(numel(vHoles),1);
%4 cadran gauche avant gauche arrière droit avant droit arrière
% for p = 1:numel(vHoles)
%     label = vHoles(p).Label;
%     if label(1) == 'D'
%         indletter = double(label(2))-64;
%         nb = str2num(label(3:end));
%         if indletter >= 1 & indletter <= 26
%             if nb <  maxnb_by_letter_D(indletter)/2
%                 typemuxsrs(p) = 1 ;
%             else
%                 typemuxsrs(p) = 3;
%             end
%         end
%     else
%         indletter = double(label(2))-64;
%         nb = str2num(label(3:end));
%         if indletter >= 1 & indletter <= 26
%             if nb <  maxnb_by_letter_G(indletter)/2
%                 typemuxsrs(p) = 2 ;
%             else
%                 typemuxsrs(p) = 4;
%             end
%         end
%     end
% end



% cadran a droit a et gauche c 
for p = 1:numel(vHoles)
    label = vHoles(p).Label;
    if label(1) == 'D'
        indletter = double(label(2))-64;
        nb = str2num(label(3:end));
        if indletter >= 4 & indletter <= 26  % de d à ...
            typemuxsrs(p) = 1 ; %a1 %Right
        else
            typemuxsrs(p) = 1 ; %a1
        end
    else
        indletter = double(label(2))-64;
        nb = str2num(label(3:end));
        if indletter >= 4 & indletter <= 26
            typemuxsrs(p) = 2 ; %c1 %Left
        else
            typemuxsrs(p) = 2; %c1
        end
    end
end
%Fin de la priorisation des trous


%Attribute source
%les distances les plus proches
%les contaminations
listSrstxt =[{'b2a1'},{'b1a2'},{'b4a3'},{'b3a4'},...
    {'b6a5'},{'b5a6'},{'b8a7'},{'b7a8'},...
    {'b10a9'},{'b9a10'},{'b12a11'},{'b11a12'},...
    {'b14a13'},{'b13a14'},{'b16a15'},{'b15a16'};...
    {'d2c1'},{'d1c2'},{'d4c3'},{'d3c4'},...
    {'d6c5'},{'d5c6'},{'d8c7'},{'d7c8'},...
    {'d10c9'},{'d9c10'},{'d12c11'},{'d11c12'},...
    {'d14c13'},{'d13c14'},{'d16c15'},{'d15c16'}
    {'f2e1'},{'f1e2'},{'f4e3'},{'f3e4'},...
    {'f6e5'},{'f5e6'},{'f8e7'},{'f7e8'},...
    {'f10e9'},{'f9e10'},{'f12e11'},{'f11e12'},...
    {'f14e13'},{'f13e14'},{'f16e15'},{'f15e16'}
    {'h2g1'},{'h1g2'},{'h4g3'},{'h3g4'},...
    {'h6g5'},{'h5g6'},{'h8g7'},{'h7g8'},...
    {'h10g9'},{'h9g10'},{'h12g11'},{'h11g12'},...
    {'h14g13'},{'h13g14'},{'h16g15'},{'h15g16'}];
listSrstxt = listSrstxt';
listSrsnum = [];
for ind = 1:numel(listSrstxt)
    listSrsnum = [listSrsnum;src_s2src_n( listSrstxt{ind})];
end
listSrsnum = reshape(listSrsnum,16,4);

%Attribution 4
%Maximum d'essai en function de l'ordre des sources il est possible
%que toute les sources ne peuvent pas être utiliser, on utilise un
%random pour modifier l'ordre et faire différent test on conserve
%celui qui permet de garder le plus de source à 4 sans
%contamination
indsrs = find(sMtg.v_HolesMtg<999999);
sMtg.v_HolesMtg(indsrs) = 0;
maxintouse = 0;
maxintouse3 = 0;
hwaitbar = waitbar(0);
max_ind_source_attribue = 0;
listSrs_use = listSrsnum;

%initialisation des paramètres
ind_source_attribue_total=[];
tempind_source_attribueold = 1;
tempind_source_attribue = tempind_source_attribueold ;
tempsrsdispoold = ones(nbsrc,1);
tempsrsdispo =  tempsrsdispoold;

indsrs = find(sMtg.v_HolesMtg<999999);
sMtg.v_HolesMtg(indsrs) = 0;
tempMtgfinalold = sMtg.v_HolesMtg;
tempMtgfinal =  tempMtgfinalold;

level = 'L4'; %questdlg('How many mux level you want to place','title','L4','L3','L3');

if strcmp(level,'L4')
    indtousefinal2 = [];
    for ind_essai = 1:max_essai
        waitbar(ind_essai/max_essai*0.5,hwaitbar)
        nb_listsrs = size(listsrslevel4,1);
        randind = randperm(nb_listsrs);
        listsrslevel4rnd = listsrslevel4( randind,:);
        ind_source_attribue = tempind_source_attribueold;
        srsdispo = tempsrsdispoold;
        muxdispo = ones(16,1);
        sMtg.v_HolesMtg = tempMtgfinalold;
        for ind_possible = 1:size(listsrslevel4rnd,1)
            srs1 = listsrslevel4rnd(ind_possible,1);
            srs2 = listsrslevel4rnd(ind_possible,2);
            srs3 = listsrslevel4rnd(ind_possible,3);
            srs4 = listsrslevel4rnd(ind_possible,4);
            if srsdispo(srs1)==1 & srsdispo(srs2)==1  & srsdispo(srs3)==1 & srsdispo(srs4)==1
                srsdispo(srs1) = 0;
                srsdispo(srs2) = 0;
                srsdispo(srs3) = 0;
                srsdispo(srs4) = 0;
                p1 = sMtg.v_pSrc(srs1);
                p2 = sMtg.v_pSrc(srs2);
                p3 = sMtg.v_pSrc(srs3);
                p4 = sMtg.v_pSrc(srs4);
                [priorite muxchosen muxdispo] = auto_attribution_srs_label_4([p1,p2,p3,p4], typemuxsrs,muxdispo,ind_source_attribue);
%                 sMtg.v_HolesMtg(p1) = listSrs_use(ind_source_attribue,priorite(1));
%                 sMtg.v_HolesMtg(p2) = listSrs_use(ind_source_attribue,priorite(2));
%                 sMtg.v_HolesMtg(p3) = listSrs_use(ind_source_attribue,priorite(3));
%                 sMtg.v_HolesMtg(p4) = listSrs_use(ind_source_attribue,priorite(4));
                sMtg.v_HolesMtg(p1) = listSrs_use(muxchosen,priorite(1));
                sMtg.v_HolesMtg(p2) = listSrs_use(muxchosen,priorite(2));
                sMtg.v_HolesMtg(p3) = listSrs_use(muxchosen,priorite(3));
                sMtg.v_HolesMtg(p4) = listSrs_use(muxchosen,priorite(4));
                ind_source_attribue = ind_source_attribue +1;
                if  ind_source_attribue >= size(listSrs_use,1);
                    break
                end
            end
        end
        ind_source_attribue_total =[ind_source_attribue_total,ind_source_attribue];
        %Combinaison de 2 possible après
        ind =  find(srsdispo);
        %Peut-on retrouver une combinaison des sources dans la liste actuel
        indok2 =[];
        for i = 1:numel(ind)
            indok = find(listsrslevel2(:,1)==ind(i));
            for j = 1:numel(ind)
                indok2srs = find(listsrslevel2(indok,2)==ind(j));
                if ~isempty(indok2srs)
                    indok2 = [indok2, indok2srs+indok(1)-1];
                end
            end
        end
        tot_ind_sourceattribue = ind_source_attribue * 4  + numel(indok2)*0.01;
        if max_ind_source_attribue < tot_ind_sourceattribue
            max_ind_source_attribue= tot_ind_sourceattribue;
            tempMtgfinal = sMtg.v_HolesMtg;
            tempsrsdispo = srsdispo;
            tempmuxdispo = muxdispo; 
            tempind_source_attribue = ind_source_attribue;
            indtousefinal2 = indok2;
            essaisfinal = ind_essai; 
        end
    end
    listsrslevel2temp =  listsrslevel2;
    if ~isempty(indtousefinal2)
        listsrslevel2 =  listsrslevel2temp(indtousefinal2,:);
    else
        listsrslevel2 = [];
    end
    figure
    plot(ind_source_attribue_total-1,'x');
    title(['Nb mux placed by combinaison of 4 sources, final trial used ', num2str(essaisfinal)]);
    xlabel(['Trial'])
    ylabel(['Nb mux without contamination with a,c,e,g'])
    
end
%***********Fin possibilité de 4 contaminations
tempMtgfinalold = tempMtgfinal;
tempsrsdispoold = tempsrsdispo;
tempind_source_attribueold = tempind_source_attribue;
ind_source_attribue_total=[];
%***********Possibilité de 3 contaminations
%On fait de même avec les combinaisons à 3
if strcmp(level,'L3')
    indtousefinal2 = [];
    for  ind_essai = 1:max_essai
        waitbar(ind_essai/max_essai*0.25+0.5,hwaitbar)
        nb_listsrs = size(listsrslevel3,1);
        randind = randperm(nb_listsrs);
        sMtg.v_HolesMtg = tempMtgfinalold;
        srsdispo =  tempsrsdispoold;
        ind_source_attribue = tempind_source_attribueold;
        listsrslevel3rnd = listsrslevel3( randind,:);
        for ind_possible = 1:size(listsrslevel3rnd,1)
            srs1 = listsrslevel3rnd(ind_possible,1);
            srs2 = listsrslevel3rnd(ind_possible,2);
            srs3 = listsrslevel3rnd(ind_possible,3);
            if srsdispo(srs1)==1 & srsdispo(srs2)==1  & srsdispo(srs3)==1
                srsdispo(srs1) = 0;
                srsdispo(srs2) = 0;
                srsdispo(srs3) = 0;
                p1 = sMtg.v_pSrc(srs1);
                p2 = sMtg.v_pSrc(srs2);
                p3 = sMtg.v_pSrc(srs3);
                priorite = auto_attribution_srs_label([p1,p2,p3], typemuxsrs );
                sMtg.v_HolesMtg(p1) = listSrs_use(ind_source_attribue,priorite(1));
                sMtg.v_HolesMtg(p2) = listSrs_use(ind_source_attribue,priorite(2));
                sMtg.v_HolesMtg(p3) = listSrs_use(ind_source_attribue,priorite(3));
                ind_source_attribue = ind_source_attribue +1;
                flag_plus_de_srs = ind_source_attribue >= size(listSrs_use,1);
                if  flag_plus_de_srs
                    break
                end
            end
        end
        ind_source_attribue_total =[ind_source_attribue_total,ind_source_attribue];
        %Combinaison de 2 possible après
        ind =  find(srsdispo);
        %Peut-on retrouver une combinaison des sources dans la liste actuel
        indok2 =[];
        for i = 1:numel(ind)
            indok = find(listsrslevel2(:,1)==ind(i));
            for j = 1:numel(ind)
                indok2srs = find(listsrslevel2(indok,2)==ind(j));
                if ~isempty(indok2srs)
                    indok2 = [indok2, indok2srs+indok(1)-1];
                end
            end
        end
         tot_ind_sourceattribue = ind_source_attribue * 3  %+ numel(indok2)*0.25;  
        if max_ind_source_attribue < tot_ind_sourceattribue 
            max_ind_source_attribue=ind_source_attribue;
            tempMtgfinal = sMtg.v_HolesMtg;
            tempsrsdispo = srsdispo;
            tempind_source_attribue = ind_source_attribue;
            indtousefinal2 = indok2;
        end
    end

  listsrslevel2temp =  listsrslevel2;
    if ~isempty(indtousefinal2)
        listsrslevel2 =  listsrslevel2temp(indtousefinal2,:);
    else
        listsrslevel2 = [];
    end
    figure
    plot(ind_source_attribue_total-1,'x');       
    title(['Nb mux placed by combinaison of 4 sources']);
    xlabel(['Trial'])
    ylabel(['Nb mux without contamination with a,c,e'])
    
end

%***********Fin Possibilité de 3 contaminations
tempMtgfinalold = tempMtgfinal;
tempsrsdispoold = tempsrsdispo;
tempind_source_attribueold = tempind_source_attribue;

place2 = 0
if ~isempty(listsrslevel2)
    %***********Possibilité de 2
    for  ind_essai = 1:max_essai
        waitbar(ind_essai/max_essai*0.25+0.75,hwaitbar)
        nb_listsrs = size(listsrslevel2,1);
        randind = randperm(nb_listsrs);
        listsrslevel2rnd = listsrslevel2( randind,:);
        sMtg.v_HolesMtg = tempMtgfinalold;
        srsdispo =  tempsrsdispoold;
        muxdispo = tempmuxdispo;
        ind_source_attribue =  tempind_source_attribueold;
        for ind_possible = 1:size(listsrslevel2,1)
            srs1 = listsrslevel2rnd(ind_possible,1);
            srs2 = listsrslevel2rnd(ind_possible,2);
            if srsdispo(srs1)==1 & srsdispo(srs2)==1
                srsdispo(srs1) = 0;
                srsdispo(srs2) = 0;
                p1 = sMtg.v_pSrc(srs1);
                p2 = sMtg.v_pSrc(srs2);
                [priorite muxchosen muxdispo] = auto_attribution_srs_label_2([p1,p2], typemuxsrs ,muxdispo);
                sMtg.v_HolesMtg(p1) = listSrs_use(muxchosen,priorite(1));
                sMtg.v_HolesMtg(p2) = listSrs_use(muxchosen,priorite(2));
                ind_source_attribue = ind_source_attribue +1;
                flag_plus_de_srs = ind_source_attribue >= size(listSrs_use,1);
                if  flag_plus_de_srs
                    break
                end
            end
        end
        if max_ind_source_attribue < ind_source_attribue
             max_ind_source_attribue=ind_source_attribue;
             place2.tempMtgfinal_2 = sMtg.v_HolesMtg;
             place2.tempsrsdispo_2 = srsdispo;
             place2.tempmuxdispo_2  = muxdispo; 
             place2.tempind_source_attribue_2 = ind_source_attribue;
            break
        end
    end
end
%***********Fin possibilité de 2 contaminations

%Source seule
if isfield(place2,'tempind_source_attribue_2')
    ind_source_attribue =  place2.tempind_source_attribue_2;
    sMtg.v_HolesMtg =  place2.tempMtgfinal_2;
    srsdispo =  place2.tempsrsdispo_2;
    muxdispo =  place2.tempmuxdispo_2;
end
if ~isempty(find(srsdispo)) %verifier si toute les sources sont ok s'il reste des sources physiques indépendante disponible les utiliser
    srs1 = find(srsdispo);
    for i=1:numel(srs1)
        p1 = sMtg.v_pSrc(srs1(i));
        if ind_source_attribue <= size(listSrs_use,1)
            [priorite muxchosen, muxdispo]= auto_attribution_srs_label_1(p1,typemuxsrs, muxdispo);
            sMtg.v_HolesMtg(p1) = listSrs_use(ind_source_attribue,priorite(1));
            ind_source_attribue = ind_source_attribue +1;
            srsdispo(srs1(i))=0;
        else
            disp(['We are not able to place all srs'])%it miss srs ',num2str(find(srsdispo))])
            break
        end
    end
end

if isempty(find(srsdispo))
    not_all = 0;
else
    not_all = 1;
end
 close(hwaitbar);



%Verifier les croisement gauche et droite
% List srs gauche 
ListLeft = [listSrs_use(:,2); listSrs_use(9:16,3); listSrs_use(9:16,4)]; %c1àc16, e9àe16, g9à16
ListRight = [listSrs_use(:,1); listSrs_use(1:8,3); listSrs_use(1:8,4)]; %a1àa16, e1àe8, g1à8
reportwrong = [' Cross between fiber On right : '];
nbwrong = 0
for i=1:numel(sMtg.v_HolesMtg)
    if typemuxsrs(i) == 1; % RIGHT a
        idwrongside = find(sMtg.v_HolesMtg(i)==ListLeft); 
        if ~isempty(idwrongside)
        reportwrong = [reportwrong,' ',src_n2src_s(sMtg.v_HolesMtg(i))'];
        nbwrong = nbwrong+1;
        end
    end
end 
reportwrong = [reportwrong,'  On left : '];
for i=1:numel(sMtg.v_HolesMtg)  
    if typemuxsrs(i) == 2; % LEFT c
        idwrongside = find(sMtg.v_HolesMtg(i)==ListRight); 
        if ~isempty(idwrongside)
        reportwrong = [reportwrong,' ',src_n2src_s(sMtg.v_HolesMtg(i))'];
        nbwrong = nbwrong+1;
        end
    end
end
reportwrong = [' wrong side = ', num2str(nbwrong), reportwrong];
%
msg =  [ ' unplaced = ' num2str(numel(find(srsdispo))), reportwrong];

%        for ind = 1:numel(rightindsrs)
%            sMtg.v_HolesMtg(rightindsrs(ind)) =listSrsnum (indlistdet);
%            indlistsrs =indlistsrs + 1;
%        end
%
%        for ind = 1:numel(leftindsrs)
%             sMtg.v_HolesMtg(leftindsrs(ind)) = listSrsnum (indlistdet);
%            indlistsrs=indlistsrs + 1;
%        end
%
%
%   Non utiliser pour l'instant
%        rightindsrs = [];
%        leftindsrs =[];
%        for ind_srs = 1:nbsrc
%             labelsrs = vHoles( sMtg.v_pSrc(ind_srs)).Label;
%             if labelsrs(1)=='D'
%                 rightindsrs = [rightindsrs,sMtg.v_pSrc(ind_srs) ];
%             elseif labelsrs(2) == 'G' | labelsrs(2) == 'Z'
%                 leftindsrs=  [leftindsrs, sMtg.v_pSrc(ind_srs)];
%             end
%        end

oHelmet= set_Mtg(oHelmet,sMtg);