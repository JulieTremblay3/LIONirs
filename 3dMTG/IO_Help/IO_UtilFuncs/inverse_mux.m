%Permet de changer un mux a1 par un autre a2, tout les autres sources
%du même mux x(sont automatiquement changer
function [oHelmet] = inverse_mux(oHelmet)

sMtg = get_Mtg(oHelmet);
mux = inputdlg('Enter mux to change on each line example a1 and a2 write as 1:2','Personnal mux',1,{'1:2'})
srs = mux{1};
[muxa,muxb]= strtok(srs,':');
muxa = str2num(muxa);
muxb = str2num(muxb(2:end));


mtg = sMtg.v_HolesMtg(:);

%a1, c1, e1, g1
lista = [muxa, muxa+32, muxa+64, muxa+96];
listb = [muxb, muxb+32, muxb+64, muxb+96];

listsrsnum = [18001,
                20003
                22005
                24007
                26009
                28011
                30013
                32015
                17002
                19004
                21006
                23008
                25010
                27012
                29014
                31016
                50033
                52035
                54037
                56039
                58041
                60043
                62045
                64047
                49034
                51036
                53038
                55040
                57042
                59044
                61046
                63048
                82065
                84067
                86069
                88071
                90073
                92075
                94077
                96079
                81066
                83068
                85070
                87072
                89074
                91076
                93078
                95080
                114097
                116099
                118101
                120103
                122105
                124107
                126109
                128111
                113098
                115100
                117102
                119104
                121106
                123108
                125110
                127112];
          

%listsrsname 
listsrs = mod(listsrsnum,1000);
mtg1 = mod(mtg,1000);
p1 = zeros(1,4);
p2 = zeros(1,4);
listsrsA = zeros(1,4);
listsrsB = zeros(1,4);
for i = 1:numel(lista)
    p = find(mtg1==lista(i));  
    if ~isempty(p)
    p1(i) = p;
    end
    p = find(listsrs ==lista(i))
    listsrsA(i) = listsrsnum(p)
    p = find(mtg1==listb(i));        
    if ~isempty(p)
    p2(i) = p;
    end 
    p = find(listsrs ==listb(i))
    listsrsB(i) = listsrsnum(p)
end
mtgtemp = mtg;

for i = 1:4
    if p2(i)~= 0
        sMtg.v_HolesMtg(p2(i))= listsrsA(i); 
    end
    if p1(i)~= 0
        sMtg.v_HolesMtg(p1(i))= listsrsB(i);       
    end              
end
oHelmet= set_Mtg(oHelmet,sMtg);