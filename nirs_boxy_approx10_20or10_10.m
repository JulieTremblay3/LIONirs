function SD = nirs_boxy_approx10_20or10_10(SD,use10_10system)
%Approximate channels by 10-20 electrode system

%this function can be simplified as much of the info is in SD
%test 1212
%v_pSrc = mtg.v_pSrc; %alias - positions sur montage
%v_pDet = mtg.v_pDet;
%v_pEle = mtg.v_pEle; %Electrodes
%v_HolesEle = mtg.v_HolesEle; %Noms des électrodes
%canaux
Cidx=0;
n_s=size(SD.SrcPos,1); %nombre de sources
n_d=size(SD.DetPos,1); %nombre de détecteurs

%If SD.ElePos not defined or of size less than 19, then use 
%montage from one patient (PEV) as the default.

 ElePosPEV =[0.503635980484236,-6.42109842933709,4.52769093551111;
    4.87042118935769,-5.90352152277156,2.75687653968204;
    -5.14295598304066,-4.78946764994081,4.55389809361421;
    5.20471181432292,-4.65633690653410,7.63022287598458;
    0.480332535817484,-4.68576325897243,9.50624027238281;
    -3.25228712922376,-4.15631667963710,9.72167812006713;
    10.2505497383691,-2.52434255554954,2.39444260110667;
    -7.13844916128474,-1.60179708791733,5.76800266637413;
    8.20976246987310,-0.309413602779901,9.05673688754269;
    2.36402374251312,-0.730086711610575,11.4884040728641;
    -3.85101541886160,-0.585225718767958,10.9030687956262;
    10.0523770901632,3.17933431493396,2.91271002512250;
    -6.74489475628845,2.81324613587249,6.55340691642697;
    5.16627056139820,4.32749222610975,9.09834754533281;
    0.432123235493804,3.80993884143206,10.9192205462793;
    -3.24176689763683,3.60126520281447,10.4535878083220;
    -4.44197773122473,6.07855033999104,5.99562036942143;
    5.34301849451918,6.55805824811752,4.00117453208372;
    0.793768743730182,6.81160469302323,6.12446029378749;];
  EleNamesPEV = {'T8';'F8';'P8';'F4';'C4';'P4';'Fp2';'O2';'Fz';'Cz';'Pz';'FP1';'O1';'F3';'C3';'P3';'P7';'F7';'T7';};
  
try 
    SD.ElePos; 
catch    
  SD.ElePos = ElePosPEV; 
  SD.EleNames = EleNamesPEV;
end

n_e=length(SD.ElePos); %number of electrodes  

%First replace electrode names by slightly more standard names
for Eidx=1:n_e
   if strcmp(SD.EleNames{Eidx,1}, 'T3'), SD.EleNames{Eidx,1} = 'T7'; end
   if strcmp(SD.EleNames{Eidx,1}, 'T4'), SD.EleNames{Eidx,1} = 'T8'; end
   if strcmp(SD.EleNames{Eidx,1}, 'T5'), SD.EleNames{Eidx,1} = 'P7'; end
   if strcmp(SD.EleNames{Eidx,1}, 'T6'), SD.EleNames{Eidx,1} = 'P8'; end
end

%Loop through electrodes to fill in the missing ones
if n_e < SD.MaxElectrodes %Was previously hard-coded number of electrodes
    diff_e = setdiff(EleNamesPEV,SD.EleNames);
    for Eidx=1:length(diff_e)
        for Eidx2=1:length(EleNamesPEV)
            if strcmp(diff_e{Eidx,1}, EleNamesPEV{Eidx2,1})
                %add an electrode and its name
                SD.ElePos = [SD.ElePos; ElePosPEV(Eidx2,:)];
                SD.EleNames{n_e+Eidx,1} = EleNamesPEV{Eidx2};
            end
        end
    end
end

%Sort: Change order in SD.ElePos to make it consistent with the template: PEV 
%see required ordering in interpolate10_10system
for Eidx=1:length(EleNamesPEV) %SD.EleNames and EleNamesPEV should be the same length 
    for Eidx2=1:length(SD.EleNames)
        if strcmp(SD.EleNames{Eidx2,1}, EleNamesPEV{Eidx,1})  
            tmp_pos(Eidx,:) = SD.ElePos(Eidx2,:);
            tmp_names{Eidx,1} = SD.EleNames{Eidx2,1};
        end
    end
end
%copy arrays
for Eidx=1:length(EleNamesPEV)
    SD.ElePos(Eidx,:) = tmp_pos(Eidx,:);
    SD.EleNames{Eidx,1} = tmp_names{Eidx,1};
end
        
%find number of channels
for Sidx = 1:n_s %64 est le nombre maximal de sources alloué???
   Cidx=Cidx+length(SD.Pairs{Sidx});       
end %end for
n_c=Cidx; 
%allocate space
SD.ChnPos=zeros(n_c,3);
SD.ChnDist=zeros(n_c,1);
%loop over channels, to calculate their position
Cidx=0;
for Sidx = 1:n_s 
   dets = SD.Pairs{Sidx};
   for j=1:length(dets)
        Cidx=Cidx+1;
        SD.cn{Cidx,1}=[Sidx dets(j)]; %source detector pairs
        %linear interpolation between source and detector
        SD.ChnPos(Cidx,:)=0.5*(SD.SrcPos(Sidx,:)+SD.DetPos(dets(j),:));
        SD.ChnDist(Cidx)=SD.Dist{Sidx,1}(j);
   end %end for        
end %end for
 
if use10_10system
    %define 10_10 system
    [SD.ElePosExt, SD.EleNamesExt] = interpolate10_10system(SD.ElePos);
    [SD.ChnNames, SD.ChnMinDist, SD.EleChnCount, SD.EleChnPairs]...
        = associate_names_to_channels(SD.ChnPos,SD.ElePosExt, SD.EleNamesExt);
else %use 10_20 system
    %find electrode closest to each channel
    [SD.ChnNames, SD.ChnMinDist, SD.EleChnCount, SD.EleChnPairs]...
        = associate_names_to_channels(SD.ChnPos,SD.ElePos, SD.EleNames);   
end %end if use10_10system
end

function [ChnNames, ChnMinDist, EleChnCount, EleChnPairs]= ...
    associate_names_to_channels(ChnPos,ElePos,EleNames)
%For each channel, find standard electrode position which is closest, and
%assign the electrode name; keep track of how many channels are closest to
%each electrode
n_c=size(ChnPos,1);
n_e=size(ElePos,1);

dist=zeros(n_e,1);
EleChnCount=zeros(n_e,1);
ChnNames{n_c,1}='';
ChnMinDist=zeros(n_c,1);
EleChnPairs{n_e,1}='';

for Cidx=1:n_c
    c_x=ChnPos(Cidx,1);
    c_y=ChnPos(Cidx,2);
    c_z=ChnPos(Cidx,3);
    for Eidx=1:n_e 
        e_x= ElePos(Eidx,1);
        e_y= ElePos(Eidx,2);
        e_z= ElePos(Eidx,3);
        dist(Eidx)=((c_x-e_x)^2+(c_y-e_y)^2+(c_z-e_z)^2)^(0.5);
        %faster: dist(Eidx) = abs(c_x-e_x)+abs(c_y-e_y)+abs(c_z-e_z);
    end %end for Eidx
    %find minimum distance
    [min_value,min_e]=min(dist);
    %assign label to closest channel
    ChnNames{Cidx,1}=EleNames{min_e};
    %count how many channels have been assigned to each electrode
    EleChnCount(min_e)=EleChnCount(min_e)+1; 
    ChnMinDist(Cidx)=min_value; %keep track of minimum value found
    %add to pairing of 10-10 or 10-20 positions with channels
    EleChnPairs{min_e,1}=[EleChnPairs{min_e,1}; Cidx]; 
end %end for Cidx
end



function [ElePosExt, EleNamesExt] = interpolate10_10system(ElePos)
%uses position of 10_20 electrodes to get coordinates of 10_10 electrodes
%Note: a few electrode names are unusual in the HSJ montage
%HSJ T3, T4 are more usually T7, T8 and
%HSJ T5, T6 are more usually P7, P8
%SD.EleNames 
%1: T4 -> T8
%2: F8
%3: T6 -> P8
%4: F4
%5: C4
%6: P4
%7: Fp2
%8: O2
%9: Fz
%10: Cz
%11: Pz
%12: Fp1 (inconsistently written FP1)
%13: O1
%14: F3
%15: C3
%16: P3
%17: T5 -> P7
%18: F7
%19: T3 -> T7

%Use order of channels of standard montage in EEG fMRI recordings
% Ch1=Fp1,,0.5,µV
% Ch2=Fp2,,0.5,µV
% Ch3=F3,,0.5,µV
% Ch4=F4,,0.5,µV
% Ch5=C3,,0.5,µV
% Ch6=C4,,0.5,µV
% Ch7=P3,,0.5,µV
% Ch8=P4,,0.5,µV
% Ch9=O1,,0.5,µV
% Ch10=O2,,0.5,µV
% Ch11=F7,,0.5,µV
% Ch12=F8,,0.5,µV
% Ch13=T7,,0.5,µV
% Ch14=T8,,0.5,µV
% Ch15=P7,,0.5,µV
% Ch16=P8,,0.5,µV
% Ch17=Fz,,0.5,µV
% Ch18=Cz,,0.5,µV
% Ch19=Pz,,0.5,µV
% Ch20=Oz,,0.5,µV
% Ch21=FC1,,0.5,µV
% Ch22=FC2,,0.5,µV
% Ch23=CP1,,0.5,µV
% Ch24=CP2,,0.5,µV
% Ch25=FC5,,0.5,µV
% Ch26=FC6,,0.5,µV
% Ch27=CP5,,0.5,µV
% Ch28=CP6,,0.5,µV
% Ch29=TP9,,0.5,µV
% Ch30=TP10,,0.5,µV
% Ch31=EOG,,0.5,µV
% Ch32=ECG,,0.5,µV
% Ch33=F1,,0.5,µV
% Ch34=F2,,0.5,µV
% Ch35=C1,,0.5,µV
% Ch36=C2,,0.5,µV
% Ch37=P1,,0.5,µV
% Ch38=P2,,0.5,µV
% Ch39=AF3,,0.5,µV
% Ch40=AF4,,0.5,µV
% Ch41=FC3,,0.5,µV
% Ch42=FC4,,0.5,µV
% Ch43=CP3,,0.5,µV
% Ch44=CP4,,0.5,µV
% Ch45=PO3,,0.5,µV
% Ch46=PO4,,0.5,µV
% Ch47=F5,,0.5,µV
% Ch48=F6,,0.5,µV
% Ch49=C5,,0.5,µV
% Ch50=C6,,0.5,µV
% Ch51=P5,,0.5,µV
% Ch52=P6,,0.5,µV
% Ch53=AF7,,0.5,µV
% Ch54=AF8,,0.5,µV
% Ch55=FT7,,0.5,µV
% Ch56=FT8,,0.5,µV
% Ch57=TP7,,0.5,µV
% Ch58=TP8,,0.5,µV
% Ch59=PO7,,0.5,µV
% Ch60=PO8,,0.5,µV
% Ch61=Fpz,,0.5,µV
% Ch62=AFz,,0.5,µV
% Ch63=CPz,,0.5,µV
% Ch64=POz,,0.5,µV
ElePosExt =zeros(64,3);
EleNamesExt{64,1}='';

ElePosExt(1,:) = ElePos(12,:); EleNamesExt{1,1}= 'Fp1';
ElePosExt(2,:) = ElePos(7,:); EleNamesExt{2,1}= 'Fp2';
ElePosExt(3,:) = ElePos(14,:); EleNamesExt{3,1}= 'F3';
ElePosExt(4,:) = ElePos(4,:); EleNamesExt{4,1}= 'F4';
ElePosExt(5,:) = ElePos(15,:); EleNamesExt{5,1}= 'C3';
ElePosExt(6,:) = ElePos(5,:); EleNamesExt{6,1}= 'C4';
ElePosExt(7,:) = ElePos(16,:); EleNamesExt{7,1}= 'P3';
ElePosExt(8,:) = ElePos(6,:); EleNamesExt{8,1}= 'P4';
ElePosExt(9,:) = ElePos(13,:); EleNamesExt{9,1}= 'O1';
ElePosExt(10,:) = ElePos(8,:); EleNamesExt{10,1}= 'O2';
ElePosExt(11,:) = ElePos(18,:); EleNamesExt{11,1}= 'F7';
ElePosExt(12,:) = ElePos(2,:); EleNamesExt{12,1}= 'F8';
ElePosExt(13,:) = ElePos(19,:); EleNamesExt{13,1}= 'T7';
ElePosExt(14,:) = ElePos(1,:); EleNamesExt{14,1}= 'T8';
ElePosExt(15,:) = ElePos(17,:); EleNamesExt{15,1}= 'P7';
ElePosExt(16,:) = ElePos(3,:); EleNamesExt{16,1}= 'P8';
ElePosExt(17,:) = ElePos(9,:); EleNamesExt{17,1}= 'Fz';
ElePosExt(18,:) = ElePos(10,:); EleNamesExt{18,1}= 'Cz';
ElePosExt(19,:) = ElePos(11,:); EleNamesExt{19,1}= 'Pz';
ElePosExt(20,:) = 0.5*(ElePos(13,:)+ElePos(13,:)); EleNamesExt{20,1}= 'Oz'; %between O1 & O2
ElePosExt(21,:) = 0.5*(ElePos(14,:)+ElePos(10,:)); EleNamesExt{21,1}= 'FC1'; %between F3 & Cz
ElePosExt(22,:) = 0.5*(ElePos(4,:)+ElePos(10,:)); EleNamesExt{22,1}= 'FC2'; %between F4 & Cz
ElePosExt(23,:) = 0.5*(ElePos(16,:)+ElePos(10,:)); EleNamesExt{23,1}= 'CP1'; %between P3 & Cz
ElePosExt(24,:) = 0.5*(ElePos(6,:)+ElePos(10,:)); EleNamesExt{24,1}= 'CP2'; %between P4 & Cz
ElePosExt(25,:) = 0.25*(ElePos(18,:)+ElePos(14,:)+ElePos(19,:)+ElePos(15,:)); EleNamesExt{25,1}= 'FC5'; %between F7,F3,T7,C3 
ElePosExt(26,:) = 0.25*(ElePos(2,:)+ElePos(4,:)+ElePos(5,:)+ElePos(1,:)); EleNamesExt{26,1}= 'FC6'; %between F8,F4,T8,C4
ElePosExt(27,:) = 0.25*(ElePos(19,:)+ElePos(15,:)+ElePos(17,:)+ElePos(16,:)); EleNamesExt{27,1}= 'CP5'; %between T7,C3,P7,P3
ElePosExt(28,:) = 0.25*(ElePos(1,:)+ElePos(5,:)+ElePos(3,:)+ElePos(6,:)); EleNamesExt{28,1}= 'CP6'; %between T8,C4,P8,P4
%ElePosExt(29,:) see below 'TP9'; %beyond TP7
%ElePosExt(30,:) see below 'TP10'; %beyond TP8
ElePosExt(31,:) = [0 0 0]; EleNamesExt{31,1}= 'EOG';
ElePosExt(32,:) = [0 0 0]; EleNamesExt{32,1}= 'ECG';
ElePosExt(33,:) = 0.5*(ElePos(14,:)+ElePos(9,:)); EleNamesExt{33,1}= 'F1'; %between F3 & Fz
ElePosExt(34,:) = 0.5*(ElePos(4,:)+ElePos(9,:)); EleNamesExt{34,1}= 'F2'; %between F4 & Fz
ElePosExt(35,:) = 0.5*(ElePos(15,:)+ElePos(10,:)); EleNamesExt{35,1}= 'C1'; %between C3 & Cz
ElePosExt(36,:) = 0.5*(ElePos(5,:)+ElePos(10,:)); EleNamesExt{36,1}= 'C2'; %between C4 & Cz
ElePosExt(37,:) = 0.5*(ElePos(16,:)+ElePos(11,:)); EleNamesExt{37,1}= 'P1'; %between P3 & Pz
ElePosExt(38,:) = 0.5*(ElePos(6,:)+ElePos(11,:)); EleNamesExt{38,1}= 'P2'; %between P4 & Pz
ElePosExt(39,:) = 0.5*(ElePos(12,:)+ElePos(14,:)); EleNamesExt{39,1}= 'AF3'; %between Fp1 & F3
ElePosExt(40,:) = 0.5*(ElePos(7,:)+ElePos(4,:)); EleNamesExt{40,1}= 'AF4'; %between Fp2 & F4
ElePosExt(41,:) = 0.5*(ElePos(14,:)+ElePos(15,:)); EleNamesExt{41,1}= 'FC3'; %between F3 & C3
ElePosExt(42,:) = 0.5*(ElePos(4,:)+ElePos(5,:)); EleNamesExt{42,1}= 'FC4'; %between F4 & C4
ElePosExt(43,:) = 0.5*(ElePos(15,:)+ElePos(16,:)); EleNamesExt{43,1}= 'CP3'; %between C3 & P3
ElePosExt(44,:) = 0.5*(ElePos(5,:)+ElePos(6,:)); EleNamesExt{44,1}= 'CP4'; %between C4 & P4
ElePosExt(45,:) = 0.5*(ElePos(16,:)+ElePos(13,:)); EleNamesExt{45,1}= 'PO3'; %between P3 & O1
ElePosExt(46,:) = 0.5*(ElePos(6,:)+ElePos(8,:)); EleNamesExt{46,1}= 'PO4'; %between P4 & O2
ElePosExt(47,:) = 0.5*(ElePos(18,:)+ElePos(14,:)); EleNamesExt{47,1}= 'F5'; %between F7 & F3
ElePosExt(48,:) = 0.5*(ElePos(2,:)+ElePos(4,:)); EleNamesExt{48,1}= 'F6'; %between F8 & F4
ElePosExt(49,:) = 0.5*(ElePos(19,:)+ElePos(15,:)); EleNamesExt{49,1}= 'C5'; %between T7 & C3
ElePosExt(50,:) = 0.5*(ElePos(1,:)+ElePos(5,:)); EleNamesExt{50,1}= 'C6'; %between T8 & C4
ElePosExt(51,:) = 0.5*(ElePos(17,:)+ElePos(16,:)); EleNamesExt{51,1}= 'P5'; %between P7 & P3
ElePosExt(52,:) = 0.5*(ElePos(3,:)+ElePos(6,:)); EleNamesExt{52,1}= 'P6'; %between P8 & P4
ElePosExt(53,:) = 0.5*(ElePos(12,:)+ElePos(18,:)); EleNamesExt{53,1}= 'AF7'; %between Fp1 & F7
ElePosExt(54,:) = 0.5*(ElePos(7,:)+ElePos(2,:)); EleNamesExt{54,1}= 'AF8'; %between Fp2 & F8
ElePosExt(55,:) = 0.5*(ElePos(18,:)+ElePos(19,:)); EleNamesExt{55,1}= 'FT7'; %between F7 & T7 
ElePosExt(56,:) = 0.5*(ElePos(2,:)+ElePos(1,:)); EleNamesExt{56,1}= 'FT8'; %between F8 & T8
ElePosExt(57,:) = 0.5*(ElePos(19,:)+ElePos(17,:)); EleNamesExt{57,1}= 'TP7'; %between T7 & P7
ElePosExt(58,:) = 0.5*(ElePos(1,:)+ElePos(3,:)); EleNamesExt{58,1}= 'TP8'; %between T8 & P8
ElePosExt(59,:) = 0.5*(ElePos(17,:)+ElePos(13,:)); EleNamesExt{59,1}= 'PO7'; %between P7 & O1
ElePosExt(60,:) = 0.5*(ElePos(3,:)+ElePos(8,:)); EleNamesExt{60,1}= 'PO8'; %between P8 & O2
ElePosExt(61,:) = 0.5*(ElePos(12,:)+ElePos(7,:)); EleNamesExt{61,1}= 'Fpz'; %between Fp1 & Fp2
ElePosExt(62,:) = 0.333*(ElePos(12,:)+ElePos(7,:)+ElePos(9,:)); EleNamesExt{62,1}= 'AFz'; %between  Fp1,Fp2,Fz
ElePosExt(63,:) = 0.5*(ElePos(10,:)+ElePos(11,:)); EleNamesExt{63,1}= 'CPz'; %between Cz & Pz
ElePosExt(64,:) = 0.333*(ElePos(13,:)+ElePos(14,:)+ElePos(11,:)); EleNamesExt{64,1}= 'POz'; %between O1,O2,Pz

%beyond TP7: take TP7 = 0.5*(TP9+CP5), i.e. TP9 = 2*TP7-CP5
ElePosExt(29,:) = 2*ElePosExt(57,:)-ElePosExt(27,:); EleNamesExt{29,1}= 'TP9'; 
%beyond TP7: take TP8 = 0.5*(TP10+CP6), i.e. TP10 = 2*TP8-CP6
ElePosExt(30,:) = 2*ElePosExt(58,:)-ElePosExt(28,:); EleNamesExt{30,1}= 'TP10'; 

end
