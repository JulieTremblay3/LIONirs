function [ElePosExt EleNamesExt] = nirs_boxy_interpolate10_10system(ElePos)
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
