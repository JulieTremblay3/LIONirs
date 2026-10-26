%CAS OHR_L > OHR_R
DatGrup1 = [0.2620 %OHR_R
-0.2382
0.1690
-0.1588
-0.4561
0.3160
-0.0104
0.2612
0.1329
0.2962
0.1034
0.6396
0.5144
0.1001
0.1107
-0.1967
0.5798];

DatGrup2 = [2.2952 % OHR_L 
0.1718
1.4997
1.4920
1.2108
0.3213
0.0752
0.0535
0.0823
0.2820
0.7002
0.1169
0.2026
0.1887
0.1333
-0.3808
0.5604];



%CAS OHR_L > OHR_R 
disp('OHR_L > OHR_R')
INDEP = 0; %- DEPENDENT paired ttest
 ESTAD = 1;
 NPERM =10000;
[FSupSup,FSupDeriv,FSupTime,FUniv,T0ij] = TestPermut2Grupos(ESTAD,INDEP,DatGrup2,DatGrup1,NPERM);
 d = computeCohen_d(DatGrup2,DatGrup1);
[h,p,ci,stats] = ttest(DatGrup2,DatGrup1)
% Pval= 1-FUniv
disp(['Permutation ' num2str(NPERM),' PAIRED ttest p= ', num2str(FUniv),' t=', num2str(stats.tstat),' CohenD=',num2str(d)])


[h,p,ci,stats] =ttest(DatGrup2,DatGrup1 ,  'tail' , 'right' );
disp(['PAIRED ttest p= ' num2str(p),' t=', num2str(stats.tstat),' CohenD=',num2str(d)])


%example OHL_R
DatGrup1 = [0.4377 
-0.1061
0.3066
-0.4523
0.3799
0.0964
-0.1166
0.0773
0.1676
0.0530
0.0836
0.4019
0.4646
0.0203
0.1188
0.4261
0.0912
]
%example OHL_L

DatGrup2 = [1.0142
0.1004
-0.9048
0.5733
-1.4889
0.1619
-0.1199
0.0098
-0.0888
0.6142
-0.0742
-0.0040
0.0682
-0.0805
0.0650
-0.1553
-0.0979]



disp('OHL_R > OHL_L')
INDEP = 0; %- DEPENDENT paired ttest
 ESTAD = 1;
 NPERM =10000;
[FSupSup,FSupDeriv,FSupTime,FUniv,T0ij] = TestPermut2Grupos(ESTAD,INDEP,DatGrup1,DatGrup2,NPERM);
 d = computeCohen_d(DatGrup1,DatGrup2);

% Pval= 1-FUniv
disp(['Permutation ' num2str(NPERM),' PAIRED ttest p= ', num2str(FUniv),' t=', num2str(stats.tstat),' CohenD=',num2str(d)])
[h,p,ci,stats] =ttest(DatGrup1 ,DatGrup2,  'tail' , 'right' );
disp(['PAIRED ttest p= ' num2str(p),' t=', num2str(stats.tstat),' CohenD=',num2str(d)])