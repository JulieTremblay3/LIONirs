function strFib = ISSIDnumber2Imaginc(nFib)
%CONVERTI LES NUMEROS DE SOURCE ISS
% a1b2 018001 -> E1
% a3b4 020003 -> E2
% a5b6 022005 -> E3
% a7b8 024007 -
% a9b10
% a11b12
% a13b14
% a15b16
% a2b1 
switch nFib
    case 018001
        strFib = 'S1';
    case 020003
        strFib = 'S2';
    case 022005
        strFib = 'S3'; 
    case 024007
        strFib = 'S4';
    case 026009
        strFib = 'S5';
    case 028011
        strFib = 'S6';
    case 030013
        strFib = 'S7';
    case 032015
        strFib = 'S8';
    case 017002
        strFib = 'S9';
    case 019004
        strFib = 'S10';
    case 021006
        strFib = 'S11';
    case 023008
        strFib = 'S12';
    case 025010
        strFib = 'S13';
    case 027012
        strFib = 'S14';    
    case 029014
        strFib = 'S15';   
    case 031016
        strFib = 'S16';   
    case 050033
        strFib = 'S17';  
    case 052035
        strFib = 'S18';   
    case 054037
        strFib = 'S19'; 
    case 056039
        strFib = 'S20'; 
    case 058041
        strFib = 'S21';
    case 060043
        strFib = 'S22';
    case 062045
        strFib = 'S23';
    case 064047
        strFib = 'S24';
    case 049034
        strFib = 'S25';
    case 051036
        strFib = 'S26';
    case 053038
        strFib = 'S27';
    case 055040
        strFib = 'S28';
    case 057042
        strFib = 'S29';
    case 059044
        strFib = 'S30';
    case 061046
        strFib = 'S31';
    case 063048 
        strFib = 'S32';  
    case 082065
        strFib = 'S33';
    case 084067
        strFib = 'S34';
    case 086069
        strFib = 'S35';
    case 088071
        strFib = 'S36';
    case 090073
        strFib = 'S37';
    case 092075
        strFib = 'S38';
    case 094077
        strFib = 'S39';   
    case 096079
        strFib = 'S40';
   case 081066
        strFib = 'S41';
    case 083068
        strFib = 'S42';
    case 085070
        strFib = 'S43';
    case 087072
        strFib = 'S44';
    case 089074
        strFib = 'S45';
    case 091076
        strFib = 'S46'; 
    case 093078
        strFib = 'S47';    
    case 095080
        strFib = 'S48';   
    case 114097
        strFib = 'S49';        
    case 116099
        strFib = 'S50';   
    case 118101
        strFib = 'S51'; 
    case 120103
        strFib = 'S52';
    case 122105
        strFib = 'S53';
    case 124107
        strFib = 'S54';
    case 126109
        strFib = 'S55';
    case 128111
        strFib = 'S56';        
    case 113098
        strFib = 'S57';
    case 115100
        strFib = 'S58';
    case 117102
        strFib = 'S59';
    case 119104
        strFib = 'S60';
    case 121106
        strFib = 'S61';
    case 123108
        strFib = 'S62';
    case 125110
        strFib = 'S63';
    case 127112
        strFib = 'S64';
    otherwise
        strFib = 'ND';
end