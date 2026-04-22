function strFib = ISSIDnumber2NIRx(nFib)
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
        strFib = 'E1';
    case 020003
        strFib = 'E2';
    case 022005
        strFib = 'E3'; 
    case 024007
        strFib = 'E4';
    case 026009
        strFib = 'E5';
    case 028011
        strFib = 'E6';
    case 030013
        strFib = 'E7';
    case 032015
        strFib = 'E8';
    case 017002
        strFib = 'E9';
    case 019004
        strFib = 'E10';
    case 021006
        strFib = 'E11';
    case 023008
        strFib = 'E12';
    case 025010
        strFib = 'E13';
    case 027012
        strFib = 'E14';    
    case 029014
        strFib = 'E15';   
    case 031016
        strFib = 'E16';   
    case 050033
        strFib = 'E17';  
    case 052035
        strFib = 'E18';   
    case 054037
        strFib = 'E19'; 
    case 056039
        strFib = 'E20'; 
    case 058041
        strFib = 'E21';
    case 060043
        strFib = 'E22';
    case 062045
        strFib = 'E23';
    case 064047
        strFib = 'E24';
    case 049034
        strFib = 'E25';
    case 051036
        strFib = 'E26';
    case 053038
        strFib = 'E27';
    case 055040
        strFib = 'E28';
    case 057042
        strFib = 'E29';
    case 059044
        strFib = 'E30';
    case 061046
        strFib = 'E31';
    case 063048 
        strFib = 'E32';  
    case 082065
        strFib = 'E33';
    case 084067
        strFib = 'E34';
    case 086069
        strFib = 'E35';
    case 088071
        strFib = 'E36';
    case 090073
        strFib = 'E37';
    case 092075
        strFib = 'E38';
    case 094077
        strFib = 'E39';   
    case 096079
        strFib = 'E40';
   case 081066
        strFib = 'E41';
    case 083068
        strFib = 'E42';
    case 085070
        strFib = 'E43';
    case 087072
        strFib = 'E44';
    case 089074
        strFib = 'E45';
    case 091076
        strFib = 'E46'; 
    case 093078
        strFib = 'E47';    
    case 095080
        strFib = 'E48';   
    case 114097
        strFib = 'E49';        
    case 116099
        strFib = 'E50';   
    case 118101
        strFib = 'E51'; 
    case 120103
        strFib = 'E52';
    case 122105
        strFib = 'E53';
    case 124107
        strFib = 'E54';
    case 126109
        strFib = 'E55';
    case 128111
        strFib = 'E56';        
    case 113098
        strFib = 'E57';
    case 115100
        strFib = 'E58';
    case 117102
        strFib = 'E59';
    case 119104
        strFib = 'E60';
    case 121106
        strFib = 'E61';
    case 123108
        strFib = 'E62';
    case 125110
        strFib = 'E63';
    case 127112
        strFib = 'E64';
    otherwise
        strFib = 'ND';
end 