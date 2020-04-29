%Recoit un nom de détecteur A et le transfert en numéro de détecteur de
%Homer
function num = StrBoxy2SDDet(str)
if strcmp(str, 'D01');
    num = 1;
elseif strcmp(str, 'D02');
    num = 2;
elseif strcmp(str, 'D03');
    num = 3;
elseif strcmp(str, 'D04');
    num = 4;
elseif strcmp(str, 'D05');
    num = 5;
elseif strcmp(str, 'D06');
    num = 6;
elseif strcmp(str, 'D07');
    num = 7;
elseif strcmp(str, 'D08');
   num = 8;
elseif strcmp(str, 'D09');
    num = 9;
elseif strcmp(str, 'D10');
    num = 10;
elseif strcmp(str, 'D11');
    num = 11;
elseif strcmp(str, 'D12');
    num = 12;
elseif strcmp(str, 'D13');
    num = 13;
elseif strcmp(str, 'D14');
    num = 14;
elseif strcmp(str, 'D15');
    num = 15;
elseif strcmp(str, 'D16')
    num = 16;
elseif strcmp(str, 'D17')
    num=17;
elseif strcmp(str, 'D18')
    num=18;
elseif strcmp(str, 'D19')
    num=19;   
elseif strcmp(str, 'D20')
    num=20;   
elseif strcmp(str, 'D21') 
    num=21; 
elseif strcmp(str, 'D22') 
    num=22;    
elseif strcmp(str, 'D23') 
    num=23;    
elseif strcmp(str, 'D24') 
    num=24;    
elseif strcmp(str, 'D25') 
    num=25;
elseif strcmp(str, 'D26') 
    num=26;
elseif strcmp(str, 'D27')
    num=27;
elseif strcmp(str, 'D28')
    num=28;
elseif strcmp(str, 'D29') 
    num=29;
elseif strcmp(str, 'D30')
    num=30;
elseif strcmp(str, 'D31')
    num=31;
elseif strcmp(str, 'D32')
    num=32; 
else    
    num = 0; 
end

   