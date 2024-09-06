%Recoit un nom de détecteur A et le transfert en numéro de détecteur de
%Homer
function num = StrBoxy2SDDet_ISS(str)
if strcmp(str,'A');
    num = 1;
elseif strcmp(str,'B');
    num = 2;
elseif strcmp(str, 'C');
    num = 3;
elseif strcmp(str,'D');
    num = 4;
elseif strcmp(str, 'E');
    num = 5;
elseif strcmp(str, 'F');
    num = 6;
elseif strcmp(str, 'G');
    num = 7;
elseif strcmp(str, 'H');
   num = 8;
elseif strcmp(str, 'I');
    num = 9;
elseif strcmp(str, 'J');
    num = 10;
elseif strcmp(str, 'K');
    num = 11;
elseif strcmp(str, 'L');
    num = 12;
elseif strcmp(str, 'M');
    num = 13;
elseif strcmp(str, 'N');
    num = 14;
elseif strcmp(str, 'O');
    num = 15;
elseif strcmp(str, 'P')
    num = 16;
elseif strcmp(str, 'Q')
    num=17;
elseif strcmp(str, 'R')
    num=18;
elseif strcmp(str, 'S')
    num=19;   
elseif strcmp(str, 'T')
    num=20;   
elseif strcmp(str, 'U') 
    num=21;          
elseif strcmp(str, 'V') 
    num=23;    
elseif strcmp(str, 'W') 
    num=24;    
elseif strcmp(str, 'X') 
    num=25;
elseif strcmp(str, 'Y') 
    num=26;
elseif strcmp(str, 'Z')
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

   