%association faite pour boxy a homer
function  num = StrBoxy2SDPairs(str)

if strcmp(str,'a1b2')
    num = 1;
elseif strcmp(str, 'a3b4') 
    num = 2;    
elseif strcmp(str, 'a5b6')
    num = 3;    
elseif  strcmp(str, 'a7b8')
    num = 4;   
elseif strcmp(str, 'a9b10') 
    num = 5;    
elseif strcmp(str, 'a11b12') 
    num = 6;   
elseif strcmp(str, 'a13b14')
    num = 7;
elseif strcmp(str, 'a15b16')
    num = 8;    
elseif strcmp(str, 'a2b1')
    num = 9;    
elseif strcmp(str, 'a4b3')
    num = 10;    
elseif strcmp(str, 'a6b5')
    num = 11;    
elseif strcmp(str, 'a8b7')
    num = 12;    
elseif strcmp(str, 'a10b9')
    num = 13;    
elseif strcmp(str, 'a12b11')
    num = 14;    
elseif strcmp(str, 'a14b13')
    num = 15;    
elseif strcmp(str, 'a16b15') 
    num = 16; 
elseif  strcmp(str, 'c1d2')
    num = 17;   
elseif strcmp(str, 'c3d4') 
    num = 18;    
elseif strcmp(str, 'c5d6')
    num = 19;    
elseif strcmp(str, 'c7d8')
    num = 20;    
elseif  strcmp(str, 'c9d10')
    num = 21;     
elseif strcmp(str, 'c11d12')
    num = 22;
elseif strcmp(str, 'c13d14')
    num = 23;
elseif strcmp(str, 'c15d16')
    num = 24;
elseif strcmp(str, 'c2d1')
    num = 25;   
elseif strcmp(str, 'c4d3')
    num = 26;    
elseif  strcmp(str, 'c6d5')
    num = 27;   
elseif  strcmp(str, 'c8d7')
    num = 28;   
elseif  strcmp(str, 'c10d9')
    num = 29;   
elseif   strcmp(str, 'c12d11')
    num = 30;  
elseif strcmp(str, 'c14d13')
    num = 31;    
elseif strcmp(str, 'c16d15')
    num = 32;    
elseif strcmp(str,'e1f2')
    num = 33;
elseif strcmp(str, 'e3f4') 
    num = 34;    
elseif strcmp(str, 'e5f6')
    num = 35;    
elseif  strcmp(str, 'e7f8')
    num = 36;   
elseif strcmp(str, 'e9f10') 
    num = 37;    
elseif strcmp(str, 'e11f12') 
    num = 38;   
elseif strcmp(str, 'e13f14')
    num = 39;
elseif strcmp(str, 'e15f16')
    num = 40;    
elseif strcmp(str, 'e2f1')
    num = 41;    
elseif strcmp(str, 'e4f3')
    num = 42;    
elseif strcmp(str, 'e6f5')
    num = 43;    
elseif strcmp(str, 'e8f7')
    num = 44;    
elseif strcmp(str, 'e10f9')
    num = 45;    
elseif strcmp(str, 'e12f11')
    num = 46;    
elseif strcmp(str, 'e14f13')
    num = 47;    
elseif strcmp(str, 'e16f15') 
    num = 48;
elseif  strcmp(str, 'g1h2')
    num = 49;   
elseif strcmp(str, 'g3h4') 
    num = 50;    
elseif strcmp(str, 'g5h6')
    num = 51;   
elseif strcmp(str, 'g7h8')
    num = 52;    
elseif  strcmp(str, 'g9h10')
    num = 53;     
elseif strcmp(str, 'g11h12')
    num = 54;
elseif strcmp(str, 'g13h14')
    num = 55;
elseif strcmp(str, 'g15h16')
    num = 56;
elseif strcmp(str, 'g2h1')
    num = 57;   
elseif strcmp(str, 'g4h3')
    num = 58;    
elseif  strcmp(str, 'g6h5')
    num = 59;   
elseif  strcmp(str, 'g8h7')
    num = 60;   
elseif  strcmp(str, 'g10h9')
    num = 61;   
elseif   strcmp(str, 'g12h11')
    num = 62;  
elseif strcmp(str, 'g14h13')
    num = 63;    
elseif strcmp(str, 'g16h15')
    num = 64;    
elseif 1
    num = 0;
end