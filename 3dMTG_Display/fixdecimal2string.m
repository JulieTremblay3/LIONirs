function strnb = fixdecimal2string(nb,leftdecimal,rightdecimal)
leftdecimal;
if rightdecimal==0
    nbleftdecimal = round(nb);
else
    nbleftdecimal = fix(nb);
end

if nbleftdecimal < 0
    nbleft =sprintf(['-%0',num2str(leftdecimal),'.0f'],abs(nbleftdecimal));
else
    nbleft =sprintf(['+%0',num2str(leftdecimal),'.0f'],abs(nbleftdecimal));
end

nbrightdecimal = nb-nbleftdecimal;
nbright = sprintf(['%0',num2str(rightdecimal),'.',num2str(rightdecimal),'f'],abs(nbrightdecimal));

strnb = [nbleft,nbright(2:end)];