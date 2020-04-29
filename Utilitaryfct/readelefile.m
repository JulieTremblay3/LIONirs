function [name_ele,x,y,z]=readelefile(elefile)
fid = fopen([elefile]);
if fid == -1
    msgbox('The ele file is not valid : ',[elefile]);
    return
end
line = 0;
name_ele = [];
x = [];
y = [];
z = [];
while line ~= -1
    line =  fgetl(fid);
    if line == -1 
        break
    end
    if strcmp(strtrim(line),'')
        break
    end
    [name, rem] = strtok(line);
    name_ele = [name_ele, {name}];
    [i, rem] = strtok(rem);
    [j, rem] = strtok(rem);
    [k, rem] = strtok(rem);
    x = [x,str2num(i)];
    y = [y,str2num(j)];
    z = [z,str2num(k)];
end
fclose(fid);