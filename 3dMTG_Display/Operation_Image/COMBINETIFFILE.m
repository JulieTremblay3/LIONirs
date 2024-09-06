[file,path] = uigetfile('.tif','multiselect','on')
Total = [];
for i = 1:numel(file)
    [X] = imread([path,file{i}],'tif'); 
    Total = [Total,X];
end
imwrite(Total,['all', file{1},'.tif'],'tif');

[file,path] = uigetfile('.tif','multiselect','on')
[X] = imread([path,file],'tif'); 



[file,path] = uigetfile('.tif','multiselect','on')
Total = [];
for i = 1:numel(file)
    [X] = imread([path,file{i}],'tif'); 
    Total = [Total;X];
end
imwrite(Total,'all.tif','tif');