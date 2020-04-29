function [pathout,nameout] = combinetiffileV_only(path,file)
Total = [];
for i = 1:numel(file)
    [X] = imread([path,file{i}],'tif'); 
    Total = [Total;X];
end
figure
imagesc(Total)
[nbxT,nby,z]= size(Total)

[nameout,pathout]=uiputfile('all.tif');
imwrite(Total,[pathout,nameout],'tif');

