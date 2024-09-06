function TrimImage()
%Load tif file and trim all empty field around. 

[file,path] = uigetfile('.tif','multiselect','on')
if ~iscell(file)
    file = {file};
end

for i = 1:numel(file)
    [X] = imread([path,file{i}],'tif'); 
    X2d = sum(X,3);
    idx = find(765~=sum(X2d,1)./size(X2d,1));
    idy = find(765~=sum(X2d,2)./size(X2d,2));
    Total = X(idy,idx,:);
   imwrite(Total,[path,'tr', file{i},'.tif'],'tif');
end

