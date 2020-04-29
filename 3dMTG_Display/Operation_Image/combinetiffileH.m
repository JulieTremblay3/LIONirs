function [nameout,pathout]=combinetiffileH(path,file,limitex,limitey,time)
Total = [];
for i = 1:numel(file)
    x =limitex(1,1):limitex(1,2);
    y =limitey(1,1):limitey(1,2); 
    [img] = imread([path,file{i}],'tif');
%     figure
%     imagesc(img)
    img=img(y,x,:);
    
    %ADD number time
    if ~isempty(time)
        h = figure    
        text(0.5,0.5,[num2str(round(time(i)))],'FontSize',50,'HorizontalAlignment','left','VerticalAlignment','bottom')
        saveas(h,'temp','tif');
        close(h)    
        imgall = imread('temp','tif');     
        nb = imgall(330:420,610:794,:); % attention y,x sur la figure pour couper le chiffre
        [nbx,nby,x] = size(nb);
        img(1:nbx,1:nby,:) = nb;
    end
%      figure
%     imagesc(img);
    Total = [Total,img];
end
[nameout,pathout]=uiputfile([path,'all',file{1},'.tif']);
imwrite(Total,[pathout,nameout],'tif');
