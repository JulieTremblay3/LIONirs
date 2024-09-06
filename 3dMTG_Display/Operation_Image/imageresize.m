function img=imageresize(img,x,y)
%resize to equalize dimension
if ~isempty(x)
    m = size(img,1)/x;
    img2 = imresize(img,1/m);
    clear img
    if size(img2,1)< x
          img = [img2;255*ones(x-size(img2,1),size(img2,2),3)];     
    elseif size(img2,1)> x   
          img = img2(1:x,:,:);     
    else
          img=img2;
    end   
end

if ~isempty(y)
    m = size(img,2)/y;
    img2 = imresize(img,1/m);
    clear img
    if size(img2,2)<y
          img = [img2,255*ones(size(img2,1),y-size(img2,2),3)];     
    elseif size(img2,2)> y   
          img = img2(:,1:y,:);     
    else
          img=img2;
    end      
end