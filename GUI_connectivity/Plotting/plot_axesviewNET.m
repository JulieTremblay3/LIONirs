function handles = plot_axesviewNET(varargin)


handles = varargin{1};
DATA = get(handles.GUI_LookMat,'UserData');

p = mfilename('fullpath');
[folder,name,ext]=fileparts(p);
img = load_nii([folder,'/MRI_BB.img']);
img = load_nii([folder,'/template.img']);
A = double(permute(img.img,[2 1 3]));
nbx = size(A,1); %de avant arrière
nby = size(A,2); %de gauche à droite
nbz = size(A,3); %de bas en haut
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
ylabel('Z','fontsize',14)
FID = [112 15 75;%LPA
    8 11 80;%RPA
    57 12 8]; %NAS
%AJUSTEMENT EN Z
zoffset = 5;

displaymap = [get(handles.radio_Axi_Down,'value'),...
    get(handles.radio_Axi_Up,'value'),...
    get(handles.radio_Sag_Left,'value'),...
    get(handles.radio_Sag_Right,'value'),...
    get(handles.radio_Cor_Front,'value'),...
    get(handles.radio_Cor_Back,'value')];
imap = 1;
haxes = zeros(numel(displaymap),1);
hfig = figure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Axial z bas à haut (view haut)
if displaymap(1)
    figurehaxes(1)=subplot(1,sum(displaymap),imap); imap = imap+1;hold on
    slice =round(nbz/2);
    [xx,yy,zz] = meshgrid((1:nbx)-round(nbx/2),(1:nby)-round(nby/2),slice-zoffset);
    surface(xx,yy,zz,squeeze(A(:,:, slice)'),...
        'FaceColor','texturemap',...
        'EdgeColor','none',...
        'CDataMapping','direct');
    map = gray(255);
    colormap(map);
    set(gca,'color',[0 0 0]);
    set(gca,'CameraTarget',[0,0,0]);
    set(gca,'CameraPosition',[0,0,800]);
    set(gca,'CameraUpVector',[1,0,0]); %(en y on regarde en haut)
    xlabel('X','fontsize',14);
    ylabel('Y','fontsize',14);
    ylabel('Z','fontsize',14);
    plotchannel(handles,DATA);
    title('Top view')

end

%Axial z bas à haut (view bas)
if displaymap(2)
    haxes(2)=subplot(1,sum(displaymap),imap);imap = imap+1;hold on
    slice =round(nbz/2);
    [xx,yy,zz] = meshgrid((1:nbx)-round(nbx/2),(1:nby)-round(nby/2),slice-zoffset);
    surface(xx,yy,zz,squeeze(A(:,:, slice)'),...
        'FaceColor','texturemap',...
        'EdgeColor','none',...
        'CDataMapping','direct');
    map = gray(255);
    colormap(map);
    set(gca,'CameraTarget',[0,0,0]);
    set(gca,'CameraPosition',[0,0,-800]);
    set(gca,'CameraUpVector',[1,0,0]); %(en y on regarde en haut)
    set(gca,'color',[0 0 0]);
    xlabel('X','fontsize',14);
    ylabel('Y','fontsize',14);
    ylabel('Z','fontsize',14);
    plotchannel(handles,DATA);
    title('Bottom view')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Sagital y droite gauche
if displaymap(3)
    haxes(3) = subplot(1,sum(displaymap),imap);imap = imap+1;
    hold on
    slice = round(nby/2);
    [xx,yy,zz] = meshgrid((1:nbx)-round(nbx/2),slice-round(nby/2),(1:nbz)-zoffset);
    surface(squeeze(xx),squeeze(yy),squeeze(zz),squeeze(A(:,slice,: )),...
        'FaceColor','texturemap',...
        'EdgeColor','none',...
        'CDataMapping','direct');
    set(gca,'CameraTarget',[0,0,0]);
    set(gca,'CameraPosition',[0,800,0]);
    set(gca,'CameraUpVector',[0,0,1]);
    set(gca,'color',[0 0 0]);
    map = gray(255);
    colormap(map);
    xlabel('X','fontsize',14);
    ylabel('Y','fontsize',14);
    ylabel('Z','fontsize',14);
    plotchannel(handles,DATA);
    title('Sagital Left ');
end

%Sagital right
if displaymap(4)
    haxes(4)= subplot(1,sum(displaymap),imap);;imap = imap+1;hold on
    slice = round(nby/2);
    [xx,yy,zz] = meshgrid((1:nbx)-round(nbx/2),slice-round(nby/2),(1:nbz)-zoffset)
    surface(squeeze(xx),squeeze(yy),squeeze(zz),squeeze(A(:,slice,: )),...
        'FaceColor','texturemap',...
        'EdgeColor','none',...
        'CDataMapping','direct');
    set(gca,'CameraTarget',[0,0,0]);
    set(gca,'CameraPosition',[0,-800,0]);
    set(gca,'CameraUpVector',[0,0,1]);
    set(gca,'color',[0 0 0]);
    map = gray(255);
    colormap(map);
    xlabel('X','fontsize',14);
    ylabel('Y','fontsize',14);
    ylabel('Z','fontsize',14);
    plotchannel(handles,DATA);
    title('Sagital Right');
end


%Coronal x (avant arrière)
if displaymap(5)
    haxes(5)= subplot(1,sum(displaymap),imap);;imap = imap+1;hold on
    slice = round(nbx/2);
    [xx,yy,zz] = meshgrid(slice-round(nbx/2),(1:nby)-round(nby/2),(1:nbz)-zoffset)
    surface(squeeze(xx),squeeze(yy),squeeze(zz),squeeze(A(slice,:,: )),...
        'FaceColor','texturemap',...
        'EdgeColor','none',...
        'CDataMapping','direct');
    set(gca,'CameraTarget',[0,0,0]);
    set(gca,'CameraPosition',[800,0,0]);
    set(gca,'CameraUpVector',[0,0,1]);
    set(gca,'color',[0 0 0]);
    map = gray(255);
    colormap(map);
    xlabel('X','fontsize',14);
    ylabel('Y','fontsize',14);
    zlabel('Z','fontsize',14);
    plotchannel(handles,DATA);
    title('Coronal Front');
end
if displaymap(6)
    haxes(6)= subplot(1,sum(displaymap),imap);        imap = imap+1;hold on
    slice = round(nbx/2);
    [xx,yy,zz] = meshgrid(slice-round(nbx/2),(1:nby)-round(nby/2),(1:nbz)-zoffset);
    surface(squeeze(xx),squeeze(yy),squeeze(zz),squeeze(A(slice,:,: )),...
        'FaceColor','texturemap',...
        'EdgeColor','none',...
        'CDataMapping','direct');
    set(gca,'CameraTarget',[0,0,0]);
    set(gca,'CameraPosition',[-800,0,0]);
    set(gca,'CameraUpVector',[0,0,1]);
    set(gca,'color',[0 0 0]);
    map = gray(255);
    colormap(map);
    xlabel('X','fontsize',14);
    ylabel('Y','fontsize',14);
    zlabel('Z','fontsize',14);
    plotchannel(handles,DATA);
    title('Coronal Back');
    
end
end

function plotchannel(handles,DATA)
isubject = get(handles.popup_listsujet,'value');
ML = DATA{isubject}.zone.ml;
pos = DATA{isubject}.zone.pos*10;
for p=1:size(pos,1)
    x = pos(p,1);
    y = pos(p,2);
    z = pos(p,3);
    h = plot3(x,y,z,'xy')
    strDet = SDDet2strboxy_ISS(ML(p,2));
    strSrs = SDPairs2strboxy_ISS(ML(p,1));
    set(h,'linewidth',5,'displayname',[strDet, ' ',strSrs ]);
end
 colormaplist =jet(200)
idzone = get(handles.listbox_subzone_selected,'string');
for idlistzone = 1:numel(idzone)
    [str1,str2]=strtok(idzone{idlistzone },'/')
    idlabelall = [];
    for izone = 1:numel(DATA{isubject}.zone.plotLst)
        idlabelall = [idlabelall, {DATA{isubject}.zone.label{izone}}];
    end
    str = {str1}
    x = strmatch(str ,idlabelall ,  'exact')
    if ~isempty(x)
        idlistx = DATA{isubject}.zone.chMAT{x};
    end
    str = {str2(2:end)};
    y = strmatch(str,idlabelall ,  'exact')
    if ~isempty(y)
        idlisty = DATA{isubject}.zone.chMAT{y};
    end
    thr = str2num(get(handles.edit_threshold,'string'));
    Mat = DATA{isubject}.MAT;
    Mat = threshold_absolute(Mat,thr);
    id = find(isnan(Mat));
    if ~isempty(id);Mat(id)=0;end;
    subMAT = Mat([idlistx,idlisty],[idlistx,idlisty]);
    submatint = DATA{isubject}.MAT([idlistx,idlisty],[idlistx,idlisty]);
    maskval = [zeros(numel(idlistx ),numel(idlistx )),zeros(numel(idlistx ),numel(idlisty ));
        ones(numel(idlisty ),numel(idlistx )),zeros(numel(idlisty ),numel(idlisty )) ];
    subMAT = ceil(subMAT.*maskval);
    subpos = pos([idlistx,idlisty],1:3)
    idlistxy = [idlistx,idlisty];
    [row,col] = find(subMAT)
    if isempty(row)
        return
    end
    hold on
    for i=1:numel(row)
        val(i) = submatint(row(i),col(i))
    end
   [val,idsort]=sort(val);

   cmax = str2num(get(handles.edit_cmax,'string'));
   cmin = str2num(get(handles.edit_cmin,'string'));
    for i=1:size(row,1)
        iorder =idsort(i);
        h = plot3([subpos(row(iorder),1); subpos(col(iorder),1)],[subpos(row(iorder),2); subpos(col(iorder),2)],[subpos(row(iorder),3); subpos(col(iorder),3)],'y' );
        if 1          
          norm_data = (submatint(row(iorder),col(iorder)) - cmin) / ( cmax -  cmin );
          id = round(norm_data*200);
          if id > 200; id=200;end
          if id <1; id=1;end
             
          set(h,'color', colormaplist(id,:));
          set(h,'linewidth', 4)
         
            strDetA = SDDet2strboxy_ISS(ML(idlistxy(row(iorder)),2));
            strSrsA = SDPairs2strboxy_ISS(ML(idlistxy(row((iorder))),1));
            strDetB = SDDet2strboxy_ISS(ML(idlistxy(col((iorder))),2));
            strSrsB = SDPairs2strboxy_ISS(ML(idlistxy(col((iorder))),1));
             set(h,'linewidth',5,'displayname',[strDetA, ' ',strSrsA , 'to',strDetB, ' ',strSrsB  ]);
      
        end
    end
    % [X,Y,Z] = adjacency_plot_und(subMAT ,pos)
    % figure;plot3(X,Y,Z);
end
end
%
%
%
%     hold on
%     idsujet = get(handles.popup_listsujet,'value')
%     ML = DATA{idsujet}.zone.ml
%     pos = DATA{idsujet}.zone.pos
%     title('view up')
%     vColor = [0 0 1]
%     for p=1:size(pos,1)
%         x = pos(p,1)*10;
%         y = pos(p,2)*10;
%         z = pos(p,3)*10;
%         h = plot3(x,y,z,'xy')
%         strDet = SDDet2strboxy_ISS(ML(p,2));
%         strSrs = SDPairs2strboxy_ISS(ML(p,1))
%         set(h,'linewidth',5,'displayname',[strDet, ' ',strSrs ])
%     end