function handles=plot_axes_linkConnectogramm(varargin)
%plot_axes_linkConnectogramm(handles, newfigure) 

%plot_axes_linkConnectogramm(handles, newfigure) 
handles = varargin{1};
if numel(varargin)>1; newfigure = varargin{2};else; newfigure = 0;end 
if newfigure ==0    
    axes(handles.axes_viewlink);
    cla;
    hold on;
else
    figure;hold on;
end
 
  fileorderconnectogramme  = get(handles.edit_linkSettingmapConnectogram,'string');
    
  [filepath,name,ext] = fileparts( fileorderconnectogramme); 
  try
    if strcmp(ext,'.txt')
         [num,txt,rawzonecolor] = readtxtfile_asxlsread(fileorderconnectogramme);     
         if size(rawzonecolor,2)>2
            listselected= rawzonecolor(2:end,1);
         else
            listselected= rawzonecolor(1:end,1);      
         end
         
    elseif strcmp(ext,'.xls')|strcmp(ext,'.xlsx')
       %  [num,txt,raw] = xlsread(fileorderconnectogramme);
           [num,txt,rawzonecolor] = xlsread(fileorderconnectogramme);
           listselected= rawzonecolor(2:end,1);
    elseif 0  %just list without color    
        fid = fopen( fileorderconnectogramme,'r');
        id = 1;
        get(handles.listbox_selectedzone,'string');
        while ~feof(fid)
            tline = fgetl(fid);
            listselected{id} = tline;
            id = id+1;
        end
            fclose(fid);

    
    end
     if size(rawzonecolor,2)>2
        maskcolor = 1;
     else
        maskcolor = 0; %could not perform mask color options
     end
  catch
       listselected = get(handles.listbox_selectedzone,'string');
       maskcolor = 0;
       disp('Advice: enter a list order for connectogram');
  end
  
DATA = get(handles.GUI_LookMat,'UserData');
id = get(handles.popup_listsujet, 'value');
MAT = DATA{id}.MAT;
if get(handles.radio_fisher,'value');
    MAT =1/2*(log((1+MAT )./(1-MAT )));
end


List=   strvcat(DATA{id}.ZoneList);
ML =  DATA{id}.zone.ml;
idzone = [];
idlist = [];
idlabelall = [];
%POUR IDENTIFICATION DE TOUTE LES ZONES
for izone = 1:numel( DATA{id}.zone.plotLst)
    chzone =  DATA{id}.zone.plotLst{izone};
    for ichzone = 1:numel(chzone)
        ich = chzone(ichzone);
        strDet = SDDet2strboxy_ISS(ML(ich,2));
        strSrs = SDPairs2strboxy_ISS(ML(ich,1));
        idch = strmatch([strDet, ' ',strSrs ],List,'exact');    
        idlist = [idlist, idch];
        if ichzone==1
            idzone =[idzone, izone];
        else
            idzone =[idzone, 0];
        end
    end
    idlabelall = [idlabelall, {[ DATA{id}.zone.label{izone}]}];
end
   try
    listok = listselected;
    catch
        disp('Please enter a list for the connectogram')
    end
    idlist = [];
    idlabel=[];
    idzone =[];
    idzonecolor = [];
    for ilistzone = 1:numel(listok)
    for izone = 1:numel(DATA{id}.zone.plotLst)
        chzone = DATA{id}.zone.plotLst{izone};
        labelzone = DATA{id}.zone.label{izone};
         x = strmatch({labelzone} , {listok{ilistzone}}, 'exact');
         if ~isempty(x)
%             for ichzone = 1:numel(chzone)
%                 ich = chzone(ichzone);
%                 if strcmp(DATA{id}.System,'ISS')
%                     strDet = SDDet2strboxy_ISS(ML(ich,2));
%                     strSrs = SDPairs2strboxy_ISS(ML(ich,1));
%                     idch = strmatch([strDet, ' ',strSrs ],List,'exact');
%                 end
                 idch=DATA{id}.zone.chMAT{izone};
                 idlist = [idlist, idch];
               % if ichzone==1
               %     idzone =[idzone, izone];
               % else
                    idzone =[idzone,izone, zeros(1,numel(idch)-1)];
                    idzonecolor = [idzonecolor, ilistzone.*ones(1,numel(idch))];
                    %  end
%             end
            idlabel = [idlabel, {[DATA{id}.zone.label{izone}, sprintf('_%03.0f',ilistzone)]}];
         end
    end
    end
    idline = [find(idzone)-0.5,numel(idzone)+0.5];
%FIN
% Create custom colormap
if get(handles.radio_negativemap,'value')
    MAT = -MAT;
else
   MAT = MAT;
end
if get(handles.popup_ConnectogramColor,'value')==1 %jet colormap list only no mask need
    if get(handles.radio_negativemap,'value')
        colorMap = flipud(jet(100));
    else
        colorMap = jet(100);
    end
    cmin = str2num(get(handles.edit_cmin,'string'));
    cmax = str2num(get(handles.edit_cmax,'string'));
    cstep = (cmax-cmin)/100;
    cf=cmin:cstep:cmax-cstep;
    for i=1:size(MAT,1)
        for j = 1:size(MAT,2)
        colorMatrix(i,j) = sum(cf<MAT(i,j));
        end
    end
     colorMatrix = colorMatrix(idlist,idlist);
elseif get(handles.popup_ConnectogramColor,'value')==2 %jet colormap mask/color needed  
    if maskcolor
    cmin = str2num(get(handles.edit_cmin,'string'));
    cmax = str2num(get(handles.edit_cmax,'string'));
    cstep = (cmax-cmin)/100;
    cf=cmin:cstep:cmax-cstep;
    for i=1:size(MAT,1)
        for j = 1:size(MAT,2)
        colorMatrix(i,j) = sum(cf<MAT(i,j));
        end
    end
     colorMatrix = colorMatrix(idlist,idlist);
     
     %use definition xls color to find all not define ROI combinaision
    colorMatrixmask= ones(numel(idlist));
    idMap = 1;
    colorMap = zeros((size(rawzonecolor,1)-1)^2,3) ;
    for i=2:size(rawzonecolor,1)
       for j=2:size(rawzonecolor,1)
           try
                colorMap(idMap,:) = str2num(rawzonecolor{i,j})./255;    
           catch
               colorMap(idMap,:) = nan;
           end
        i_id = find(idzonecolor==(i-1));
        j_id = find(idzonecolor==(j-1));
        if isnan(colorMap(idMap,:))
            colorMatrix(i_id,j_id) = nan;
        end
        idMap = idMap+1;
        end
    end
     if get(handles.radio_negativemap,'value')
        colorMap = flipud(jet(100));
    else
        colorMap = jet(100);
     end
    else
        msgbox('Define zone x zone matrix color or mask to use this mode')
        return
    end
   
elseif get(handles.popup_ConnectogramColor,'value')==3   %color by zone mask/color needed
    if maskcolor
    colorMatrix= ones(numel(idlist));
    colorMap = zeros((size(rawzonecolor,1)-1)^2,3); %-(size(rawzonecolor,1)-1)
    idMap = 1;
    for i=2:size(rawzonecolor,1)
       for j=2:size(rawzonecolor,1)
           try
                colorMap(idMap,:) = str2num(rawzonecolor{i,j})./255;    
           catch
               colorMap(idMap,:) = nan;
           end
        i_id = find(idzonecolor==(i-1));
        j_id = find(idzonecolor==(j-1));
        colorMatrix(i_id,j_id)=idMap;
        idMap = idMap+1;
        end
    end
    else
        msgbox('Define zone x zone matrix color or mask to use this mode')
        return
    end
   
elseif get(handles.popup_ConnectogramColor,'value')==4 %black color
    colorMap = zeros(100,3);
    colorMatrix = ones(size(MAT));
    colorMatrix = colorMatrix(idlist,idlist);
end
%figure;imagesc(colorMatrix)


colorlistline = zeros(size(MAT,1),3);

tmp = [find(idzone), numel(idzone)];
tmphald = floor(find(idzone) + (tmp(2:end) - tmp(1:end-1))/2);
idmiddle = zeros(1,numel(idzone));
idmiddle(tmphald ) = idzone(find(idzone));
% Create custom node labels
try
myLabel = cell(length(x));
catch
    return
end
%line(numel(listok))
idcolor = 1;
if get(handles.popup_connectogramlabel,'value')==1;
for i = 1:length(idzone)
    if idmiddle(i) %idzone(i)
        myLabel{i,1} =  DATA{id}.zone.label{idmiddle(i)};
        idcolor = idcolor + 1;
    else
        myLabel{i,1} =  ' '; 
    end
    colorhomemade(i,:) =  [0,0,0]; % colorlistline(idcolor,:);
end

for i = 1:length(idzone)    
    colorhomemade(i,:) =  [0,0,0];
   % idcolor = idcolor + 1
end

elseif get(handles.popup_connectogramlabel,'value')==2
for i=1:numel(idzone)
    try
     myLabel{i,1} = [sprintf('%03.0f', idlist(i)), List(idlist(i),:)];
    catch
        myLabel{i,1} = 'nd'
    end
end
end
%
%%

x = MAT(idlist,idlist);
id0 = find(isnan(x));
x(id0)=0;
id0 = find(isinf(x));
x(id0)=0;
thresh = str2num(get(handles.edit_threshold,'string'));
x(x < thresh) = 0;
x(x >=  thresh) = 1;

% %myconnectogram(x,myLabel,colorMatrix(idlist,idlist),colorMap,idlist);
%  myconnectogram(adjacencyMatrix,label,colorMatrix,colorMap,idlist)
 
 adjacencyMatrix = x;
 label = myLabel;

t = linspace(-pi+0.01,pi-0.01,(length(adjacencyMatrix)+1)); %to avoid 0 tan
factor = 1.1;
coordinate = zeros(length(adjacencyMatrix),2);
for i=1:length(adjacencyMatrix)  
  plot(cos(t(i)),sin(t(i)),'o','color','k');
  edge(i)=text(cos(t(i))*factor, sin(t(i))*factor,label{i},'fontsize',14) ; 
  if abs(t(i)) > pi/2
    edge(i).Rotation=(180*(t(i)/pi + 1));
    edge(i).HorizontalAlignment = 'right';
  else
   	edge(i).Rotation=(180*t(i)/pi+1);
   
  end
   edge(i).Color = [0 0 0];
end 
[row,col,v] = find(adjacencyMatrix);
if isempty(row)
    disp('No value for Connectogram')
    return
end
for i=1:numel(row)
neworder(i) = colorMatrix(row(i),col(i));
end
[val,idchrono] = sort(neworder); 
% 
% iddoublelist = [];
% for i = 1:numel(row)/2
%     iddouble= find([col==row(i)&row==col(i)]);
%        iddoublelist = [iddoublelist,iddouble];
% end
  
for iorder = 1:length(idchrono)
       i = idchrono(iorder);
       max(idlist);
       if isnan(colorMatrix(row(i),col(i)))
   
       else
      if  colorMatrix(row(i),col(i))==0;  newcolor = [0,0,0];   
      else;  newcolor = colorMap(colorMatrix(row(i),col(i)),:);   
      end
     if  abs(row(i) - col(i)) - length(adjacencyMatrix)/2 == 0;
         if ~isnan(newcolor)
           u = [cos(t(row(i)));sin(t(row(i)))];
            v = [cos(t(col(i)));sin(t(col(i)))];
           h =  line(...
              [u(1);v(1)],...
              [u(2);v(2)],...
              'LineWidth',4,...
              'Color', newcolor);     
         end
     else
   
           u  = [cos(t(row(i)));sin(t(row(i)))];
            v  = [cos(t(col(i)));sin(t(col(i)))];
            x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
            y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
            r  = sqrt(x0^2 + y0^2 - 1);
            thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
            thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
% theta = [linspace(max(thetaLim),pi,50),...
%                        linspace(-pi,min(thetaLim),50)].';
    u1 = round( u(1)*1000);
    v1 = round( v(1)*1000);
    if u1>=0 &&  v1>=0 
        theta = [linspace(max(thetaLim),pi,50),...
                           linspace(-pi,min(thetaLim),50)].';       
    else 
       theta = linspace(thetaLim(1),thetaLim(2)).';
    end
     %DRAW ARC OF CIRCLEnum2str( u1),',',num2str( v1),
   LINKNAME =  ['(i,j)=(',num2str(idlist(row(i))),',',num2str(idlist(col(i))),')=(',List(idlist(row(i)),:),',',List(idlist(col(i)),:) ,')'];
 
   if ~isnan(newcolor)
        h = line( r*cos(theta)+x0,...
                  r*sin(theta)+y0,...
                  'LineWidth', 4,...
                  'Color', newcolor,'displayname',LINKNAME );
              %set(gca,'visible','off')

       
     if newfigure ==0
        set(h,'uicontextmenu',handles.context_link);
     end
     end
   end
       end
end
