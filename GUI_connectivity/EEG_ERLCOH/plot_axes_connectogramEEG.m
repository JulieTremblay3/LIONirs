function handles=plot_axes_connectogramEEG(varargin)
%AFFICHAGE DES MATRICES DE CONNECTIVITE PAR ZONE DANS L'INTERFACE
% varargin1 handles
% varargin2 newfigure
 
handles = varargin{1};
if numel(varargin)>1; newfigure = varargin{2};else; newfigure = 0;end 
if newfigure ==0
    axes(handles.axes_connectogram);cla;hold on
else
    figure;hold on
    colormap(jet)
end
axis ij
 xlszoneorder = get(handles.edit_XLSconnectogram,'string');
 if isempty(xlszoneorder)
    return
 end
idsubject = get(handles.popupmenu_subject,'value');
DATA = get(handles.GUI_ERLCOH_VIEW_Kmean,'UserData');
[idlabelall,idzone,idlist]=fmatorderdisplay_CONNECTOGRAM(xlszoneorder,DATA{idsubject}.Args.ZoneList);
[option,idhalf,idoption] = flinkelelist2matrix_review(DATA{idsubject}.Args.ZoneList);
cmin =str2num(get(handles.edit_cmin,'string'));
cmax =str2num(get(handles.edit_cmax,'string'));
ctr = str2num(get(handles.edit_treshold,'string'));

nlink = size(DATA{idsubject}.MATCORR,2);
nlayer = numel(DATA{idsubject}.Args.layer);
ntime = numel(DATA{idsubject}.Args.time);
nele = numel(DATA{idsubject}.Args.ZoneList);
if get(handles.popup_modedisplay,'value')==1; %time freq window
twindow  = str2num(get(handles.edit_definetime,'string'));
lwindow = str2num(get(handles.edit_frequency,'string'));
idtwindow = find(DATA{idsubject}.Args.time>=twindow(1) & DATA{idsubject}.Args.time<=twindow(end));
idlwindow = find(DATA{idsubject}.Args.layer>=lwindow(1) & DATA{idsubject}.Args.layer<=lwindow(end));
    if get(handles.popup_displayoption,'value')==1 
        A =  reshape(DATA{idsubject}.MATCORR,nlayer,ntime,nlink);
        Awindow = squeeze(nanmean(nanmean(A(idlwindow ,idtwindow,:),1),2));
         matgr1 = zeros(nele,nele);
        matgr1(idhalf)=  Awindow;
        matgr1 = matgr1 +flipud(rot90(matgr1));
    elseif get(handles.popup_displayoption,'value')==2
          datgr1 = [];
        datgr2 = [];
        for isubject = 1:numel(DATA)
            if DATA{isubject}.GR==1
                A =  reshape(DATA{isubject}.MATCORR,nlayer,ntime,nlink);
                Awindow = squeeze(nanmean(nanmean(A(idlwindow ,idtwindow,:),1),2));
                 datgr1 =  [datgr1;Awindow'];
            end
              if DATA{isubject}.GR==2
                 A = reshape(DATA{isubject}.MATCORR,nlayer,ntime,nlink);
                Awindow = squeeze(nanmean(nanmean(A(idlwindow ,idtwindow,:),1),2));
                 datgr2 =  [datgr2;Awindow'];
            end
        end
        ESTAD = 4;
        INDEP = 1;
        NPERM = 500 
        [FSupSup,FSupDeriv,FSupTime,FUniv,  Toij] = TestPermut2Grupos(ESTAD,INDEP,datgr1,datgr2 ,NPERM); 
        Awindow =  nanmean(datgr1,1)- nanmean(datgr2,1)
        ptr = str2num(get(handles.edit_ptr,'string'));
        maskp = FUniv<ptr ;
        
        Awindow = Awindow.*maskp;
         matgr1 = zeros(nele,nele);
        matgr1(idhalf)=  Awindow;
        matgr1 = matgr1 +flipud(rot90(matgr1));
    elseif get(handles.popup_displayoption,'value')==3 %kmean View subject cluster
        nlink = size(DATA{idsubject}.MATCORR,2);
        nlayer = numel(DATA{idsubject}.Args.layer);
        ntime = numel(DATA{idsubject}.Args.time);
        
        A =  reshape(DATA{idsubject}.MATCORR,nlayer,ntime,nlink);
        Awindow = squeeze(nanmean(nanmean(A(idlwindow ,idtwindow,:),1),2));
        matgr1 = zeros(nele,nele);
        matgr1(idhalf)=  Awindow;
        matgr1 = matgr1 +flipud(rot90(matgr1));
        [matgr1,idzone,idlist,idlabelall] =createROImat(matgr1,idzone,idlist,idlabelall)

    elseif get(handles.popup_displayoption,'value')==4 
         datgr1 = [];
        datgr2 = [];
        for isubject = 1:numel(DATA)
            if DATA{isubject}.GR==1
                 A =  reshape(DATA{isubject}.MATCORR,nlayer,ntime,nlink);; 
                Awindow = squeeze(nanmean(nanmean(A(idlwindow ,idtwindow,:),1),2));
                 matgr1 = zeros(nele,nele);
                 matgr1(idhalf)=  Awindow;
                 matgr1 = matgr1 +flipud(rot90(matgr1));
                  [matgr1,idzonenew,idlistnew,idlabelallnew] =createROImat(matgr1,idzone,idlist,idlabelall)
                 
                   matgr1 = reshape(matgr1,1,size(matgr1,1),size(matgr1,2));
                   datgr1 =  cat(1,datgr1,matgr1);
            end
              if DATA{isubject}.GR==2
                    A =  reshape(DATA{isubject}.MATCORR,nlayer,ntime,nlink);; 
                Awindow = squeeze(nanmean(nanmean(A(idlwindow ,idtwindow,:),1),2));
                 matgr1 = zeros(nele,nele);
                 matgr1(idhalf)=  Awindow;
                 matgr1 = matgr1 +flipud(rot90(matgr1));
                  [matgr1,idzonenew,idlistnew,idlabelallnew] =createROImat(matgr1,idzone,idlist,idlabelall)
                 matgr1 = reshape(matgr1,1,size(matgr1,1),size(matgr1,2));
                  datgr2 =  cat(1,datgr2,matgr1);
            end
        end
        ESTAD = 4;
        INDEP = 1;
        NPERM = 500 
        
             [FSupSup,FSupDeriv,FSupTime,FUniv,  Toij] = TestPermut2Grupos(ESTAD,INDEP,datgr1,datgr2 ,NPERM); 
        matgr1=  squeeze(nanmean(datgr1,1)- nanmean(datgr2,1))
        ptr = str2num(get(handles.edit_ptr,'string'));
        maskp = FUniv<ptr ;
         matgr1 = matgr1.*maskp
        idzone = idzonenew;
        idlist = idlistnew;
        idlabelall=idlabelallnew;
        
    end
elseif get(handles.popup_modedisplay,'value')==2
     icluster = get(handles.popup_clusternb,'value');
    clusterplot = get(handles.popupmenu_clusterID,'value');
    idsubject = get(handles.popupmenu_subject,'value');
    C = DATA{idsubject}.C;
    if isempty(C)
        msgbox('No cluster available for this subject')
    end
    tfIDX = find(C{icluster}.IDX==clusterplot);
    if get(handles.popup_displayoption,'value')==1 %JUST PLOT SELECTED FOR THIS BLOB WINDOW SUBJECT 
        nlink = size(DATA{idsubject}.MATCORR,2);
        nlayer = numel(DATA{idsubject}.Args.layer);
        ntime = numel(DATA{idsubject}.Args.time);
        A =  DATA{idsubject}.MATCORR; 
        Awindow = nanmean(A( tfIDX,:),1)
         matgr1 = zeros(nele,nele);
        matgr1(idhalf)=  Awindow;
        matgr1 = matgr1 +flipud(rot90(matgr1));
    elseif get(handles.popup_displayoption,'value')==2
           datgr1 = [];
        datgr2 = [];
        for isubject = 1:numel(DATA)
            if DATA{isubject}.GR==1
                 A =  DATA{isubject}.MATCORR; 
                Awindow = nanmean(A( tfIDX,:),1);
                 datgr1 =  [datgr1;Awindow];
            end
              if DATA{isubject}.GR==2
                 A =  DATA{isubject}.MATCORR; 
                Awindow = nanmean(A( tfIDX,:),1);
                 datgr2 =  [datgr2;Awindow];
            end
        end
        ESTAD = 4;
        INDEP = 1;
        NPERM = 500 
        [FSupSup,FSupDeriv,FSupTime,FUniv,  Toij] = TestPermut2Grupos(ESTAD,INDEP,datgr1,datgr2 ,NPERM); 
        Awindow =  nanmean(datgr1)- nanmean(datgr2)
        ptr = str2num(get(handles.edit_ptr,'string'));
        maskp = FUniv<ptr ;
        
        Awindow = Awindow.*maskp;
        matgr1 = zeros(nele,nele);
        matgr1(idhalf)=  Awindow;
        matgr1 = matgr1 +flipud(rot90(matgr1));
    elseif get(handles.popup_displayoption,'value')==3 %kmean View subject cluster

        nlink = size(DATA{idsubject}.MATCORR,2);
        nlayer = numel(DATA{idsubject}.Args.layer);
        ntime = numel(DATA{idsubject}.Args.time);
        A =  DATA{idsubject}.MATCORR; 
        Awindow = nanmean(A( tfIDX,:),1)
        matgr1 = zeros(nele,nele);
        matgr1(idhalf)=  Awindow;
        matgr1 = matgr1 +flipud(rot90(matgr1));
       [matgr1,idzone,idlist,idlabelall] =createROImat(matgr1,idzone,idlist,idlabelall)     
    elseif get(handles.popup_displayoption,'value')==4 %kmean View subject cluster
           datgr1 = [];
        datgr2 = [];
        for isubject = 1:numel(DATA)
            if DATA{isubject}.GR==1
                 A =  DATA{isubject}.MATCORR; 
                Awindow = nanmean(A( tfIDX,:),1);
                 matgr1 = zeros(nele,nele);
                 matgr1(idhalf)=  Awindow;
                 matgr1 = matgr1 +flipud(rot90(matgr1));
                  [matgr1,idzonenew,idlistnew,idlabelallnew] =createROImat(matgr1,idzone,idlist,idlabelall)
                 
                   matgr1 = reshape(matgr1,1,size(matgr1,1),size(matgr1,2));
                   datgr1 =  cat(1,datgr1,matgr1);
            end
              if DATA{isubject}.GR==2
                   A =  DATA{isubject}.MATCORR; 
                Awindow = nanmean(A( tfIDX,:),1);
                 matgr1 = zeros(nele,nele);
                 matgr1(idhalf)=  Awindow;
                 matgr1 = matgr1 +flipud(rot90(matgr1));
                  [matgr1,idzonenew,idlistnew,idlabelallnew] =createROImat(matgr1,idzone,idlist,idlabelall)
                 matgr1 = reshape(matgr1,1,size(matgr1,1),size(matgr1,2));
                  datgr2 =  cat(1,datgr2,matgr1);
            end
        end
        ESTAD = 4;
        INDEP = 1;
        NPERM = 500 
        
             [FSupSup,FSupDeriv,FSupTime,FUniv,  Toij] = TestPermut2Grupos(ESTAD,INDEP,datgr1,datgr2 ,NPERM); 
        matgr1=  squeeze(nanmean(datgr1,1)- nanmean(datgr2,1))
        ptr = str2num(get(handles.edit_ptr,'string'));
        maskp = FUniv<ptr ;
         matgr1 = matgr1.*maskp
        idzone = idzonenew;
        idlist = idlistnew;
        idlabelall=idlabelallnew;
     end
    

   
end

if get(handles.checkbox_FLIPSIGN,'value')
    colorMap = flipud(jet(100));
    matgr1 = -matgr1;
else
    colorMap = jet(100);
end
mask =    matgr1 >ctr ;
matgr1 = matgr1.*mask

%PLOT
try

x = zeros(numel(idlist));
%POSITIVE LINK
x(matgr1(idlist,idlist) >  ctr) = 1;
x(matgr1(idlist,idlist) <= ctr) = 0;
colorlistline =zeros(size(x,1),3);% % 
cstep = (cmax-cmin)/100;
cf=cmin:cstep:cmax-cstep;


clear idc
for i=1:size(matgr1,1)
    for j=1:size(matgr1,2)
       if 0 == sum(cf<matgr1(i,j));
           matc(i,j)= 1;
       else
           matc(i,j)= sum(cf<matgr1(i,j));
       end      
    end
end
 colorMatrix =  matc(idlist,idlist);
 adjacencyMatrix = x
 label = idlabelall


t = linspace(-pi+0.01,pi-0.01,(length(adjacencyMatrix)+1)) %to avoid 0 tan
factor = 1.1
coordinate = zeros(length(adjacencyMatrix),2)
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
       max(idlist)
      if  colorMatrix(row(i),col(i))==0;  newcolor = [0,0,0];    else;  newcolor = colorMap(colorMatrix(row(i),col(i)),:);    end
     if  abs(row(i) - col(i)) - length(adjacencyMatrix)/2 == 0
           u = [cos(t(row(i)));sin(t(row(i)))];
            v = [cos(t(col(i)));sin(t(col(i)))];
           h =  line(...
              [u(1);v(1)],...
              [u(2);v(2)],...
              'LineWidth',4,...
              'Color', newcolor);          
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
   List=  DATA{idsubject}.Args.ZoneList;
   LINKNAME =  ['(i,j)=(',num2str(idlist(row(i))),',',num2str(idlist(col(i))),')=(',List{idlist(row(i))},',',List{idlist(col(i))} ,')'];

   
    h = line( r*cos(theta)+x0,...
              r*sin(theta)+y0,...
              'LineWidth', 4,...
              'Color', newcolor,'displayname',LINKNAME );
          %set(gca,'visible','off')
    
     end
     if newfigure ==0
        %set(h,'uicontextmenu',handles.context_link);
    end
end

axis xy


catch
end
