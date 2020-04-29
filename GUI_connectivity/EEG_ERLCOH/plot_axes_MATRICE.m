function handles=plot_axes_MATRICE(varargin)
%AFFICHAGE DES MATRICES DE CONNECTIVITE PAR ZONE DANS L'INTERFACE
% varargin1 handles
% varargin2 newfigure

handles = varargin{1};
if numel(varargin)>1; newfigure = varargin{2};else; newfigure = 0;end
if newfigure ==0
    axes(handles.axes_Matrice);cla;hold on
else
    figure;hold on
    colormap(jet);
end
axis ij
xlszoneorder = get(handles.edit_xlszoneorder,'string');
if isempty(xlszoneorder)
    return
end
idsubject = get(handles.popupmenu_subject,'value');
DATA = get(handles.GUI_ERLCOH_VIEW_Kmean,'UserData');

[idlabelall,idzone,idlist]=fmatorderdisplay(xlszoneorder,DATA{idsubject}.Args.ZoneList);
[option,idhalf,idoption] = flinkelelist2matrix_review(DATA{idsubject}.Args.ZoneList);
cmin =str2num(get(handles.edit_cmin,'string'));
cmax =str2num(get(handles.edit_cmax,'string'));
ctr = str2num(get(handles.edit_treshold,'string'));

nlink = size(DATA{idsubject}.MATCORR,2);
nlayer = numel(DATA{idsubject}.Args.layer);
ntime = numel(DATA{idsubject}.Args.time);
nele = numel(DATA{idsubject}.Args.ZoneList);
A =  reshape(DATA{idsubject}.MATCORR,nlayer,ntime,nlink);% DATA = get(handles.GUI_LookMat,'UserData');
if get(handles.popup_modedisplay,'value')==1; %time freq window
    twindow  = str2num(get(handles.edit_definetime,'string'));
    lwindow = str2num(get(handles.edit_frequency,'string'));
    idtwindow = find(DATA{idsubject}.Args.time>=twindow(1) & DATA{idsubject}.Args.time<=twindow(end));
    idlwindow = find(DATA{idsubject}.Args.layer>=lwindow(1) & DATA{idsubject}.Args.layer<=lwindow(end));
    if get(handles.popup_displayoption,'value')==1
        Awindow = squeeze(nanmean(nanmean(A(idlwindow ,idtwindow,:),1),2));
        matgr1 = zeros(nele,nele);
        matgr1(idhalf)=Awindow;
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
                A =  reshape(DATA{isubject}.MATCORR,nlayer,ntime,nlink);
                Awindow = squeeze(nanmean(nanmean(A(idlwindow ,idtwindow,:),1),2));
                datgr2 =  [datgr2;Awindow'];
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
        matgr1(idhalf)=Awindow;
        matgr1 = matgr1 +flipud(rot90(matgr1));
    elseif get(handles.popup_displayoption,'value')==3
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
    
    
    
elseif   get(handles.popup_modedisplay,'value')==2 %kmean
    
    icluster = get(handles.popup_clusternb,'value');
    clusterplot = get(handles.popupmenu_clusterID,'value');
    idsubject = get(handles.popupmenu_subject,'value');
    C = DATA{idsubject}.C;
    if isempty(C)
        msgbox('No cluster for this subject')
        return
    end
    tfIDX = find(C{icluster}.IDX==clusterplot);
    if get(handles.popup_displayoption,'value')==1 %JUST PLOT SELECTED FOR THIS BLOB WINDOW SUBJECT
        nlink = size(DATA{idsubject}.MATCORR,2);
        nlayer = numel(DATA{idsubject}.Args.layer);
        ntime = numel(DATA{idsubject}.Args.time);
        A =  DATA{idsubject}.MATCORR;
        Awindow = nanmean(A( tfIDX,:),1)
        mask = abs(Awindow)>ctr ;
        matgr1 = zeros(nele,nele);
        matgr1(idhalf)=Awindow.*mask  ;
        matgr1 = matgr1 +flipud(rot90(matgr1));
        imagesc(matgr1(idlist,idlist));
        caxis([cmin, cmax]);
        
    elseif get(handles.popup_displayoption,'value')==2 %PERMUTATION
        %G1 vs G2
        1
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
        
        Awindow =  nanmean(datgr1)- nanmean(datgr2);
        ptr = str2num(get(handles.edit_ptr,'string'));
        maskp = FUniv<ptr  ;
        Awindow = Awindow.*maskp;
        matgr1 = zeros(nele,nele);
        matgr1(idhalf)=Awindow ;
        matgr1 = matgr1 +flipud(rot90(matgr1));
    elseif get(handles.popup_displayoption,'value')==3 %
          nlink = size(DATA{idsubject}.MATCORR,2);
        nlayer = numel(DATA{idsubject}.Args.layer);
        ntime = numel(DATA{idsubject}.Args.time);
         A =  DATA{idsubject}.MATCORR; 
        Awindow = nanmean(A( tfIDX,:),1)
        matgr1 = zeros(nele,nele);
        matgr1(idhalf)=  Awindow;
        matgr1 = matgr1 +flipud(rot90(matgr1));
        [matgr1,idzone,idlist,idlabelall] =createROImat(matgr1,idzone,idlist,idlabelall)
   
    elseif get(handles.popup_displayoption,'value')==4
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
if  get(handles.checkbox_FLIPSIGN,'value')
    matgr1= -matgr1;
    colorMap = flipud(jet(100));
else
    colorMap = jet(100);
end
colormap( colorMap);
mask = matgr1>ctr;
matgr1 = matgr1.*mask;
imagesc(matgr1(idlist,idlist));
caxis([cmin, cmax]);
