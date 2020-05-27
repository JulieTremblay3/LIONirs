function out = nirs_run_E_exportcomponent_list(job)
% Look in the list of selected factor to export the value
% Channel list determine the order of the channel to export
% missing channel in the subject or in the component list will be put as missing value. 

figureflag = 0;

 [~,~,ext] =fileparts(job.f_component_list{1});
    if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
       [data, text, rawData] = xlsread(job.f_component_list{1});
    elseif strcmp(ext,'.txt')   
        [data, text, rawData] = readtxtfile_asxlsread(job.f_component_list{1});
    end
    
    for icol=1:size(rawData,2)  
        if strcmp(upper(deblank(rawData{1,icol})),deblank(upper('NIRS.mat folder')))
            id.NIRSDtp = icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('Type'))
            id.typeDtp =  icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('Label'))
            id.labelDtp =  icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('Channel List'))
            id.ListDtp  =  icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('Zone List'))
            %id.labelDtp =  icol; todo
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('Name'))
            id.labelout=  icol;
        end      
    end
NIRSDtp = rawData(2:end,id.NIRSDtp);
typeDtp = rawData(2:end,id.typeDtp);
COMPDtp = rawData(2:end,id.labelDtp); 
ListDtp = rawData(2:end,id.ListDtp);
labeloutDtp =  rawData(2:end,id.labelout);
[pathoutlist, namelist, ext] = fileparts(job.f_component_list{1});
removeid = job.v_component_list_remove;
allsubject = [];
allsubjectHBO = [];
allsubjectHBR = [];
alllabel = [];
allt  = [];
nameall = [];
ncomp = 1;
Tglmall = [];
Rall = [];
labelall = [{'Detector'},{'Source'} ];
for filenb=1:size(NIRSDtp,1)
    try
    srsfile = [NIRSDtp{filenb},filesep,'NIRS.mat'];
    alltb =[]; 
    allb = [];
    allbHbO = [];
    allbHbR = [];
    Tglm = [];
    Rval = [];
    load(fullfile(NIRSDtp{filenb},'NIRS.mat'));
    load(fullfile(NIRSDtp{filenb},'SelectedFactors.mat'));
    disp(['Load ',NIRSDtp{filenb},filesep,'NIRS.mat'])
    try
        labelout =  labeloutDtp{filenb};
    catch
        labelout=  COMPDtp{filenb};
    end
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    label = COMPDtp{filenb}; 
    type = typeDtp{filenb};
    componenthemo =  'HbO'; %HemoDtp{filenb}
    NC = NIRS.Cf.H.C.N;
    for icomp=1:numel(PARCOMP) 
     PARCOMP(icomp).label;
    end
    scomp = ncomp;   
    disp(['Look for  ', upper(label)])
    for icomp = 1:numel(PARCOMP)
         if ~isempty(findstr(upper(PARCOMP(icomp).label),upper(label)))
            if strcmp(PARCOMP(icomp).type,type)            
              if strcmp(upper(PARCOMP(icomp).type),'PARAFAC')

                     listgood=[PARCOMP(icomp).listgood;PARCOMP(icomp).listgood+NC/2];
                     ChannelListfile = fullfile(pathoutlist,ListDtp{filenb});                     
                    [listHBOch, listHBRch, listnameHbO, listnameHbR , zonelist]= findchinChannelList(NIRS, ChannelListfile,listgood);
                                          
                    % HBO parafac 
                      k = PARCOMP(icomp).ComponentToKeep;
                     topo=  PARCOMP(icomp).FacB(:,k) * PARCOMP(icomp).FacC(1,k);
                      A = NaN(numel( listHBOch),1);
                      idok = find(~isnan(listHBOch));
                      A(idok,1)=topo(listHBOch(idok));
                     alllabel = [alllabel,{PARCOMP(icomp).label}];                  
                     allbHbO = [allbHbO,A];
                     save(fullfile(pathoutlist,['TopoHbO',labelout,PARCOMP(icomp).label,'.mat']),'A' ,'zonelist','srsfile' );
                     clear A
                     %HBR parafac
                     topo=  PARCOMP(icomp).FacB(:,k) * PARCOMP(icomp).FacC(2,k);
                     A = NaN(numel( listHBOch),1);
                      idok = find(~isnan(listHBOch));
                      A(idok,1)=topo(listHBOch(idok));
                     alllabel = [alllabel,{PARCOMP(icomp).label}];                  
                     allbHbR = [allbHbR,A];    
                     save(fullfile(pathoutlist,['TopoHbR',labelout,PARCOMP(icomp).label,'.mat']),'A','zonelist','srsfile' );
                     clear A
                    alllabel = [alllabel,{PARCOMP(icomp).label}];
                    ncomp = ncomp +1
              elseif strcmp(upper(PARCOMP(icomp).type),'GLM')
                     listgood=PARCOMP(icomp).listgood;
                    
                       [pathtxt,filetxt,ext]= fileparts(ListDtp{filenb});
                        if isempty(pathtxt)
                            ChannelListfile = fullfile(pathoutlist,ListDtp{filenb});
                        else
                              ChannelListfile = ListDtp{filenb};
                        end
                    [listHBOch, listHBRch, listnameHbO, listnameHbR, zonelist]= findchinChannelList(NIRS, ChannelListfile,listgood);
                     topo=PARCOMP(icomp).topo;
                     labelall = [labelall,{PARCOMP(icomp).label}];
  
                      A= nan(size(listHBOch,1),1);
                      idok = find(~isnan(listHBOch));
                      A(idok,1)=topo(listHBOch(idok));
                      alllabel = [alllabel,{PARCOMP(icomp).label}];                  
                      allbHbO = [allbHbO,A];
                      save(fullfile(pathoutlist,['TopoHbO',labelout,PARCOMP(icomp).label,'.mat']),'A' ,'zonelist','srsfile' );
                      clear A
                   
                      A= nan(size(listHBRch,1),1);
                      idok = find(~isnan(listHBRch));
                      A(idok,1)=topo(listHBRch(idok));               
                      allbHbR = [allbHbR,A];
                      save(fullfile(pathoutlist,['TopoHbR',labelout,PARCOMP(icomp).label,'.mat']),'A' ,'zonelist','srsfile' );
                      clear A        
                      ncomp = ncomp +1;
              elseif strcmp(upper(PARCOMP(icomp).type),'AVG')
                  listgood=PARCOMP(icomp).listgood; %check canaux pour séparer.... 
                     listgood=PARCOMP(icomp).listgood;
                     ChannelListfile = fullfile(pathoutlist,ListDtp{filenb});                     
                    [listHBOch, listHBRch, listnameHbO, listnameHbR, zonelist]= findchinChannelList(NIRS, ChannelListfile,listgood);
                    
                   alllabel = [alllabel,{PARCOMP(icomp).label}];

                      topo=PARCOMP(icomp).topo;              
                      A= nan(size(listHBOch,1),1);
                      idok = find(~isnan(listHBOch));
                      A(idok,1)=topo(listHBOch(idok));
                      alllabel = [alllabel,{PARCOMP(icomp).label}];                  
                      allbHbO = [allbHbO,A];
                      save(fullfile(pathoutlist,['TopoHbO',labelout,PARCOMP(icomp).label,'.mat']),'A' ,'zonelist','srsfile' );
                      clear A
                   
                      A= nan(size(listHBRch,1),1);
                      idok = find(~isnan(listHBRch));
                      A(idok,1)=topo(listHBRch(idok));               
                      allbHbR = [allbHbR,A];
                      save(fullfile(pathoutlist,['TopoHbR',labelout,PARCOMP(icomp).label,'.mat']),'A','zonelist','srsfile' );
                      clear A            
                      ncomp = ncomp +1;
              end             
            end            
        end
    end
    idkeep = [1:size(allb,2)];
   tmp= scomp:ncomp;
   idremovebloc = 0;
   id = find(removeid>=scomp & removeid<ncomp);
   if isempty(id)
       idremovebloc = 0;
   else
       idremovebloc = (removeid(id) - scomp+1);
       idkeep(idremovebloc) = [];
   end


    Amean = nanmean(allbHbO,2);
    allsubjectHBO =   [allsubjectHBO,allbHbO];
    A =[ Amean] ;
    save(fullfile(pathoutlist,['HBOmean',labelout,'.mat']),'A' ,'zonelist','srsfile' );  
     
    Amean = nanmean(allbHbR,2);
    A =[Amean];
    allsubjectHBR =  [allsubjectHBR,allbHbR];
    save(fullfile(pathoutlist,['HBRmean',labelout,'.mat']),'A' ,'zonelist','srsfile' );

    catch
        disp(['Could not export data from: '   srsfile])
    end
end

if 0
    NC = size(allsubject,1);
    idkeep = [1:size(allsubject,2)];
    if removeid==0   
    else
       idkeep(removeid) = [];
    end
    if isempty(allsubject)
        msgbox(['No component',label,' found'])
        return
    end
    Amean = nanmean(allsubject(:,idkeep),2);
    Astd = nanstd(allsubject(:,idkeep),1,2);
    Anb  = sum(~isnan(allsubject(:,idkeep)),2)-1;
    tmp = [alllabel;num2cell(allsubject)];
    if ismac
        try
            writetxtfile_asxlsfile(fullfile(pathoutlist,['Component' type,label,'.txt']),tmp);
        catch
            msgbox(['File ', fullfile(pathoutlist,['Component' type,label,'.txt']),' could not be written'])
        end
    else 
        try
            xlswrite(fullfile(pathoutlist,['Component' type,label,'.xls']),tmp);
        catch
            msgbox(['File ', fullfile(pathoutlist,['Component' type,label,'.xls']),' could not be written'])
        end
    end
    
    try
        
        tmp = num2cell(Tglmall);
        if ismac
            writetxtfile_asxlsfile(fullfile(pathoutlist,['Component Tval' type,label,'.txt']),tmp);
        else
            xlswrite(fullfile(pathoutlist,['Component Tval' type,label,'.xls']),tmp);
        end 
    catch
    end
    try
       tmp = [chlist{1}, chlist{2},num2cell(Rall)];
       tmp = [labelall;tmp ];
       if ismac
           writetxtfile_asxlsfile(fullfile(pathoutlist,['Component R2' type,label,'.txt']),tmp);
       else
           xlswrite(fullfile(pathoutlist,['Component R2' type,label,'.xls']),tmp);
       end
    catch
    end
    A = Amean;
    save(fullfile(pathoutlist,[labelout,'mean',type,'.mat']),'A' );

    A = Astd;
    save(fullfile(pathoutlist,[labelout,'std',label,type,'.mat']),'A' );
    nameall = [nameall, {[type,label]}];
  %  if strcmp(upper(componenthemo),'HBO')
    for ich = 1:size(allsubject,1)
        
        try
        tmp = allsubject(ich ,idkeep);
        idbad = find(isnan(tmp)|isinf(tmp));
        if isempty(idbad)
            tmp(idbad)=[];
        end
        [h,p,ci,stats] = ttest(tmp);
        T(ich) = stats.tstat;
        catch
            T(ich) = 0
        end
    end
    idbad = find(isnan(T)|isinf(T));
    T(idbad)=0;
    A = T';
    save(fullfile(pathoutlist,[labelout,'Tall','.mat']),'A' );
   
figure;
subplot(2,2,2);

imagesc(allsubject)
[val, id]=sort(nanmean(allsubject,1),'descend')
title('HbO component')
ylabel('Channels')
xlabel('Blocs')


%plot(nanmean(all(1:end/2,id),1),'x')

fs = NIRS.Cf.dev.fs;
time=(1/fs:1/fs:1/fs*size(allt,1));
subplot(2,2,1);hold on
for i=1:size(allt,2)
    if sum(i==removeid)
        plot(allt(:,i),'displayname',num2str(i),'color', 'r')
    else
        plot(allt(:,i),'displayname',num2str(i))
    end
end
title('All event, please write index to remove')
alltbad = allt(:,idkeep); 
subplot(2,2,3);
plot(time,mean(alltbad,2))

title('MEAN parafac time component')

subplot(2,2,4);
imagesc(Tglmall)
title(NIRSDtp{1})
    
end
out.NIRSmat = {fullfile(NIRSDtp{1},'NIRS.mat')};