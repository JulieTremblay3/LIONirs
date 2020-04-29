function out = nirs_run_E_exportcomponent_zone(job)
% Look in the list of selected factor to export the value
% Channel list determine the order of the channel to export
% missing channel in the subject or in the component list will be put as missing value. 
 
figureflag = 0;
[data, text, rawData] = xlsread(job.f_component_zone{1});
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
            id.ZoneDtp  =  icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('Name'))
            id.labelout=  icol;
        end      
    end
    
    if isfield(job.f_zonetoextract_file,'f_zonetoextract_file')
        fid = fopen(job.f_zonetoextract_file.f_zonetoextract_file{1})
        ROI = textscan(fid,'%s')
        zonelabel =  ROI{1}
        fclose(fid)
    elseif isfield(job.f_zonetoextract_file,'tocheck')
    
    
    end
NIRSDtp = rawData(2:end,id.NIRSDtp);
typeDtp = rawData(2:end,id.typeDtp);
COMPDtp = rawData(2:end,id.labelDtp); 
%ListDtp = rawData(2:end,id.ListDtp);
ZoneDtp = rawData(2:end,id.ZoneDtp );
labeloutDtp =  rawData(2:end,id.labelout);
[pathoutlist, namelist, ext] = fileparts(job.f_component_zone{1});

    allbHbO = [];
    allbHbR = [];
    labelall = [{'Zone'}];
    ncomp = 1;
for filenb=1:size(NIRSDtp,1)
    srsfile = [NIRSDtp{filenb},'\NIRS.mat'];  
    load(fullfile(NIRSDtp{filenb},'NIRS.mat'));
    load(fullfile(NIRSDtp{filenb},'SelectedFactors.mat'));
    disp(['Load ',NIRSDtp{filenb},'\NIRS.mat'])
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
                     ChannelListfile = fullfile(pathoutlist,ZoneDtp{filenb});                     
                     
                    [listHBOch, listHBRch, listnameHbO, listnameHbR , zonelist]= findchinChannelList(NIRS, ChannelListfile,listgood);
                                          
                    % HBO parafac 
                      k = PARCOMP(icomp).ComponentToKeep;
                     topo=  PARCOMP(icomp).FacB(:,k) * PARCOMP(icomp).FacC(1,k);
                      A = NaN(numel( listHBOch),1);
                      idok = find(~isnan(listHBOch));
                      A(idok,1)=topo(listHBOch(idok));
                     alllabel = [alllabel,{PARCOMP(icomp).label}];                  
                     allbHbO = [allbHbO,A];
                     save(fullfile(pathoutlist,['TopoHbO',PARCOMP(icomp).label,'.mat']),'A' ,'zonelist','srsfile' );
                     clear A
                     %HBR parafac
                     topo=  PARCOMP(icomp).FacB(:,k) * PARCOMP(icomp).FacC(2,k);
                     A = NaN(numel( listHBOch),1);
                      idok = find(~isnan(listHBOch));
                      A(idok,1)=topo(listHBOch(idok));
                     alllabel = [alllabel,{PARCOMP(icomp).label}];                  
                     allbHbR = [allbHbR,A];    
                     save(fullfile(pathoutlist,['TopoHbR',PARCOMP(icomp).label,'.mat']),'A','zonelist','srsfile' );
                     clear A
                    alllabel = [alllabel,{PARCOMP(icomp).label}];
                    ncomp = ncomp +1
              elseif strcmp(upper(PARCOMP(icomp).type),'GLM')
                     listgood=PARCOMP(icomp).listgood;
                               
                     load(fullfile(pathoutlist,ZoneDtp{filenb}),'-mat');   
                       AHbO = nan(numel(zonelabel),1)
                       AHbR = nan(numel(zonelabel),1)
                      topo=PARCOMP(icomp).topo;
                     for izone = 1:numel(zonelabel)
                            [listHBOch, listHBRch]= findzoneinlistgood(NIRS, zone,listgood,zonelabel{izone});
                             AHbO(izone)= nanmean(topo(listHBOch));
                             AHbR(izone)= nanmean(topo(listHBRch));
                             
                     end
                     
                    
                     labelall = [labelall,{PARCOMP(icomp).label}];               
                     allbHbO = [allbHbO,AHbO];                    
                     allbHbR = [allbHbR,AHbR];
                   
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
                      save(fullfile(pathoutlist,['TopoHbO',PARCOMP(icomp).label,'.mat']),'A' ,'zonelist','srsfile' );
                      clear A
                   
                      A= nan(size(listHBRch,1),1);
                      idok = find(~isnan(listHBRch));
                      A(idok,1)=topo(listHBRch(idok));               
                      allbHbR = [allbHbR,A];
                      save(fullfile(pathoutlist,['TopoHbR',PARCOMP(icomp).label,'.mat']),'A','zonelist','srsfile' );
                      clear A            
                      ncomp = ncomp +1;
              end             
            end            
        end
    end
end 
1
  
    rawxls = [labelall',  [zonelabel';   num2cell( allbHbO)']]; 
    try
    xlswrite(fullfile(pathoutlist, ['zone',namelist, ext]),rawxls)
    catch
      [file,pathout] =uiputfile( fullfile(pathoutlist, ['zone',namelist, ext]))
       xlswrite(fullfile(pathout,file),rawxls)
    end
    out.NIRSmat = {fullfile(NIRSDtp{1},'NIRS.mat')};
