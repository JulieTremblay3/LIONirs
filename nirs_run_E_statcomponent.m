function out = nirs_run_E_statcomponent(job)
alphatr = job.e_statcomponent_alpha;
 dir1 = job.e_STATCOMPPath{1}; 
 if ~isdir(dir1)
     mkdir(dir1);
 end
try
   ttest(ones(20,1),zeros(20,1));   
catch
   disp('Uncomplete Stats components, Please Install Matlab Statistics and Machine Learning Toolbox')
   return
end
 


 if isfield(job.c_statcomponent,'b_TtestOneSample')
    AllC = [];
    for i=1:numel(job.c_statcomponent.b_TtestOneSample.f_component)
        load(job.c_statcomponent.b_TtestOneSample.f_component{i},'-mat');
        if size(A,1)> size(A,2)
            AllC = [AllC,A];
        else
             AllC = [AllC,A'];
        end
        try
        info.zonelistG1{i} = zonelist;    
        catch
            zonelist = [];
             info.zonelistG1{i} = [];
        end
    end  
    save(fullfile(dir1,'MatricebyChannel.mat'), 'AllC','-mat')
    disp(fullfile(dir1,'MatricebyChannel.mat'))
    info.AllG1 = AllC;
    info.G1 = job.c_statcomponent.b_TtestOneSample.f_component;
     for ich=1:size(AllC,1)       
        if 1
            if job.c_statcomponent.b_TtestOneSample.m_TtestOneSample == 1
                 info.tail = 'both';
            elseif job.c_statcomponent.b_TtestOneSample.m_TtestOneSample == 2                
                info.tail = 'left';
            elseif job.c_statcomponent.b_TtestOneSample.m_TtestOneSample == 3 ;
                info.tail = 'right';
            end      
            [h,p,ci,stats] = ttest(AllC(ich,:),0,'tail',info.tail); %Use from Matlab Statistics and Machine Learning Toolbox�  
            df(ich) = stats.df;
            tval(ich,1) = stats.tstat;
            pval(ich,1)= p;        
            mval(ich,1) = nanmean(AllC(ich,:));
            
        end
     end
     N_events=sum(~isnan(AllC),2)
     NBSUJET = numel(job.c_statcomponent.b_TtestOneSample.f_component);
     A = mval;
     save(fullfile(dir1,['ONESAMPLE_Mean n=',num2str(NBSUJET),'.mat']),'A','zonelist')  
     disp(['Save: ',fullfile(dir1,['ONESAMPLE_Mean n=',num2str(NBSUJET),'.mat'])]);
     A = tval;
     save(fullfile(dir1,['ONESAMPLE_Tmap n=',num2str(NBSUJET),'.mat']),'A','zonelist')
     disp(['Save: ',fullfile(dir1,['ONESAMPLE_Tmap n=',num2str(NBSUJET),'.mat'])]);
     A = mval.*double(pval<alphatr);
     save(fullfile(dir1,['ONESAMPLE_mean_',  num2str(alphatr),'unc.mat']), 'A','zonelist')  
     disp(['Save: ',fullfile(dir1,['ONESAMPLE_mean_',  num2str(alphatr),'unc.mat'])]);
     A = pval
      save(fullfile(dir1,['ONESAMPLE_pval_',  num2str(alphatr),'unc.mat']),'A','zonelist')
     disp(['Save: ',fullfile(dir1,['ONESAMPLE_pval_',  num2str(alphatr),'unc.mat'])]);
     [FDR,Q] = mafdr(pval); %establish fdr among all channels results         
     A = mval.*double(Q<alphatr);
     save(fullfile(dir1,['ONESAMPLE_mean_',  num2str(alphatr),'fdr.mat']),'A','zonelist')
     disp(['Save: ',fullfile(dir1,['ONESAMPLE_mean_',  num2str(alphatr),'fdr.mat'])]);
     A = Q;
     save(fullfile(dir1,['ONESAMPLE_Q_',  num2str(alphatr),'fdr.mat']),'A','zonelist')
     disp(['Save: ',fullfile(dir1,['ONESAMPLE_mean_',  num2str(alphatr),'fdr.mat'])]);
     A = N_events;
     save(fullfile(dir1,['ONESAMPLE_N_events_',  num2str(alphatr),'.mat']),'A','zonelist')
     disp(['Save: ',fullfile(dir1,['ONESAMPLE_N_events_',  num2str(alphatr),'.mat'])]);
elseif isfield(job.c_statcomponent,'b_Ttestpaired')   
    if job.c_statcomponent.b_Ttestpaired.m_TtestOneSample == 1
        info.tail = 'both';
    elseif job.c_statcomponent.b_Ttestpaired.m_TtestOneSample == 2  
        info.tail = 'left';
    elseif job.c_statcomponent.b_Ttestpaired.m_TtestOneSample == 3
        info.tail = 'right';
    end
    AllG1 = [];
    for i = 1:numel(job.c_statcomponent.b_Ttestpaired.f_componentG1)
        load(job.c_statcomponent.b_Ttestpaired.f_componentG1{i},'-mat');
        AllG1 = [AllG1, A];   
        info.zonelistG1{i} = zonelist;
    end 
 
    AllG2 = [];
    for i = 1:numel(job.c_statcomponent.b_Ttestpaired.f_componentG2)
        load(job.c_statcomponent.b_Ttestpaired.f_componentG2{i},'-mat');    
        AllG2 = [AllG2, A];
        info.zonelistG2{i} = zonelist;
    end
        info.G1 = job.c_statcomponent.b_Ttestpaired.f_componentG1;
        info.G2 = job.c_statcomponent.b_Ttestpaired.f_componentG2;
        info.AllG1 = AllG1;
        info.AllG2 = AllG2;
    for ich = 1:size(AllG1,1) 
       [h,p,ci,stats] = ttest(AllG1(ich,:),AllG2(ich,:), 'tail',info.tail); % from Matlab Statistics and Machine Learning Toolbox�
        tval(ich,1)   = stats.tstat;
        df            = stats.df;
        pval(ich,1)   = p;      
        mval(ich,1)   = nanmean(AllG1(ich,:)) - nanmean(AllG2(ich,:));
    end
    [FDR,Q] = mafdr(pval); %establish fdr among all channels results   
    
    NBSUJET = numel(job.c_statcomponent.b_Ttestpaired.f_componentG1); 
    A = nanmean(AllG1,2);
    save(fullfile(dir1,['PAIRED_Mean G1 n=',num2str(NBSUJET),'.mat']),'A','zonelist');    
    disp(['Save: ',fullfile(dir1,['PAIRED_Mean G1 n=',num2str(NBSUJET),'.mat'])]);
    A = nanmean(AllG2,2);
    save(fullfile(dir1,['PAIRED_Mean G2 n=',num2str(NBSUJET),'.mat']),'A','zonelist') ;    
    disp(['Save: ',fullfile(dir1,['PAIRED_Mean G2 n=',num2str(NBSUJET),'.mat'])]);
    A = mval;
    save(fullfile(dir1,['PAIRED_Mean G1-G2','.mat']),'A','zonelist');
    disp(['Save: ',fullfile(dir1,['PAIRED_Mean G1-G2','.mat'])]);
    A = tval;
    save(fullfile(dir1,['PAIRED_Tmap.mat']),'A','zonelist');    
    disp(['Save: ',fullfile(dir1,['PAIRED_Tmap.mat'])]);
    A = tval.*double(pval<alphatr);
    save(fullfile(dir1,['PAIRED_Tmap_',num2str(alphatr),'unc.mat']), 'A','zonelist');    
    disp(['Save: ',fullfile(dir1,['PAIRED_Tmap_',num2str(alphatr),'unc.mat'])]);
    A = mval.*double(pval<alphatr);
    save(fullfile(dir1,['PAIRED_mean G1-G2_',num2str(alphatr),'unc.mat']),'A','zonelist');     
    disp(['Save: ',fullfile(dir1,['PAIRED_mean G1-G2_',num2str(alphatr),'unc.mat'])]);
    A = mval.*double(Q<alphatr);
    save(fullfile(dir1,['PAIRED_mean G1-G2_',  num2str(alphatr),'fdr.mat']),'A','zonelist');
    disp(['Save: ',fullfile(dir1,['PAIRED_mean G1-G2_',  num2str(alphatr),'fdr.mat'])]);
elseif isfield(job.c_statcomponent,'b_TtestUnpaired')    
    if job.c_statcomponent.b_TtestUnpaired.m_TtestOneSample == 1
        info.tail = 'both';
    elseif job.c_statcomponent.b_TtestUnpaired.m_TtestOneSample == 2  
        info.tail = 'left';
    elseif job.c_statcomponent.b_TtestUnpaired.m_TtestOneSample == 3
        info.tail ='right'; 
    end
    AllG1 = [];
    for i = 1:numel(job.c_statcomponent.b_TtestUnpaired.f_componentG1)
        load(job.c_statcomponent.b_TtestUnpaired.f_componentG1{i},'-mat');        
         if size(A,1)>=size(A,2)
             AllG1 = [AllG1,A];
        elseif size(A,1)<=size(A,2)
             AllG1 = [AllG1,A'];
        end
        info.zonelistG1{i} = zonelist;
    end 
    
     AllG2 = [];
    for i = 1:numel(job.c_statcomponent.b_TtestUnpaired.f_componentG2)
        load(job.c_statcomponent.b_TtestUnpaired.f_componentG2{i},'-mat');      
        if size(A,1)>=size(A,2)
             AllG2 = [AllG2,A];
        elseif size(A,1)<=size(A,2)
             AllG2 = [AllG2,A'];
        end
        info.zonelistG2{i}=zonelist;
    end
        try
        if isfield(job.c_statcomponent.b_TtestUnpaired.c_statpermutation,'b_permutation')
            %not all channel could be use in the final matrice get the position
            %from one subjet and ensure you associate label in the Data.ZoneList to
            %a position in the zone.ml first col source second col detector
            load(job.c_statcomponent.b_TtestUnpaired.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.b_MCT_ZoneBased{1},'-mat')
            chanpos = zeros(numel(zonelist),3);
            for id=1:size(zonelist)
                name = zonelist{id};
                [DetL,SrsL] =strtok(name,' ');
                SDdetL = StrBoxy2SDDet(DetL);
                SDsrsL =str2num(SrsL(3:end));
                idml = find(zone.ml(:,1)==  SDsrsL & zone.ml(:,2)==  SDdetL & zone.ml(:,4)== 1);
                if ~isempty(idml)
                chanpos(id,:) =  zone.pos(idml,1:3);
                end
            end
        end
        catch
             disp('Multiple comparaison cluster could not find channel position, use all cluster')
        end
    
    
    
     info.AllG1=AllG1;
     info.G1=job.c_statcomponent.b_TtestUnpaired.f_componentG1;
     info.AllG2=AllG2; 
     info.G2=job.c_statcomponent.b_TtestUnpaired.f_componentG2;
 
     for ich = 1:size(AllG1,1)
        [h,p,ci,stats] = ttest2(AllG1(ich,:),AllG2(ich,:),'tail', info.tail ); 
        tval(ich,1)   = stats.tstat;
        df          = stats.df;
        pval(ich,1)   = p;      
        mval(ich,1)   = nanmean(AllG1(ich,:)) - nanmean(AllG2(ich,:));
     end 
     NBSUJET = numel(job.c_statcomponent.b_TtestUnpaired.f_componentG2)+numel(job.c_statcomponent.b_TtestUnpaired.f_componentG1);
     A = nanmean(AllG1,2);
     save(fullfile(dir1,['TWOSAMPLE_Mean G1 n=',num2str(numel(job.c_statcomponent.b_TtestUnpaired.f_componentG1)),'.mat']),'A','zonelist');
     disp(['Save: ',fullfile(dir1,['TWOSAMPLE_Mean G1 n=',num2str(numel(job.c_statcomponent.b_TtestUnpaired.f_componentG1)),'.mat'])]);
     A = nanmean(AllG2,2);
     save(fullfile(dir1,['TWOSAMPLE_Mean G2 n=',num2str(numel(job.c_statcomponent.b_TtestUnpaired.f_componentG2)),'.mat']),'A','zonelist');
     disp(['Save: ',fullfile(dir1,['TWOSAMPLE_Mean G2 n=',num2str(numel(job.c_statcomponent.b_TtestUnpaired.f_componentG2)),'.mat'])]);
     A = mval;
     save(fullfile(dir1,['TWOSAMPLE_Mean G1-G2_','.mat']),'A','zonelist');
     disp(['Save: ',fullfile(dir1,['TWOSAMPLE_Mean G1-G2_','.mat'])]);         
     A = tval;
     save(fullfile(dir1,['TWOSAMPLE_Tmap','.mat']),'A','zonelist');
     disp(['Save: ',fullfile(dir1,['TWOSAMPLE_Tmap','.mat'])]);
     A = tval.*double(pval<alphatr);
     save(fullfile(dir1,['TWOSAMPLE_Tmap_', num2str(alphatr),'unc.mat']), 'A','zonelist');
     disp(['Save: ',fullfile(dir1,['TWOSAMPLE_Tmap_', num2str(alphatr),'unc.mat'])]);
     A = mval.*double(pval<alphatr);
     save(fullfile(dir1,['TWOSAMPLE_mean_', num2str(alphatr),'unc.mat']),'A','zonelist');    
     disp(['Save: ',fullfile(dir1,['TWOSAMPLE_mean_', num2str(alphatr),'unc.mat'])]); 
     try
     [FDR,Q] = mafdr(pval);         
     A = mval.*double(Q<alphatr);
     save(fullfile(dir1,['TWOSAMPLE_mean_',num2str(alphatr),'fdr.mat']),'A','zonelist');
     disp(['Save: ',fullfile(dir1,['TWOSAMPLE_mean_',num2str(alphatr),'fdr.mat'])]);
     catch
         
     end
       if  isfield(job.c_statcomponent.b_TtestUnpaired.c_statpermutation,'b_permutation') 
           halfMAT = [AllG1,AllG2]
          Randgroupe = [ones(size( AllG1,2),1);ones(size( AllG2,2),1)*2] ;          
          iduse = [find(Randgroupe==1);find(Randgroupe==2)];
          nperm = str2num((job.c_statcomponent.b_TtestUnpaired.c_statpermutation.b_permutation.e_npermutation))
         for iperm=1:nperm
            permfix(iperm,:) =  randperm(numel(Randgroupe));
           
         end
         for idperm=1:nperm
         for ich = 1:size(AllG1,1)
                GR = permfix(idperm,:);
                [h,p,ci,stats] = ttest2(halfMAT(ich,GR(1:numel(find(Randgroupe==1)))), halfMAT(ich,GR((numel(find(Randgroupe==1))+1):end)));
                tvalperm(ich,idperm) =stats.tstat;
         end
         end
             t_dist_all =  tvalperm ;
        [fecdf,xecdf] = ecdf(t_dist_all(:));
        clustercritical=xecdf(sum(fecdf<(1-alphatr)));
        disp(['critical for clustering p=', num2str(alphatr),' T=',num2str(clustercritical)])
         neighbourdist = job.c_statcomponent.b_TtestUnpaired.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_neighbourdist;
          [stat, matneig] = FindClusterBasedPermutationInComponent(chanpos, neighbourdist,clustercritical, tval, tvalperm );
           
          for iposcluster  = 1:numel(stat.posclusters)
            id = find(stat.posclusterslabelmat==iposcluster);
            plotall  = zone.plot{1} ;
            plotLstall = zone.plotLst{1};
            zone.label{end+1} = ['Pos Cluster', num2str(iposcluster),' ',num2str(stat.posclusters(iposcluster).prob)];
            zone.plot{end+1} = plotall(id);
            zone.plotLst{end+1} =    plotLstall( id);
          end
            for inegcluster  = 1:numel(stat.negclusters)
                id = find(stat.negclusterslabelmat==inegcluster);
                plotall  = zone.plot{1} ;
                plotLstall = zone.plotLst{1};
                zone.label{end+1} = ['Neg Cluster ', num2str(inegcluster),' ',num2str(stat.negclusters(inegcluster).prob)];
                zone.plot{end+1} = plotall(id);
                zone.plotLst{end+1} =    plotLstall( id);
              end
          
           [filepathzone,namezone,ext] =fileparts(job.c_statcomponent.b_TtestUnpaired.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.b_MCT_ZoneBased{1});

           save(fullfile(dir1, ['Cluster',namezone,ext]),'zone','-mat');
           disp(['Save ', fullfile(dir1, ['Cluster',namezone,ext]),...
               ' with cluster']);
       end  

elseif isfield(job.c_statcomponent,'c_ANOVAN')
    if isfield(job.c_statcomponent.c_ANOVAN,'b_ANOVAN') %by channel     
     [filepath,name,ext] = fileparts(job.c_statcomponent.c_ANOVAN.b_ANOVAN.f_anovan{1});
     if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
         try
            [num,txt,raw] = xlsread(job.c_statcomponent.c_ANOVAN.b_ANOVAN.f_anovan{1});
         catch
            [num,txt,raw] = readtxtfile_asxlsread(job.c_statcomponent.c_ANOVAN.b_ANOVAN.f_anovan{1});
         end

     elseif  strcmp(ext,'.txt')
          [num,txt,raw] = readtxtfile_asxlsread(job.c_statcomponent.c_ANOVAN.b_ANOVAN.f_anovan{1});
     end
     
     %Case anovan each ch save each average and pval add fdr correction! 
     %load the observation if 
    [pathstr, name, ext]=fileparts(raw{2,3});

    AllCOM = []; 
    for i=2:size(raw,1)
        tmp         = load(fullfile(raw{i,1},raw{i,2}),'-mat'); % Load data  
        [filepath,name,ext] = fileparts(raw{i,3});
        if  isempty(filepath) %try to open in same folder
            fid = fopen(fullfile(raw{i,1},raw{i,3}));        %Load channel list
        else %try to open channel file directly at the good location folder
            fid = fopen(raw{i,3});        %Load channel list 
        end
        chlist      = textscan(fid, '%s%s');
        fclose(fid); 
        %find zone list index
        srs         = chlist{1};
        det         = chlist{2};
        idactual    = [];
        %verify the order is constant with the list export 
        for ichlist = 1:numel(chlist{1})
            id = strmatch([srs{ichlist},' ', det{ichlist}],tmp.zonelist,'exact');
            if ~isempty(id)
                idactual = [idactual,id];
            end
        end
        AllCOM = [AllCOM,tmp.A(idactual)];
        info.zonelist{i-1} = tmp.zonelist;
     end
     
     zonelist = tmp.zonelist(idactual);
     
     A = AllCOM'
     save(fullfile(dir1,'databychannel.mat'), 'A', 'zonelist')
     disp(['Save ', fullfile(dir1,'databychannel.mat')])
    groupcell = raw(2:end, 4:end);
    name = raw(1, 4:end);%conditioni  could be a number or a zone label to identify region if zone use 0 to not be consider ! 
    for i = 1:size(groupcell ,2)
        groupedef{i} = ( groupcell(:,i));
    end
    
    if ~isdir(fullfile(dir1,'ch'))
        mkdir(fullfile(dir1,'ch'))
    end
    info.groupedef =  groupedef{i};
    info.AllCOMP = AllCOM;
    for ich = 1:size(AllCOM,1)
        y = AllCOM(ich,:);          
        [p,tbl,stats,terms] = anovan(y',groupedef,'model','interaction','varnames',name,'display','off' );
        for idim = 1:numel(p)
            pval(ich,idim)      = p(idim);
            Fval(ich,idim)      = tbl{idim+1,6};
%             mult  = multcompare(stats,'display','off');
            [mult,m,h,gnames]   = multcompare(stats,'dimension',[1:size(terms,2)],'display','off');
            GRFac1              = [];

            for igname = 1:size(gnames,1)
                label           = gnames{igname};
                equaltoremove   = strfind(label,'=');
                commatoremove   = strfind(label,',');
                label([equaltoremove, commatoremove]) = [];
                gnames{igname}  = label;
            end
            for i=1:size(terms,2)
                GRFac1 = [GRFac1;stats.grpnames{i}];
            end
           % eval(stats.grpnames{1})
            for icomp = 1:size(mult,1)
%                eval( ['p',GRFac1{mult(icomp ,1)},'_',GRFac1{mult(icomp,2)},'(',num2str(ich),',',num2str(idim),')','=',num2str(mult(icomp,6)),';'])
%                eval( ['mdiff',GRFac1{mult(icomp ,1)},'_',GRFac1{mult(icomp,2)},'(',num2str(ich),',',num2str(idim),')','=',num2str(mult(icomp,4)),';'])
%               c1 = gnames{mult(icomp,1)}
%               c2 = gnames{mult(icomp,2)}
                eval( ['p',gnames{mult(icomp,1)},'_',gnames{mult(icomp,2)},'(',num2str(ich),',',num2str(idim),')','=',num2str(mult(icomp,6)),';']);
               eval( ['mdiff',gnames{mult(icomp,1)},'_',gnames{mult(icomp,2)},'(',num2str(ich),',',num2str(idim),')','=',num2str(mult(icomp,4)),';']);

            end
         end 
       save(fullfile([dir1,filesep,'ch'],[ zonelist{ich},'stats.mat']), 'stats','terms','p','tbl')
        % results = multcompare(stats,'Dimension',[1 2]);
    end 
    [FDR,Q] = mafdr(pval(:));
    FDR = reshape(FDR, size(pval));

    for idim = 1:size(pval,2)
        A = double( pval(:,idim)<0.05);
        factor = name(find(terms(idim,:)));
        label = [];
        for ifactor = 1:numel(factor)
            label = [label,factor{ifactor}];
        end
        save(fullfile( dir1,['anovan_p0.05 ',label,'.mat']),'A','zonelist'); 
        disp(['Save: ',fullfile( dir1,['anovan_p0.05 ',label,'.mat'])])
         A = double(1- pval(:,idim));
        save(fullfile( dir1,['anovan_1-p',label,'.mat']),'A','zonelist'); 
        disp(['Save: ',fullfile( dir1,['anovan_1-p',label,'.mat'])])
         A = double( Fval(:,idim));
        save(fullfile( dir1,['anovan_F ',label,'.mat']),'A','zonelist'); 
        disp(['Save: ',fullfile( dir1,['anovan_F ',label,'.mat'])])
        A = double(FDR(:,idim)<alphatr);
        save(fullfile( dir1,['anovan_FDRq',num2str(alphatr),label,'.mat']),'A','zonelist') ;
        disp(['Save: ',fullfile( dir1,['anovan_FDRq',num2str(alphatr),label,'.mat'])])
        A = double(1- FDR(:,idim));
        save(fullfile( dir1,['anovan_FDR1-q ',label,'.mat']),'A','zonelist'); 
        disp(['Save: ',fullfile( dir1,['anovan_FDR1-q ',label,'.mat'])])
        A = double(1- FDR(:,idim));
        for icomp = 1:size(mult,1)
            eval(['A = p',gnames{mult(icomp,1)},'_',gnames{mult(icomp,2)},';']);           
            fileout =fullfile(dir1,['posthoc p ',gnames{mult(icomp,1)},'_',gnames{mult(icomp,2)},'.mat']);
            save(fileout,'A','zonelist'); 
             eval(['A = mdiff',gnames{mult(icomp,1)},'_',gnames{mult(icomp,2)},';']);               
            fileout =fullfile(dir1,['posthoc mdiff ',gnames{mult(icomp,1)},'_',gnames{mult(icomp,2)},'.mat']);
            save(fileout,'A','zonelist');
            eval(['ptmp = p',gnames{mult(icomp,1)},'_',gnames{mult(icomp,2)},';']);  
            A = A.*(ptmp<alphatr);
            fileout =fullfile(dir1,['posthoc mdiffp',num2str(alphatr),gnames{mult(icomp,1)},'_',gnames{mult(icomp,2)},'.mat']);
            save(fileout,'A','zonelist');  
        end
        disp(['Save: ', fullfile(dir1,'posthoc')]);
    end 
    elseif  isfield(job.c_statcomponent.c_ANOVAN,'b_ANOVANzone')
      
     [filepath,name,ext] = fileparts(job.c_statcomponent.c_ANOVAN.b_ANOVANzone.f_anovan{1});
     if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
         try
         [num,txt,raw] = xlsread(job.c_statcomponent.c_ANOVAN.b_ANOVANzone.f_anovan{1});
         catch
            [num,txt,raw] = readtxtfile_asxlsread(job.c_statcomponent.c_ANOVAN.b_ANOVANzone.f_anovan{1});
         end
     elseif  strcmp(ext,'.txt')
         [num,txt,raw] = readtxtfile_asxlsread(job.c_statcomponent.c_ANOVAN.b_ANOVANzone.f_anovan{1});
     end
    
     %Case anovan each ch save each average and pval add fdr correction! 
     %load the observation if 
    [pathstr, name, ext]=fileparts(raw{2,3});
     AllCOM = [];
   
    
     for i=2:size(raw,1)       
         try
            tmp= load(fullfile(raw{i,1},raw{i,2}),'-mat');   
         catch
             disp(['Fail to open ',fullfile(raw{i,1},raw{i,2})])
         end
      try
          zone = load(fullfile(raw{i,1},raw{i,3}),'-mat'); 
      catch
          try
          zone = load(raw{i,3},'-mat');
          catch              
          disp(['Zone expected, file ',fullfile(raw{i,1},raw{i,3}),' could not be load'])
          return
          end
      end
        labelzone = raw{i,4};
        for izone = 1:numel(zone.zone.label)
            if strcmp(zone.zone.label{izone}, labelzone)
                goodzone = izone;
                %zoneidnum{i-1,1} = izone;
                zoneid{i-1,1} = labelzone;
            end
        end
        plotLst = zone.zone.plotLst{goodzone};
        AllCOM = [AllCOM,nanmax(tmp.A(plotLst))];
     end
 
     info.AllCOMP = AllCOM; 
     info.ROI =   raw(2:end,4);
     info.groupedef = raw(2:end,2);
     info.zone =  raw(2:end,3);
     
     groupcell = raw(2:end, 5:end);
     name = raw(1, 5:end);%conditioni  could be a number or a zone label to identify region if zone ! 

    %  groupedef{1} =  zoneid;
     for i = 1:size(groupcell ,2)
        groupedef{i} = ( groupcell(:,i));
     end
     try
    for ich = 1:size(AllCOM,1)         
         y = AllCOM(ich,:);         
         [p,tbl,stats,terms] = anovan(y',groupedef,'model','interaction','varnames',name,'display','off' );
         for idim = 1:numel(p)
            pval(ich,idim) = p(idim);
            Fval(ich,idim) = tbl{idim+1,6};
         end
         save(fullfile(dir1,['anovan', 'REGION','stat.mat']), 'stats','terms','p','tbl');
         disp(['Save: ',fullfile(dir1,['anovan', 'REGION','stat.mat'])]);
%          results = multcompare(stats,'Dimension',[1 2]);

    end
     catch
         disp('ANOVAN error') 
     end
    lastcol = [{'Val'};num2cell(y')];
    tmp = [raw,lastcol];
    try
        xlswrite(fullfile(dir1,'ROI.xls'),tmp)
    catch
        writetxtfile_asxlswrite(fullfile(dir1,'ROI.xls'),tmp);
    end
        disp(['Save extract amplitude by ROI to open outside matlab : winopen ', fullfile(dir1,'ROI.xls')])
    end
 elseif isfield(job.c_statcomponent,'c_StatcomponentExport')
     if isfield(job.c_statcomponent.c_StatcomponentExport,'f_componentexportch') %by channel     
     [filepath,name,ext] = fileparts(job.c_statcomponent.c_StatcomponentExport.f_componentexportch{1});
     if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
         try
            [num,txt,raw] = xlsread(job.c_statcomponent.c_StatcomponentExport.f_componentexportch{1});
         catch
            [num,txt,raw] = readtxtfile_asxlsread(job.c_statcomponent.c_StatcomponentExport.f_componentexportch{1});
         end

     elseif  strcmp(ext,'.txt')
          [num,txt,raw] = readtxtfile_asxlsread(job.c_statcomponent.c_StatcomponentExport.f_componentexportch{1});
     end
     
     %Case anovan each ch save each average and pval add fdr correction! 
     %load the observation if 
     try
    [pathstr, name, ext]=fileparts(raw{2,3});
    
    AllCOM = []; 
    for i=2:size(raw,1)
        tmp         = load(fullfile(raw{i,1},raw{i,2}),'-mat'); % Load data  
        if ~isempty(tmp.zonelist)
        [filepath,name,ext] = fileparts(raw{i,3});
        if  isempty(filepath) %try to open in same folder
            fid = fopen(fullfile(raw{i,1},raw{i,3}));        %Load channel list
        else %try to open channel file directly at the good location folder
            fid = fopen(raw{i,3});        %Load channel list 
        end
        chlist      = textscan(fid, '%s%s');
        fclose(fid); 
        %find zone list index
        srs         = chlist{1};
        det         = chlist{2};
        idactual    = [];
        %verify the order is constant with the list export 
        for ichlist = 1:numel(chlist{1})
            id = strmatch([srs{ichlist},' ', det{ichlist}],tmp.zonelist,'exact');
            if ~isempty(id)
                idactual = [idactual,id];
            end
        end
        AllCOM = [AllCOM,tmp.A(idactual)];
        info.zonelist{i-1} = tmp.zonelist;
        zonelist = tmp.zonelist(idactual);
        else
             AllCOM = [AllCOM,tmp.A];
             %disp('Warning be carefull to verify that the montage and channel are the same among all the experiment')
            zonelist = [];
        end 
    end
     catch
     disp('No channel list available data will only be add together in the default order')
       AllCOM = []; 
       percentagenoise = [];
        for i=2:size(raw,1)
            try
            tmp         = load(fullfile(raw{i,1},raw{i,2}),'-mat'); % Load data  
            AllCOM = [AllCOM,tmp.A]; 
            if isfield( tmp,'pourcentagenoise')
                percentagenoise = [percentagenoise ,tmp.pourcentagenoise];
            end
           
            zonelist = []; 
            catch
                tmp.A(:) = nan; 
                AllCOM = [AllCOM,tmp.A];       
                zonelist = [];  
                disp(['Missing ', fullfile(raw{i,1},raw{i,2})])
                percentagenoise = [percentagenoise, nan];
            end
        end     
     end

     

     
     A = AllCOM';
     info = raw;
     save(fullfile(dir1,'databychannel.mat'), 'A', 'zonelist','info','percentagenoise')
     disp(['Save ', fullfile(dir1,'databychannel.mat')])
   elseif  isfield(job.c_statcomponent.c_StatcomponentExport,'f_componentexportzone')
      
     [filepath,name,ext] = fileparts(job.c_statcomponent.c_StatcomponentExport.f_componentexportzone{1});
     if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
         try
         [num,txt,raw] = xlsread(job.c_statcomponent.c_StatcomponentExport.f_componentexportzone{1});
         catch
            [num,txt,raw] = readtxtfile_asxlsread(job.c_statcomponent.c_StatcomponentExport.f_componentexportzone{1});
         end
     elseif  strcmp(ext,'.txt')
         [num,txt,raw] = readtxtfile_asxlsread(job.c_statcomponent.c_StatcomponentExport.f_componentexportzone{1});
     end
    
     %Case anovan each ch save each average and pval add fdr correction! 
     %load the observation if 
    [pathstr, name, ext]=fileparts(raw{2,3});
     AllCOM = [];
   
    
     for i=2:size(raw,1)       
         try
            tmp= load(fullfile(raw{i,1},raw{i,2}),'-mat');   
         catch
             disp(['Fail to open ',fullfile(raw{i,1},raw{i,2})])
         end
      try
          zone = load(fullfile(raw{i,1},raw{i,3}),'-mat'); 
      catch
          try
          zone = load(raw{i,3},'-mat');
          catch              
          disp(['Zone expected, file ',fullfile(raw{i,1},raw{i,3}),' could not be load'])
          return
          end
      end
        labelzone = raw{i,4};
        for izone = 1:numel(zone.zone.label)
            if strcmp(zone.zone.label{izone}, labelzone)
                goodzone = izone;
                %zoneidnum{i-1,1} = izone;
                zoneid{i-1,1} = labelzone;
            end
        end
        plotLst = zone.zone.plotLst{goodzone};
        AllCOM = [AllCOM,nanmax(tmp.A(plotLst))];
     end
 
     info.AllCOMP = AllCOM; 
     info.ROI =   raw(2:end,4);
     info.groupedef = raw(2:end,2);
     info.zone =  raw(2:end,3);
     
     groupcell = raw(2:end, 5:end);
     name = raw(1, 5:end);%conditioni  could be a number or a zone label to identify region if zone ! 

     A = AllCOM';
       lastcol = [{'Val'};num2cell(A)];
    info = [raw,lastcol];
     save(fullfile(dir1,'databyzone.mat'), 'A','info')
     disp(['Save ', fullfile(dir1,'databyzone.mat')])
    try
        xlswrite(fullfile(dir1,'ROI.xls'),info)
    catch
        writetxtfile_asxlswrite(fullfile(dir1,'ROI.xls'),info);
    end
        disp(['Save extract amplitude by ROI to open outside matlab : winopen ', fullfile(dir1,'ROI.xls')])
     end
 end
if job.m_statcomponent_saveoption==0
    save(fullfile(dir1,['info.mat']),'info')
end


matlabbatch{1}.spm.tools.nirsHSJ.M_dataComponent.E_statcomponent = job;
save(fullfile(dir1,'Batch_stat.mat'),'matlabbatch')
disp(['Save Batch_stat.mat'])
disp(['Results saved in: ', dir1])

out = 'no NIRS.mat dependancy';


function [fdr, q, pi0, rs] = mafdr(p, varargin)
%MAFDR estimates false discovery rates (FDR) of multiple-hypothesis testing.
%
%   FDR = MAFDR(P) estimates the positive false discovery rate (pFDR) from
%   the p-values P of multiple-hypothesis testing using the procedure
%   described by Storey (2002). P is a numeric vector or a DataMatrix
%   object.
%
%   [FDR,Q] = MAFDR(P) returns the q-values. The q-value is a hypothesis
%   testing error measure for each observation with respect to the pFDR. 
%
%   [FDR,Q,PIO] = MAFDR(P) returns the estimated a-priori probability that
%   the null hypothesis is true.
%
%   MAFDR(...,'LAMBDA',L) sets the tuning parameter lambda used in
%   estimating the a-priori probability that the null hypotheses is true. L 
%   can be a scalar or a vector with at least 4 elements. All the values in
%   L must be greater than 0 and less than 1. If a vector is specified,
%   MAFDR automatically chooses the tuning parameter lambda. L defaults to
%   [0.01:0.01:0.95].
%
%   MAFDR(...,'METHOD',METHOD) sets a method to automatically choose the
%   tuning parameter lambda. METHOD can be 'bootstrap' (default) or
%   'polynomial' (the cubic polynomial). 
%
%   [FDR,Q,PIO,RS] = MAFDR(...) returns the square of correlation
%   coefficient RS if using the 'polynomial' method to choose lambda.
%
%   FDR = MAFDR(P,'BHFDR',TF) returns the FDR adjusted p-values by
%   linear-step up (LSU) procedure originally introduced by Benjamini and
%   Hochberg (1995) if TF is set true. TF defaults to false.
%
%   MAFDR(...,'SHOWPLOT',TF) displays the estimated a-priori probability
%   that the null hypotheses is true as a function of the tuning parameter
%   lambda, and the plot of q-values vs. p-values if TF is set to true and
%   using Storey procedure. If TF is set to true, and 'BHFDR' is true,
%   MAFDR diplays the plot of FDR adjusted p-values vs. p-values. TF
%   defaults to false.
% 
%   Example:
%       load prostatecancerexpdata;
%       p = mattest(dependentData, independentData, 'permute', true);
%       [fdr, q] = mafdr(p, 'showplot', true);
%   
%   See also AFFYGCRMA, AFFYRMA, GCRMA, MAIRPLOT, MALOGLOG, MAPCAPLOT,
%   MATTEST, MAVOLCANOPLOT, RMASUMMARY. 

% Copyright 2006-2011 The MathWorks, Inc.

% References: 
% [1] J.D. Storey. "A direct approach to false discovery rates",
%     Journal of the Royal Statistical Society, B (2002), 64(3),
%     pp.479-498.
% [2] J.D. Storey, and R. Tibshirani. "Statistical significance for
%     genomewide studies", Proc.Nat.Acad.Sci. 2003, 100(16), pp. 9440-9445.
% [3] J.D. Storey, J.E. Taylor and D. Siegmund. "Strong control,
%     conservative point estimation, and simultaneous conservative
%     consistency of false discovery rates: A unified approach", Journal of
%     the Royal Statistical Society, B (2004), 66, pp. 187-205.
% [4] Y.Benjamini, Y.Hochberg. "Controlling the false discovery rate: a
%     practical and powerful approach to multiple testing. Journal of the
%     Royal Statistical Society, B (1995), 57, pp.289-300.

% Check the inputs
if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

bioinfochecknargin(nargin,1,mfilename);

if isa(p, 'bioma.data.DataMatrix')
    rownames = p.RowNames;
    p = p.(':')(':');
else
    rownames = [];
end

if ~isnumeric(p) || ~isreal(p) || ~isvector(p)
    error(message('bioinfo:mafdr:PValuesNotNumericAndRealVector'))
end

% Initialization
bhflag = false;
showplotflag = false;
lambda = (0.01:0.01:0.95);
bootflag = true;
lambdaGiven = false;
methodGiven = false;

% deal with the various inputs
if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:mafdr:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'bhfdr', 'lambda', 'method', 'showplot'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:mafdr:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:mafdr:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % BH FDR
                    bhflag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 2 % lambda
                    if ~isnumeric(pval) && isvector(pval) 
                        error(message('bioinfo:mafdr:LambdaMustBeNumericAndVector'));
                    end
                    if numel(pval)> 1 && numel(pval) <4
                        error(message('bioinfo:mafdr:badLambdaRange'))
                    end
                    lambda = sort(pval);
                    if any(lambda<=0) || any(lambda>=1) || any(diff(lambda)==0)
                        error(message('bioinfo:mafdr:badLambdaValues'))
                    end
                    lambdaGiven = true;
                case 3 %lambda method
                    bootflag = bioinfoprivate.optPartialMatch(pval,{'bootstrap', 'polynomial'},okargs{k},mfilename)==1;
                    methodGiven = true;
                case 4 % showplot
                    showplotflag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
            end
        end
    end
end

if bhflag && nargout>1
    error(message('bioinfo:mafdr:maxlhsbh'))
end
if ~bhflag && bootflag && nargout>3
    error(message('bioinfo:mafdr:maxlhspn'))
end
if  bhflag && (methodGiven || lambdaGiven)
    warning(message('bioinfo:mafdr:LambdaAndMethodIgnored'))
end
if methodGiven && isscalar(lambda)
    warning(message('bioinfo:mafdr:MethodIgnored'))
end

q=[];
pi0=[];
rs = [];
nanidx = isnan(p);
fdr = p;
if bhflag
    fdr(~nanidx) = bhfdr(p(~nanidx));
    
    if ~isempty(rownames)
       fdr = bioma.data.DataMatrix(fdr, rownames, {'FDR'});
    end
    
    if showplotflag
        if isempty(rownames)
            showBHplot(p,fdr)
        else
            showBHplot(p,fdr.(':')(':'))
        end
    end
else
    lambda = lambda(:);
    q = fdr;
    [fdr(~nanidx), q(~nanidx), pi0_all, pi0, cs, rs] = storeyFDR(p(~nanidx), lambda, bootflag);
    if ~isempty(rownames)
       fdr = bioma.data.DataMatrix(fdr, rownames, {'FDR'});
       q = bioma.data.DataMatrix(q, rownames, {'q-values'});
    end
    
    if showplotflag
        if isempty(rownames)
            showplots(lambda, p, q, pi0_all, pi0, cs)
        else
            showplots(lambda, p, q.(':')(':'), pi0_all, pi0, cs)
        end
    end
end

%--------------- Helper functions ----------------------%
function [fdr, q, pi0_all, pi0, cs, rs] = storeyFDR(p, lambda, bootflag)
% Compute pFDR and q-values using Storey and Tibshirani procedures
% estimate the proportion of genes that are truly null with a tuning
% parameter lambda.
m = numel(p);
pi0_all = estimatePI0(p, lambda);
cs = [];
rs = [];

% find the true pi0 from a range of lambda
if numel(lambda) > 1
    nan_idx = isnan(pi0_all);
    if bootflag
        pi0 = bootstrapchooser(pi0_all(~nan_idx), lambda(~nan_idx), p);
    else
        [pi0, cs, rs] = polychooser(pi0_all(~nan_idx), lambda(~nan_idx));
    end
else
    pi0 = pi0_all;
end

if pi0 <= 0
    error(message('bioinfo:mafdr:BadEstimatedPI0Value'));
end

% Estimate the positive FDR for each gene
[~, idx] = sort(p);
v = m * p;
r = zeros(size(v));
r(idx) = 1:m;

fdr = min(pi0 * v./r, 1);

% Calculate the q-value for each gene
qord = fdr(idx);

% Running min in reverse order
qord = cummin(qord, 'reverse');

q(idx) = qord;
q=q(:);

%------------------------------------------------------------
function [pi0,cs, rs] = polychooser(pi0_all, lambda)
% Choosing lambda for the estimate of pi0 using cubic polynomial of pi0 on
% lambda. Also returns the goodness of fit r^2, r-correlation coefficient.

[cs, s] = polyfit(lambda, pi0_all, 3);

% Compute pi0
d = diff(polyval(cs,lambda)) ./diff(lambda);
[mind, idx] = min(abs(d));
if mind < 0.005
    pi0 = polyval(cs, lambda(idx));
else
    pi0 = polyval(cs, max(lambda));
end

if pi0 > 1
     warning(message('bioinfo:mafdr:PoorEstimatedPI0Value'));
end
pi0 = min(pi0, 1);

% compute r^2 - rs
sd = std(pi0_all);
n = numel(lambda) - 1;
rs = 1 - s.normr^2/(n*sd^2);

if rs < 0.90 
    warning(message('bioinfo:mafdr:PoorCubicPolynomialFit',sprintf('%0.4f',rs)));
end


%--------------------------------------------------------
function pi0 = bootstrapchooser(pi0, lambda, p)
% Storey+Taylor+Siegmund(2004)
% Choosing lambda for the estimate of pi0 using bootstrap sampling on
% p-values. 
min_pi0 = min(pi0);

B = 100; % number of bootstrap replicates
m = numel(p);
n = numel(lambda);

% Get the indices of resampling, and extract these from the data
inds = unidrnd(m, m, B);
p_boot = p(inds);

mse = zeros(n,1);
mseCount = zeros(n,1);

% Bootstrap, and compute MSE
for i = 1:B
    pi0_boot = estimatePI0(p_boot(:,i), lambda);
    nnanidx = ~isnan(pi0_boot);
    mse(nnanidx) = mse(nnanidx) + (pi0_boot(nnanidx) - min_pi0).^2;
    mseCount(nnanidx) = mseCount(nnanidx)+1;
end

mse = mse ./ mseCount .* max(mseCount);
mse(mseCount<=(max(mseCount)./2)) = inf;

% Find the minimum MSE. 
[~,minmse_idx] = min(mse);
pi0 = pi0(minmse_idx);

if pi0 > 1
     warning(message('bioinfo:mafdr:PoorEstimatedPI0Value'));
     pi0 = min(pi0,1);
end

%----------------------------------------------------------
function pi0 = estimatePI0(p, lambda)
% Storey+Tibshirani (2002)
% p - p-values
% lambda - tuning parameter

[F, p0] = ecdf(p);
pi0 = interp1q(p0, 1-F, lambda) ./ (1-lambda);

%------------------------------------------------------------
function showplots(lambda, p, q, pi0_all, pi0, cs)
if numel(lambda) ==1
    max_lambda = max(0.95,lambda);
    min_lambda = min(0.01,lambda);
    lambda_x = (min_lambda:0.02:max_lambda)';
    pi0_y = estimatePI0(p, lambda_x);
    lambda_pi0 = lambda;
else
    lambda_x = lambda;
    pi0_y = pi0_all;
    max_lambda = max(lambda);
    min_lambda = min(lambda);
    lambda_pi0 = max(lambda);
end

% add the smoother line
if isempty(cs)
    idx = isnan(pi0_y);
    cs = polyfit(lambda_x(~idx),pi0_y(~idx),3);
end

lambda_xx = linspace(min_lambda, max_lambda, 101);
pi0_yy = polyval(cs, lambda_xx);

miny = min(pi0_all) - 0.1;

% % xx = [lambda_pi0, lambda_pi0;lambda_pi0,1];
% % yy = [miny+0.01, pi0; pi0, pi0];
xx = [0, 1];
yy = [pi0, pi0];

subplot(2,1,1)
l1 = plot(lambda_x, pi0_y, 'r.'); %#ok
hold on;
l2 = plot(lambda_xx, pi0_yy, 'b-', 'LineWidth', 2);
l3 = plot(lambda_pi0, pi0, 'ko', 'LineWidth',2);
ylim([miny, 1]);

plot(xx(:), yy(:), 'k:');
hx = xlabel('');
hy = ylabel('');
ht = title('');
set([hx hy ht], 'interpreter', 'latex', 'FontSize', 12)
set(hx,'string', '$\lambda$');
set(hy, 'string', '$\hat\pi_0(\lambda)$')
set(ht, 'string', ['$\hat\pi_0$=', sprintf('%.4f', pi0)]) 
hold off
set(gca,'Tag','PI0vslambdaAxes')

hl = legend([l2, l3], 'cubic polynomial fit', 'testing');
legend_str = {'cubic polynomial fit', '$\hat\pi_0$'};
set(hl, 'interpreter', 'latex','string', legend_str);

subplot(2,1,2)
[p,idx] = sort(p);
q = q(idx);
plot(p,q);
xlabel('p-value','FontSize', 12 )
ylabel('q-value', 'FontSize', 12)
set(gca,'Tag','pvsqAxes')

%------------
function fdr = bhfdr(p)
% Compute the fdr using BH Linear step-up procedure

m = numel(p);
[p_ord, idx] = sort(p);
fdr_ord =  p_ord(:) .* (m./(1:m))';
% Running min in reverse order (in-place)
fdr_ord= cummin(fdr_ord, 'reverse');
fdr(idx) = fdr_ord;

%------------------------------------------------------------
function showBHplot(p,fdr)
[p, idx]= sort(p);
fdr = fdr(idx);
plot(p, fdr);
xlabel('p-value')
ylabel('FDR adjusted p-value')
function bioinfochecknargin(numArgs,low,name)
%BIOINFOCHECKNARGIN Validate number of input arguments
%
%   BIOINFOCHECKNARGIN(NUM,LOW,FUNCTIONNAME) throws an MException if the
%   number of input arguments NUM to function FUNCTIONNAME is less than the
%   minimum number of expected inputs LOW.
% 
%   Example
%      bioinfochecknargin(nargin, 3, mfilename)
%
%   See also MFILENAME, NARGCHK, NARGIN, NARGOUT, NARGOUTCHK.

%   Copyright 2007-2012  The MathWorks, Inc.

    
if numArgs < low
    msg = getString(message('bioinfo:bioinfochecknargin:NotEnoughInputs'));
    msgId = sprintf('bioinfo:%s:NotEnoughInputs',name);
    x = MException(msgId,msg);
    x.throwAsCaller;    
end
