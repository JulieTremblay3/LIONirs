function out = nirs_run_E_statmatrix(job)

xlslistfile = job.f_matrix{1};
[~,~,ext] =fileparts([xlslistfile]);
if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
    [raw, txt, info]=xlsread([xlslistfile]);
elseif strcmp(ext,'.txt')
    [num, txt, info] = readtxtfile_asxlsread([xlslistfile])
end
%Load the matrices same as any
groupeall = [];
for isubject=2:size(info,1)
    id = isubject-1;
    MAT = load(fullfile(info{isubject,1},[ info{isubject,2},'.mat']));
    DATA{id}.ZoneList = MAT.ZoneList;
    if isfield(MAT, 'meancorr')
        matcorr = MAT.meancorr;
        if job.m_fishertransform==1 % add fisher transform  
            matcorr =  1/2*log((1+matcorr)./(1-matcorr));
            DATA{id}.MAT = matcorr;
        else
            DATA{id}.MAT = matcorr;
        end
    end
    
    DATA{id}.name = info{isubject,2};
    DATA{id}.MATtrial =  MAT.matcorr;
    DATA{id}.GR = info{isubject,4};
    infocov = [];
    if size(info,2) > 4
        %Define Additional covariable
        idcov = 1 ;
        for icolumn = 5:size(info,2)
            eval( ['DATA{id}.cov',num2str(idcov),' = ', num2str(info{isubject,icolumn}),';']);
            idcov=idcov+1;
            if isubject==2
                infocov = [infocov ;info(1,icolumn)];
            end
        end
    end
    
    
    DATA{id}.System = 'ISS';
    load(fullfile(info{isubject,1}, info{isubject,3}),'-mat');
    names = fieldnames(zone);
    for iname = 1:numel(names)
        eval(['DATA{id}.zone.',names{iname},' =zone.',names{iname},';']);
    end
    list_subject{id} =DATA{id}.name;
    groupeall = [groupeall; info{isubject,4}];
    clear MAT
end

if job.m_nodeunit==1 %channel mode
    idsubject = 1:numel(groupeall);
    for isubject = 1:numel(groupeall)
        MATall(isubject,:,:)=DATA{isubject}.MAT;
        groupid(isubject)= DATA{idsubject(isubject)}.GR;
    end
    ZONEid = [info{end,3}];
    ZoneList =  DATA{end}.ZoneList;
    labelnode = 'c';
elseif  job.m_nodeunit==2
    MATall =zeros(numel(DATA),numel(DATA{id}.zone.label),numel(DATA{id}.zone.label));
    for isubject = 1:numel(groupeall)
        List = DATA{isubject}.ZoneList;
        for izone = 1:numel(DATA{isubject}.zone.label)
            ML = DATA{isubject}.zone.ml;
            DATA{isubject}.zone.plotLst;
            idlisti = [];
            idliststr = [];
            chzone = DATA{isubject}.zone.plotLst{izone};
            for ichzone = 1:numel(chzone);
                ich = chzone(ichzone);
                if strcmp(DATA{isubject}.System,'ISS')
                    strDet = SDDet2strboxy_ISS(ML(ich,2));
                    strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                    idch = strmatch([strDet, ' ',strSrs ],List,'exact');
                end
                idliststr =[idliststr,{[strDet, ' ',strSrs ]}];
                idlisti = [idlisti, idch];
            end
            for jzone = 1:numel(DATA{isubject}.zone.label)
                idlistj = [];
                chzone = DATA{isubject}.zone.plotLst{jzone};
                for ichzone = 1:numel(chzone);
                    ich = chzone(ichzone);
                    if strcmp(DATA{isubject}.System,'ISS')
                        strDet = SDDet2strboxy_ISS(ML(ich,2));
                        strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                        idch = strmatch([strDet, ' ',strSrs ],List,'exact');
                    end
                    idlistj = [idlistj, idch];
                end
                matROI = DATA{isubject}.MAT(idlisti,idlistj);
                id = find(matROI==0);
                if isempty(id)
                    matROI(id)=nan;
                end
                MATall(isubject,izone,jzone) = nanmean(matROI(:));
                if izone==jzone
                    matnbnanbyizone(isubject,izone)=numel(find(sum(double(isnan(matROI))) ==size(matROI,1)));
                    matnbtotchbyizone(isubject,izone) = size(matROI,1);
                end
            end
            
            
            
        end
        groupid(isubject)= DATA{isubject}.GR;
        labelnode = 'z';
    end
    zoneuse=DATA{isubject}.zone;
    ZoneList = [];
    plottmp=[];
    plotLst = [];
    for izoneList = 1:size(MATall,2)
        MLfake(izoneList,1) = izoneList;%source
        MLfake(izoneList,2) = 1; %detecteur
        MLfake(izoneList,3) = 1;
        MLfake(izoneList,4) = 1;
        strDet = SDDet2strboxy_ISS(MLfake(izoneList,2));
        strSrs = SDPairs2strboxy_ISS(MLfake(izoneList,1));
        ZoneLabel{izoneList,1}=zoneuse.label{izoneList};
        ZoneList{izoneList,1} = [strDet,' ', strSrs];
        plottmp{izoneList} = [izoneList,1];
        plotLst{izoneList} = [izoneList];
    end
    %save zone list associate
    zone.plot = plottmp;
    zone.plotLst = plotLst;
    zone.label = ZoneLabel;
    zone.color = zone.color;
    zone.ml = MLfake;
    zone.chMAT = plotLst;
    save(fullfile(info{isubject,1},['avg', info{isubject,3}]),'zone','-mat');
    ZONEid = ['avg', info{isubject,3}];
end

 [filepath,name,ext] = fileparts(xlslistfile);
if isfield(job.c_statmatrix,'b_TtestOneSamplematrix') 
    AllC = []; 
    id =1;
    %Use a specific groupe
    idG1 = find(  groupeall==job.c_statmatrix.b_TtestOneSamplematrix.e_TtestOneSampleGR);
    MATallG1 = MATall( idG1,:,:);
    meanall = squeeze(nanmean( MATallG1,1));
    refval  = 0.1;
    tval = squeeze((nanmean(MATallG1(:,:,:),1)-refval)./(nanstd(MATallG1(:,:,:),1)./sqrt(sum(~isnan(MATallG1(:,:,:)),1))  ));
    dfall = squeeze(sum(~isnan(MATallG1(:,:,:)),1))-1;
    for i=1:size(tval,1)
        for j=1:size(tval,2)
            try
                % Compute the correct p-value for the test, and confidence intervals
                % if requested.
                if job.c_statmatrix.b_TtestOneSamplematrix.m_TtestOneSample_matrix == 1 % two-tailed test
                     pval(i,j) = 2 * tcdf(-abs(-tval(i,j)), dfall(i,j));
                    %     if nargout > 2
                    %         crit = tinv((1 - alpha / 2), df) .* ser;
                    %         ci = cat(dim, xmean - crit, xmean + crit);
                    %     end
                elseif job.c_statmatrix.b_TtestOneSamplematrix.m_TtestOneSample_matrix == 3 % right one-tailed test
                     pval(i,j) = tcdf(-tval(i,j), dfall(i,j));
                    %     if nargout > 2
                    %         crit = tinv(1 - alpha, df) .* ser;
                    %         ci = cat(dim, xmean - crit, Inf(size(p)));
                    %     end
                elseif job.c_statmatrix.b_TtestOneSamplematrix.m_TtestOneSample_matrix == 2 % left one-tailed test
                     pval(i,j) = tcdf(tval(i,j), dfall(i,j));
                    %     if nargout > 2
                    %         crit = tinv(1 - alpha, df) .* ser;
                    %         ci = cat(dim, -Inf(size(p)), xmean + crit);
                    %     end
                end
            catch
                pval(i,j) = nan;
            end
        end
    end
    dir1 = job.e_statmatrixPath{1};
    
    if ~isdir(dir1)
        mkdir(dir1)
    end
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];

    
    matcorr =  meanall;
    meancorr = meanall;
    totaltrialgood = mean(dfall(:));    
    file = [name,labelnode,'OneSampleTtest mean','.mat'];
    save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr','totaltrialgood');
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];
    
    
    % ZoneList = MAT.ZoneList;
    matcorr =  tval;
    meancorr = tval;
    totaltrialgood = mean(dfall(:));    
    file = [name,labelnode,'OneSampleTtest tval','.mat'];
    save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr','totaltrialgood');
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];
    
    file = [name,labelnode,'OneSampleTtest tvalp05','.mat'];
    matcorr = tval.*double(pval<0.05);
    meancorr = tval.*double(pval<0.05);
    save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr','totaltrialgood');
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];
    [FDR,Q] = mafdr(pval(:));
    Q = reshape(Q,size(pval));
    infonew = [infonew;new];
    
    file = [name,labelnode,'OneSampleTtest meanp05','.mat'];
    matcorr = meanall.*double(pval<0.05);
    meancorr = meanall.*double(pval<0.05);
    save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr','totaltrialgood');
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];
    
    file = [name,labelnode,'OneSampleTtest tvalFDR05','.mat'];
    matcorr = tval.*double(Q<0.05);
    meancorr = tval.*double(Q<0.05);
    save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr','totaltrialgood');
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];
    
    
    [FDR,Q] = mafdr(pval(:));
    Q = reshape(Q,size(pval));
    file = [name,'OneSampleTtest tvalFDR01','.mat'];
    matcorr = tval.*double(Q<0.01);
    meancorr = tval.*double(Q<0.01);
    save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr','totaltrialgood');
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];
    copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid))
    %  dir1 = job.e_statmatrixPath{1};

    if ismac
        % Code to run on Mac platform problem with xlswrite
        writetxtfile(fullfile(dir1,[name,labelnode,'SimpleTtest.txt']),infonew);
        disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'SimpleTtest.txt'])]);
    else
        xlswrite(fullfile(dir1,[name,labelnode,'SimpleTtest.xlsx']),infonew);
        disp(['Result .xlsx file saved' fullfile(dir1,[name,labelnode,'SimpleTtest.xlsx'])])
    end
elseif isfield(job.c_statmatrix,'b_PermutationTest')
       ESTAD = 4;
        INDEP = 1;
        NPERM = str2num(job.c_statmatrix.b_PermutationTest.e_npermutation);
        g1 = find(groupid==job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR);
        g2 = find(groupid==job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR2);
         
        cb1 = MATall(g1,:,:);
        pb1 = MATall(g2,:,:);   
        %initialised 
        Toij=zeros(size( cb1,2));
        ncb1=zeros(size( cb1,2));
        npb1 =zeros(size( pb1,2));
   
        if NPERM<=1 %permuation =1 do unpaired test
            mean_G1=squeeze(nanmean(cb1,1));
            mean_G2=squeeze(nanmean(pb1,1));
            var_G1 = squeeze(nanvar(cb1,0,1));
            var_G2 = squeeze(nanvar(pb1,0,1));
            n_G1 = squeeze(sum(~isnan(cb1),1));
            n_G2 = squeeze(sum(~isnan(pb1),1));
            df = n_G1+ n_G2-2;                 
            varc = (1./n_G1+ 1./n_G2 ).*((n_G1-1).*var_G1 + (n_G2-1).* var_G2 )./df;
            Toij= (mean_G1-mean_G2)*1./sqrt(varc);
            ncb1 = n_G1;
            npb1 = n_G2;
            
        for i=1:size(Toij,1)
        for j=1:size(Toij,2)
            try
                % Compute the correct p-value for the test
                if 1 % two-tailed test pval
                     FUniv(i,j) = 2 * tcdf(-abs(- Toij(i,j)), df(i,j));
                end

            catch
                FUniv(i,j) = nan;
            end
        end
    end
            
  
    [FDRmat,pcritic] = mafdr(FUniv(:));
    %figure;cdfplot(FUniv(:))
        %WRITE IN A NEW FILE
        else
        if (sum(isnan(pb1(:)))+sum(isnan(cb1(:))))>0 %if presence NAN do the stat without them...
            for i=1:size(pb1,2) % row
                for j=1:size(cb1,3) %col
                    tmpcb1 = cb1(:,i,j);
                    tmppb1 = pb1(:,i,j);
                    idnan = find(isnan(tmpcb1));
                    if ~isempty(idnan)
                        tmpcb1(idnan)=[];
                    end
                    idnan = find(isnan(tmppb1));
                    if ~isempty(idnan)
                        tmppb1(idnan)=[];
                    end
                    if Toij(i,j)==0
                    [FSupSup,FSupDeriv,FSupTime,pij,tij] = TestPermut2Grupos(ESTAD,INDEP,tmpcb1,tmppb1,NPERM);
                    %apply symetric
                    ncb1(i,j) = numel(tmpcb1);
                    npb1(i,j) = numel(tmppb1);
                    FUniv(i,j) = pij;
                    Toij(i,j) = tij;
                     ncb1(j,i) = numel(tmpcb1);
                    npb1(j,i) = numel(tmppb1);
                    FUniv(j,i) = pij;
                    Toij(j,i) = tij;
                    
                    end
                end
            end
        else
            [FSupSup,FSupDeriv,FSupTime,FUniv,Toij] = TestPermut2Grupos(ESTAD,INDEP,pb1,cb1,NPERM);
            ncb1 = ones(size(cb1,2),size(cb1,3))*numel(size(cb1,1));  
            npb1 = ones(size(pb1,2),size(pb1,3))*numel(size(pb1,1));
        end
        end
        MeanG1 = squeeze(nanmean(cb1));
        MeanG2 = squeeze(nanmean(pb1));
        
        
        dir1 = job.e_statmatrixPath{1};
        if ~isdir(dir1)
            mkdir(dir1)
        end
        
        
        %WRITE IN A NEW FILE
        infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
        file = [name,'_',labelnode,num2str(NPERM),'permutation tstat','.mat'];
        matcorr = Toij;
        meancorr = Toij;
        save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
        new = [{dir1},{file}, {ZONEid},{1} ];
        infonew = [infonew;new];
        
        file = [name,'_',labelnode,num2str(NPERM),'permutation 1-pval','.mat'];
        matcorr = 1-FUniv;
        meancorr = 1-FUniv;
        save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
        new = [{dir1},{file}, {ZONEid},{1} ];
        infonew = [infonew;new];
              
        
        file = [name,'_',labelnode,num2str(NPERM),'permutation G1-G2','.mat'];
        matcorr = real(squeeze((nanmean(cb1,1)- nanmean(pb1,1))));
        meancorr = real(squeeze((nanmean(cb1,1)- nanmean(pb1,1))));
        save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
        new = [{dir1},{file}, {ZONEid},{1} ];
        infonew = [infonew;new];
        

        file = [name,'_',labelnode,num2str(NPERM),'permutation G2-G1','.mat'];
        matcorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1)));
        meancorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1)));
        save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
        new = [{dir1},{file}, {ZONEid},{1} ];
        infonew = [infonew;new];
        
        file = [name,'_',labelnode,num2str(NPERM),'permutation G1-G2 p05','.mat'];
        matcorr = real(squeeze(nanmean(cb1,1))-squeeze(nanmean(pb1,1))).*double(FUniv<0.05);
        meancorr = real(squeeze(nanmean(cb1,1))-squeeze(nanmean(pb1,1))).*double(FUniv<0.05);
        save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
        new = [{dir1},{file}, {ZONEid},{1} ];
        infonew = [infonew;new];
        
        file = [name,'_',labelnode,num2str(NPERM),'permutation G2-G1 p05','.mat'];
        matcorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1))).*double(FUniv<0.05);
        meancorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1))).*double(FUniv<0.05);
        save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
        new = [{dir1},{file}, {ZONEid},{1} ];
        infonew = [infonew;new];
        
        file = [name,'_',labelnode,num2str(NPERM),'permutation mean G1','.mat'];
        matcorr = real(squeeze(nanmean(cb1,1)));
        meancorr = real(squeeze(nanmean(cb1,1)));
        save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
        new = [{dir1},{file}, {ZONEid},{1} ];
        infonew = [infonew;new];

        file = [name,'_',labelnode,num2str(NPERM),'permutation mean G2','.mat'];
        matcorr = real(squeeze(nanmean(pb1,1)));
        meancorr = real(squeeze(nanmean(pb1,1)));
        save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
        new = [{dir1},{file}, {ZONEid},{1} ];
        infonew = [infonew;new];
        
        
        file = [name,'_',labelnode,num2str(NPERM),'permutation N G1','.mat'];
        matcorr = ncb1;
        meancorr = ncb1;
        save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
        new = [{dir1},{file}, {ZONEid},{1} ];
        infonew = [infonew;new];
       
 
        file = [name,'_',labelnode,num2str(NPERM),'permutation N G2','.mat'];
        matcorr = npb1;
        meancorr = npb1;
        save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
        new = [{dir1},{file}, {ZONEid},{1} ];
        infonew = [infonew;new];
        
        copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid))
    if ismac
        % Code to run on Mac platform problem with xlswrite
        [filepath,name,ext] = fileparts(xlslistfile);
        writetxtfile(fullfile(dir1,[name,'_',labelnode,num2str(NPERM),'PermutationTtest.txt']),infonew);
        disp(['Result .txt file saved: ' fullfile(dir1,[name,'_',labelnode,num2str(NPERM),'PermutationTtest.txt'])])

    else
        [filepath,name,ext] = fileparts(xlslistfile);
        xlswrite(fullfile(dir1,[name,'_',labelnode,num2str(NPERM),'PermutationTtest.xlsx']),infonew);
        disp(['Result .xlsx file saved ' fullfile(dir1,[name,'_',labelnode,num2str(NPERM),'PermutationTtest.xlsx'])])

    end
    
elseif isfield(job.c_statmatrix,'b_exportNBSformat')%
    dir1 = job.e_statmatrixPath{1};    
    if ~isdir(dir1)
        mkdir(dir1)
    end

    Mat = permute(MATall,[3,2,1])
    save(fullfile(dir1,'DATA.mat'),'Mat')
    design= double([ groupeall==1,groupeall==2])
    save(fullfile(dir1,'design.mat'),'design')
    contrast= [ 1,-1]
    save(fullfile(dir1,'G1higherG2contrast.mat'),'contrast','-mat')

listelectrode = DATA{1}.ZoneList
fid = fopen(fullfile(dir1,'nodeLabels.txt'),'w')
for i=1:numel(listelectrode)
    fprintf(fid,'%s',listelectrode{i})
    if i<numel(listelectrode)
        fprintf(fid,'\r')
    end
end
fclose(fid)
%OUTOUT POUR NBSToolbox not realisitic position... to be built if needed
fid = fopen(fullfile(dir1,'COG.txt'),'w')
for i=1:numel(listelectrode)
    fprintf(fid,'%2.3f %2.3f %2.3f\r',0,0,0 )
end
fclose(fid)
elseif isfield(job.c_statmatrix,'b_PearsonCorr_Mat')
    dir1 = job.e_statmatrixPath{1};
     infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    covariableall=[];
    covariablestring = job.c_statmatrix.b_PearsonCorr_Mat.b_Covariable_Mat;
    [token,remain] =strtok(covariablestring,',');
    covariableall =  [covariableall,{token}];
    while ~isempty(remain)        
        [token,remain] =strtok(remain,',');
        covariableall =  [covariableall,{token}];
    end
    
    for icov = 1:numel(covariableall)
        
     Pearsony = covariableall{icov};
     ycol = 0; 
     for icol=1:size(info,2)  
         if ~isnan(info{1,icol})
            if strcmp(strtrim(upper(deblank(info{1,icol}))), strtrim(upper(Pearsony)))
                ycol = icol;
            end   
         end
     end
     if ycol
         id= 1;
         for i=2:size(info,1)
            score(id,1) = info{i,ycol };
            id = id+1;
         end
         notfoundstophere = 1;
     else
         disp([Pearsony,' column not found'])
         notfoundstophere = 0;
     end 
     if notfoundstophere
     PearsonCoef = zeros(size(MATall,2),size(MATall,2));
     PearsonCoefSig = zeros(size(MATall,2),size(MATall,2));
    for i=1:size(MATall,2)
        for j=1:size(MATall,2)
            iduse = find(~isnan(score)& ~isnan(MATall(:,i,j)));
            X = MATall(iduse,i,j);
            y = score(iduse,1);
            if ~isempty(iduse)
             [rho,pval] = corr(X, y);
             PearsonCoef(i,j)=rho;
             PearsonCoefSig(i,j)=rho*double(pval< 0.05);
            end            
        end        
    end
    
    %WRITE IN A NEW FILE
      
        file = [name,'_',Pearsony,'PEARSON','.mat'];
        matcorr = PearsonCoef;
        meancorr = PearsonCoef;
        save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
        new = [{dir1},{file}, {ZONEid},{1} ];
        infonew = [infonew;new];
        

        file = [name,'_',Pearsony,'PEARSONp05','.mat'];
        matcorr = PearsonCoefSig;
        meancorr = PearsonCoefSig;
        save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
        new = [{dir1},{file}, {ZONEid},{1} ];
        infonew = [infonew;new];     
     end
    end
        copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid));
    if ismac
        % Code to run on Mac platform problem with xlswrite
        [filepath,name,ext] = fileparts(xlslistfile);
        writetxtfile(fullfile(dir1,['PearsonCorrelation.txt']),infonew);
        disp(['Result .txt file saved: ' fullfile(dir1,['PearsonCorrelation.txt'])])
    else
        [filepath,name,ext] = fileparts(xlslistfile);
        xlswrite(fullfile(dir1,['PearsonCorrelation.xlsx']),infonew);
        disp(['Result .xlsx file saved ' fullfile(dir1,['PearsonCorrelationPermutationTtest.xlsx'])])
    end
elseif isfield(job.c_statmatrix,'b_GLM_Mat')
    
    dir1 = job.e_statmatrixPath{1};
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    covariableall=[];
    covariablestring = job.c_statmatrix. b_GLM_Mat.b_Covariable_Mat;
    [token,remain] =strtok(covariablestring,',');
    covariableall =  [covariableall,{token}];
    while ~isempty(remain)        
        [token,remain] =strtok(remain,',');
        covariableall =  [covariableall,{token}];
    end
    
    for icov = 1:numel(covariableall)
     notfoundstophere = 1;
     Pearsony = covariableall{icov};
     ycol = 0; 
     for icol=1:size(info,2)  
         if ~isnan(info{1,icol})
            if strcmp(strtrim(upper(deblank(info{1,icol}))), strtrim(upper(Pearsony)))
                ycol = icol;
            end   
         end
     end
     if ycol
         id= 1;
         for i=2:size(info,1)
            score(id,icov) = info{i,ycol };
            id = id+1;
         end
       
     else
         disp([Pearsony,' column not found, no regression could be compute'])
         out='Stat tests';
         return
     end 
    end
    
      
        for icov = 1:numel(covariableall)
            eval(['bCOV',num2str(icov),' = zeros(size(MATall,2),size(MATall,2));']);
            eval(['bCOV',num2str(icov),'sig = zeros(size(MATall,2),size(MATall,2));']);
        end
  
       for i=1:size(MATall,2)
        for j=1:size(MATall,2)
            iduse = find(sum(~isnan(score),2)==size(score,2)& ~isnan(MATall(:,i,j)));
            X = score(iduse,:);
            y = MATall(iduse,i,j);
            if ~isempty(iduse)
                %R2 statistic, the F-statistic and its p-value, and an estimate of the error variance.
                [b,bint,r,rint,stats] = regress(y,X);
                try
                for icov = 1:numel(covariableall)
                 eval(['bCOV',num2str(icov),'(',num2str(i),',',num2str(j),')=',num2str(b(icov)),';'])
                end
                catch
                    1
                end
                if stats(3)<0.05
                    for icov = 1:numel(covariableall)
                        eval(['bCOV',num2str(icov),'sig(',num2str(i),',',num2str(j),')=',num2str(b(icov)),';'])
                    end
                end
            end            
        end        
    end
    
    %WRITE IN A NEW FILE
      for icov = 1:numel(covariableall)
        file = [name,'_',covariableall{icov},'.mat'];
         eval(['matcorr =','bCOV',num2str(icov),';'])
         eval(['meancorr =','bCOV',num2str(icov),';'])
         save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
         new = [{dir1},{file}, {ZONEid},{1} ];
         infonew = [infonew;new];
          file = [name,'_',covariableall{icov},'p05.mat'];
         eval(['matcorr =','bCOV',num2str(icov),'sig;'])
         eval(['meancorr =','bCOV',num2str(icov),'sig;'])
         save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
         new = [{dir1},{file}, {ZONEid},{1} ];
         infonew = [infonew;new];
      end
             
        copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid));
    if ismac
        % Code to run on Mac platform problem with xlswrite
        [filepath,name,ext] = fileparts(xlslistfile);
        writetxtfile(fullfile(dir1,['GLM.txt']),infonew);
        disp(['Result .txt file saved: ' fullfile(dir1,['GLM.txt'])])
    else
        [filepath,name,ext] = fileparts(xlslistfile);
        xlswrite(fullfile(dir1,['GLM.xlsx']),infonew);
        disp(['Result .xlsx file saved ' fullfile(dir1,['GLM.xlsx'])])
    end
    
end
out='Stat tests';

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
