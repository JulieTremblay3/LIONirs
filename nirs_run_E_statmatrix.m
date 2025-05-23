function out = nirs_run_E_statmatrix(job)
%Utilise les matrices cr��es pour calcul�s les statistiques.
%Fichier xls avec colonne 1 path, colonne 2 nom de la matrices, colonne 3
%la zone ou le channel list
xlslistfile = job.f_matrix{1};
[~,~,ext] =fileparts(xlslistfile);
try
    if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
        try
            [raw, txt, info]=xlsread(xlslistfile);
        catch
            disp(['Error could not read file: ',xlslistfile]);

            disp(['Verify the file location or if the file is already open']);
        end
    elseif strcmp(ext,'.txt')
        [num, txt, info] = readtxtfile_asxlsread(xlslistfile);
    end
catch
    disp(['Error could not read file: ',xlslistfile]);
    disp(['Verify the file location or if the file is already open']);
    return
end

%TRY TO REMOVE nan SUJBECT AT THE END
iremoverow = [];
for i = 1: size(info,1);
    if isnan(info{i,1});
        iremoverow= [iremoverow,i];
    end
end
info(iremoverow,:)=[];
%Load the matrices same as any
groupeall = [];
% try
%     [filepath,name,ext] = fileparts(xlslistfile);
%     jobreal = job;
%     disp(['LOAD ',filepath, filesep,'WORKSPACE.mat'])
%     load(fullfile(filepath, 'WORKSPACE.mat'))
%     job = jobreal;
% catch

for isubject=2:size(info,1)
    id = isubject-1;
    try
        tmp = info{isubject,2};
        if strcmp(tmp(end-2:end),'mat')
            try
                MAT = load(fullfile(info{isubject,1},[ info{isubject,2}]));
                disp(['Load ',fullfile(info{isubject,1},info{isubject,2}), ' GR',num2str(info{isubject,4})])
            catch
                disp(['Error Load ',fullfile(info{isubject,1},[ info{isubject,2},'.mat'])]);
                MAT.meancorr = nan(size(DATA{1,1}.MAT));
                MAT.matcorr =  MAT.meancorr;
                MAT.ZoneList = DATA{1,1}.ZoneList;
            end
        else
            try
                MAT = load(fullfile(info{isubject,1},[ info{isubject,2},'.mat']));
                disp(['Load ',fullfile(info{isubject,1},[ info{isubject,2},'.mat G',num2str(info{isubject,4})])]);
            catch
                disp(['Error Load ',fullfile(info{isubject,1},[ info{isubject,2},'.mat G',num2str(info{isubject,4})])]);
                MAT.meancorr = nan(size(DATA{1,1}.MAT));
                MAT.matcorr =  MAT.meancorr;
                MAT.ZoneList = DATA{1,1}.ZoneList;
            end
        end

        DATA{id}.ZoneList = MAT.ZoneList;
        if isfield(MAT, 'meancorr')
            matcorr = MAT.meancorr;
            if job.m_fishertransform==1 % add fisher transform
                matcorr =  1/2*log((1+matcorr)./(1-matcorr));
                if ~isreal(matcorr)
                    disp(['WARNING: Verify if you apply Fisher transform twice, the output matrix is complex.'])
                end
                DATA{id}.MAT = matcorr;

            elseif job.m_fishertransform==2 % no fisher transform
                DATA{id}.MAT = matcorr;
            elseif job.m_fishertransform==3
                matcorr =  abs(1/2*log((1+matcorr)./(1-matcorr)));
                DATA{id}.MAT = matcorr;
            end
        end

        DATA{id}.name = info{isubject,2};
        DATA{id}.GR = info{isubject,4};
        infocov = [];
        if size(info,2) > 4
            %Define Additional covariable
            idcov = 1 ;
            for icolumn = 5:size(info,2)
                if ~ischar(info{isubject,icolumn})
                    eval( ['DATA{id}.cov',num2str(idcov),' = ', num2str(info{isubject,icolumn}),';']);
                    idcov=idcov+1;
                    if isubject==2
                        infocov = [infocov ;info(1,icolumn)];
                    end
                end
            end
        end


        DATA{id}.System = 'ISS';
        if 1
            DATA{id}.System = 'EEG';
        end
        try
            load(info{isubject,3},'-mat');  %load zone specific folder
        catch
            load(fullfile(info{isubject,1}, info{isubject,3}),'-mat');  %load zone default current directory
        end

        names = fieldnames(zone); 
        for iname = 1:numel(names)
            eval(['DATA{id}.zone.',names{iname},' =zone.',names{iname},';']);
        end
        list_subject{id} =DATA{id}.name;
        groupeall = [groupeall; info{isubject,4}];
        clear MAT
    catch em
        disp(['Error with line', num2str(isubject)])
        disp(em.identifier)
        disp(em.message)
        groupeall = [groupeall; info{isubject,4}];
    end
end
%end %load already existing export
alpha_threshold = job.e_statcomponent_alpha;
if job.m_nodeunit==1 %channel mode
    idsubject = 1:numel(groupeall);
    for isubject = 1:numel(groupeall)
        try
            MATall(isubject,:,:)=DATA{isubject}.MAT;
            groupid(isubject)= DATA{idsubject(isubject)}.GR;
        catch
            MATall(isubject,:,:)=nan;
            groupid(isubject)= 0;
        end
    end
    [~,tmpname,~]=fileparts(info{isubject,3});
    ZONEid = [tmpname,'.zone'];
    ZoneList =  DATA{end}.ZoneList;
    labelnode = 'c';
    dir1 = job.e_statmatrixPath{1};
elseif  job.m_nodeunit==2 %use by zone

    MATall =zeros(numel(DATA),numel(DATA{id}.zone.label),numel(DATA{id}.zone.label));
    for isubject = 1:numel(groupeall)
        try
            List = DATA{isubject}.ZoneList;
            name = List{1};
            if strcmp(name(1:2),'D0')
                DATA{isubject}.System = 'NIRx';

            elseif strcmp(name(1:1),'E')
                DATA{isubject}.System = 'EEG';
            else
                DATA{isubject}.System = 'ISS';

            end
            for izone = 1:numel(DATA{isubject}.zone.label)
                ML = DATA{isubject}.zone.ml;
                DATA{isubject}.zone.plotLst;
                idlisti = [];
                idliststr = [];
                chzone = DATA{isubject}.zone.plotLst{izone};
                for ichzone = 1:numel(chzone);
                    ich = chzone(ichzone);
                    switch DATA{isubject}.System
                        case 'ISS Imagent'
                            strDet = SDDet2strboxy_ISS(ML(ich,2));
                            strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                            idch = strmatch([strDet, ' ',strSrs ],List,'exact');
                        case 'NIRx'
                            strDet = SDDet2strboxy(ML(ich,2));
                            strSrs = SDPairs2strboxy(ML(ich,1));
                            idch = strmatch([strDet, ' ',strSrs ],List,'exact');
                        case 'EEG'
                            strDet = 'EEG';
                            strSrs = 'EEG';
                            idch = ich;
                    end
                    idliststr =[idliststr,{[strDet, ' ',strSrs ]}];
                    idlisti = [idlisti, idch];
                end

                for jzone = 1:numel(DATA{isubject}.zone.label)
                    idliststr =  [];
                    idlistj = [];
                    chzone = DATA{isubject}.zone.plotLst{jzone};
                    for ichzone = 1:numel(chzone);
                        ich = chzone(ichzone);
                        switch DATA{isubject}.System
                            case 'ISS Imagent'
                                strDet = SDDet2strboxy_ISS(ML(ich,2));
                                strSrs = SDPairs2strboxy_ISS(ML(ich,1));
                                idch = strmatch([strDet, ' ',strSrs ],List,'exact');
                            case 'NIRx'
                                strDet = SDDet2strboxy(ML(ich,2));
                                strSrs = SDPairs2strboxy(ML(ich,1));
                                idch = strmatch([strDet, ' ',strSrs ],List,'exact');
                            case 'EEG'
                                strDet = 'EEG';
                                strSrs = 'EEG';
                                idch = ich; %no verification take channel as it is
                        end
                        idliststr =[idliststr,{[strDet, ' ',strSrs ]}];
                        idlistj = [idlistj, idch];
                    end
                    matROI = DATA{isubject}.MAT(idlisti,idlistj);


                    id = find(matROI==0);
                    if ~isempty(id)
                        matROI(id)=nan;
                    end
                    id = find(matROI==Inf);
                    if ~isempty(id)
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
        catch
            MATall(isubject,:,:)=nan;
            groupid(isubject)= 0;
            disp(['Warning subject:', num2str(isubject), ' do not find zone it be set as nan.'])
        end
        ZoneLabel =  DATA{end}.zone.label;
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
    dir1 = job.e_statmatrixPath{1};
    [~,tmpname,~]=fileparts(info{isubject,3});
    try
        save(fullfile(dir1,['avg', tmpname,'.zone']),'zone','-mat');
    catch
        if ~isdir(dir1)
            mkdir(dir1);
            disp(['Folder ''', dir1, ''' created']);
        end
        save(fullfile(dir1,['avg', tmpname,'.zone']),'zone','-mat');
        disp(['Save ', fullfile(dir1,['avg', tmpname,'.zone'])]);
    end
    ZONEid = ['avg', tmpname,'.zone'];
    for izone = 1:numel(DATA{isubject}.zone.plotLst) %Utilise la coordonn�e moyenne pour chaque region
        idlisti = DATA{isubject}.zone.plotLst{izone};
        chanposzone(izone,:) = mean(zone.pos(idlisti,:));
    end

end

if ~isdir(dir1)
    mkdir(dir1)
    disp(['Create new folder to save results : ', dir1] )
end
matid = zeros(size(MATall,2),size(MATall,2));
id =1;
%modifier pour inclure la diagonal pour gerer le cas ou les zones sont
%utilis�es
for ielex=1:size(MATall,2)
    ielex;
    ieley = 1;
    while ieley <= ielex
        matid(ielex,ieley)=id;
        ieley = ieley + 1;
        id = id + 1;
    end
end
idhalf = find(matid);
matid(idhalf);

try

    if job.m_nodeunit==1
        sujet = 1; %need projet V5
        zonelist = DATA{sujet}.ZoneList;
        chanpos = zeros(numel(zonelist),3);
        for id=1:numel(zonelist)
            name = zonelist{id};
            [DetL,SrsL] =strtok(name,' ');
            SDdetL = StrBoxy2SDDet(DetL);
            SDsrsL =str2num(SrsL(3:end));
            idml = find(DATA{sujet}.zone.ml(:,1)==  SDsrsL & DATA{sujet}.zone.ml(:,2)==  SDdetL & DATA{sujet}.zone.ml(:,4)== 1);
            if ~isempty(idml)
                chanpos(id,:) =  DATA{1}.zone.pos(idml,1:3);
            end
        end
        disp(['Use channel position'])
    elseif    job.m_nodeunit==2 %pour les zones utilis�s les zones moyenn�es
        chanpos= chanposzone;
        disp('Use zone middle point  position')
        %figure;plot3(chanpos(:,1),chanpos(:,2),chanpos(:,3),'x')
    end
catch

    disp('Multiple comparaison cluster could not find channel position, use all cluster')
end

[filepath,name,ext] = fileparts(xlslistfile);
if  isfield(job.c_statmatrix,'m_export_matrix')
    if job.c_statmatrix.m_export_matrix == 1
        dir1 = job.e_statmatrixPath{1}
        save(fullfile(dir1,'Exportmatrice.mat'),'MATall','info','ZoneList', '-mat')
        disp(['Save export: ', fullfile(dir1,'Exportmatrice.mat')])
        save(fullfile(dir1,'WORKSPACE.mat'))
    else
        disp('xls export not available yet')
    end
elseif isfield(job.c_statmatrix,'b_TtestOneSamplematrix')
    AllC = [];
    id =1;
    %Use one or more specific groupe
    GRname = ['GR',num2str(job.c_statmatrix.b_TtestOneSamplematrix.e_TtestOneSampleGR),' '];

    idG1 = find(  sum(groupeall==job.c_statmatrix.b_TtestOneSamplematrix.e_TtestOneSampleGR,2));
    MATallG1 = MATall( idG1,:,:);
    meanall = squeeze(nanmean( MATallG1,1));
    refval  = job.c_statmatrix.b_TtestOneSamplematrix.e_TtestOneSample_meanvalue;
    tval = squeeze((nanmean(MATallG1(:,:,:),1)-refval)./(nanstd(MATallG1(:,:,:),1)./sqrt(sum(~isnan(MATallG1(:,:,:)),1))  ));
    dfall = squeeze(sum(~isnan(MATallG1(:,:,:)),1))-1;
    disp(['One sample ttest vs hypothetical mean ', num2str( job.c_statmatrix.b_TtestOneSamplematrix.e_TtestOneSample_meanvalue)])

    for i=1:size(tval,1)
        for j=1:size(tval,2)
            try
                % Compute the correct p-value for the test, and confidence intervals
                % if requested.
                if job.c_statmatrix.b_TtestOneSamplematrix.m_TtestOneSample_matrix == 1; % two-tailed test
                    pval(i,j) = 2 * tcdf(-abs(-tval(i,j)), dfall(i,j));
                    %     if nargout > 2
                    %         crit = tinv((1 - alpha / 2), df) .* ser;
                    %         ci = cat(dim, xmean - crit, xmean + crit);
                    %     end
                elseif job.c_statmatrix.b_TtestOneSamplematrix.m_TtestOneSample_matrix == 3; % right one-tailed test
                    pval(i,j) = tcdf(-tval(i,j), dfall(i,j));
                    %     if nargout > 2
                    %         crit = tinv(1 - alpha, df) .* ser;
                    %         ci = cat(dim, xmean - crit, Inf(size(p)));
                    %     end
                elseif job.c_statmatrix.b_TtestOneSamplematrix.m_TtestOneSample_matrix == 2; % left one-tailed test
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
        mkdir(dir1);
    end
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];

    try %SAVE Zone individual file
        for ifile=1:numel(idG1)
            if job.m_nodeunit==1 %channel
                file = [info{idG1(ifile)+1,2},'c','.mat'];
            else job.m_nodeunit==2 %zon
                file = [info{idG1(ifile)+1,2},'z','.mat'];
            end

            meancorr =  squeeze(MATall(idG1(ifile),:,:));
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{1} ];
            infonew = [infonew;new];
        end
    catch
    end





    meancorr = meanall;
    totaltrialgood = mean(dfall(:));
    file = [name,labelnode,GRname, 'OneSampleTtest mean','.mat'];
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];

    file = [name,labelnode,GRname, 'OneSampleTtest meanp',num2str(alpha_threshold),'.mat'];
    meancorr = meanall.*double(pval<alpha_threshold);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];

    %find fdr
    try
        [FDR,Q] = mafdr(pval(:));
        Q = reshape(Q,size(pval));
        file = [name,labelnode,GRname, 'OneSampleTtest meanFDR',num2str(alpha_threshold),'.mat'];
        meancorr = meanall.*double(Q<alpha_threshold);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    catch
        disp('Warning: could not perform FDR correction')
    end

    % ZoneList = MAT.ZoneList;

    meancorr = tval;
    totaltrialgood = mean(dfall(:));
    file = [name,labelnode,GRname, 'OneSampleTtest tval','.mat'];
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];

    file = [name,labelnode,GRname,'OneSampleTtest tvalp',num2str(alpha_threshold),'.mat'];
    meancorr = tval.*double(pval<alpha_threshold);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];

    try
        file = [name,labelnode,GRname,'OneSampleTtest tvalFDR',num2str(alpha_threshold),'.mat'];
        meancorr = tval.*double(Q<alpha_threshold);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    catch
        disp(['File : ', file ' could not be create'])
    end


    try
        file = [name,labelnode,GRname,'OneSampleTtest N.mat'];
        meancorr = squeeze(sum(~isnan(MATallG1(:,:,:)),1));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    catch
        disp(['File : ', file ' could not be create'])
    end


    try
        if ~strcmp(fullfile(info{isubject,1}, ZONEid),fullfile(dir1,  ZONEid))
            copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid));
        end
    catch
    end
    %  dir1 = job.e_statmatrixPath{1};
    try
        if ismac
            % Code to run on Mac platform problem with xlswrite
            writetxtfile(fullfile(dir1,[name,labelnode,GRname,'SimpleTtest.txt']),infonew);
            disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'SimpleTtest.txt'])]);
        else
            try
                xlswrite(fullfile(dir1,[name,labelnode,GRname,'SimpleTtest.xlsx']),infonew);
                disp(['Result .xlsx file saved: ' fullfile(dir1,[name,labelnode,'SimpleTtest.xlsx'])]);
            catch
                writetxtfile(fullfile(dir1,[name,labelnode,'SimpleTtest.txt']),infonew);
                disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,GRname,'SimpleTtest.txt'])]);
            end
        end
    catch
        disp(['Error could not save .xlsx file: ' fullfile(dir1,[name,labelnode,GRname 'SimpleTtest.xlsx'])]);
    end
elseif isfield(job.c_statmatrix,'b_UnpairedTtest')
    dir1 = job.e_statmatrixPath{1};
    if ~isdir(dir1)
        mkdir(dir1);
    end
    %WRITE IN A NEW FILE
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    GRname = ['GR ',num2str(job.c_statmatrix.b_UnpairedTtest.e_TtestOneSampleGR),' vs ',num2str(job.c_statmatrix.b_UnpairedTtest.e_TtestOneSampleGR2),' '];
    g1 = find(sum(groupeall==job.c_statmatrix.b_UnpairedTtest.e_TtestOneSampleGR,2));
    g2 = find(sum(groupeall==job.c_statmatrix.b_UnpairedTtest.e_TtestOneSampleGR2,2));
    iduse = [g1;g2];
    halfMAT = MATall(:,idhalf);
    if isfield(job.c_statmatrix.b_UnpairedTtest.c_statpermutation,'b_Nopermutation')
        for ilink = 1:size(halfMAT,2)
            [h,p,ci,stats] = ttest2(halfMAT(g1,ilink), halfMAT(g2,ilink));
            pval(:,ilink) = p;
            tval(:,ilink) = stats.tstat;
            try; cohend(:,ilink)=computeCohen_d(halfMAT(g1,ilink), halfMAT(g2,ilink),'independent');catch;end

        end
    elseif isfield(job.c_statmatrix.b_UnpairedTtest.c_statpermutation,'b_permutation')
        halfMAT = MATall(:,idhalf);
        g1 = find(sum(groupeall==job.c_statmatrix.b_UnpairedTtest.e_TtestOneSampleGR,2));
        g2 = find(sum(groupeall==job.c_statmatrix.b_UnpairedTtest.e_TtestOneSampleGR2,2));
        iduse = [g1;g2];
        nperm = str2num(job.c_statmatrix.b_UnpairedTtest.c_statpermutation.b_permutation.e_npermutation);
        for iperm=1:nperm
            permfix(iperm,:) =  randperm(numel(iduse));
        end
        disp(['Running ', num2str(nperm),' permutations on ', num2str(size(halfMAT,2)),' link :'])
        %loopslow to be improve
        for ilink = 1:size(halfMAT,2)
            fprintf('%s ',num2str(ilink));
            if mod(ilink,30)==0
                fprintf('\n ');
            end
            try
                Yval = halfMAT(iduse,ilink);
                [h,p,ci,stats] = ttest2(halfMAT(g1,ilink), halfMAT(g2,ilink));
                pval(:,ilink) = p;
                tval(:,ilink) = stats.tstat;
                try; cohend(:,ilink)=computeCohen_d(halfMAT(g1,ilink), halfMAT(g2,ilink),'independent');catch;end

                for idperm=1:nperm
                    Yval = halfMAT( iduse,ilink);
                    GR = permfix(idperm,:);
                    [h,p,ci,stats] = ttest2(Yval(GR(1:numel(g1))), Yval(GR((numel(g1)+1):end)));
                    tvalperm(:,idperm,ilink) =stats.tstat;
                end
            catch
                fprintf('error link ');
            end
        end

        t_dist_all =  tvalperm ;
        [fecdf,xecdf] = ecdf(t_dist_all(:));
        clustercritical=xecdf(sum(fecdf<(1-alpha_threshold)));
        disp(['critical for clustering p=', num2str(alpha_threshold),' T=',num2str(clustercritical)])
        %  figure;hist(t_dist_all(:),1000)
        if isfield(job.c_statmatrix.b_UnpairedTtest.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting,'b_MCT_ClusterBased')
            statrand = permute(tvalperm(1,:,:), [3 2 1]);
            statobs = permute(tval(1,:),[2,1]);
            tmp = randperm(numel(statrand));
            statrand =reshape(statrand(tmp),size(statrand,1),size(statrand,2));
            minnbchan = job.c_statmatrix.b_UnpairedTtest.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_minnbchan;
            neighbourdist = job.c_statmatrix.b_UnpairedTtest.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_neighbourdist;
            disp(['critical for clustering p=', num2str(alpha_threshold),' T=',num2str(clustercritical)])
            [stat, matneig] = FindClusterBasedPermutationInMatrix(chanpos, neighbourdist,clustercritical, statobs, statrand ,  minnbchan );
            %  [stat, matneig] = FindClusterBasedPermutationInMatrix(chanpos, 'distance',clustercritical, statobs, statrand );
            save(fullfile(dir1,'Unpairedttest_stat.mat'),'stat')
            disp(['Save cluster', fullfile(dir1,'Unpairedttest_stat.mat') ])

            % %View link for neighbor debug
            % if 0
            %     for ilink=1:300
            %         file = [name,labelnode,GRname, num2str(nperm),' neiglink',num2str(ilink),'.mat'];
            %         matcorr = zeros(size(MATall,2),size(MATall,2));
            %         matcorr(idhalf)=matneig(:,ilink);
            %         matcorr = matcorr %+flipud(rot90(matcorr));
            %         meancorr = matcorr;
            %         save(fullfile(dir1,[file]),'ZoneList','matcorr','meancorr');
            %         disp(['Save: ', fullfile(dir1,[file])]);
            %         new = [{dir1},{file}, {ZONEid},{0} ];
            %         infonew = [infonew;new];
            %     end
            % end

            if isfield(stat,'posclusterslabelmat')
                for iposcluster = 1:numel(stat.posclusters)
                    file = [name,labelnode,GRname, num2str(nperm),' posCluster', num2str(iposcluster),' MCp=',num2str(stat.posclusters(iposcluster).prob),'.mat'];
                    tmp=(stat.posclusterslabelmat==iposcluster);
                    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                    new = [{dir1},{file}, {ZONEid},{0} ];
                    infonew = [infonew;new];

                    file = [name,labelnode,GRname, num2str(nperm),' posCluster tval', num2str(iposcluster),' MCp=',num2str(stat.posclusters(iposcluster).prob),'.mat'];
                    tmp=squeeze(tval(:)).*(stat.posclusterslabelmat==iposcluster);
                    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                    new = [{dir1},{file}, {ZONEid},{0} ];
                    infonew = [infonew;new];

                end
            else
                disp('No positive cluster found')
            end
            if isfield(stat,'negclusterslabelmat') %technicaly no neg cluster in f value...
                for inegcluster = 1:numel(stat.negclusters)
                    file = [name,labelnode,GRname, num2str(nperm),' negCluster', num2str(inegcluster),' MCp=',num2str(stat.negclusters(inegcluster).prob),'.mat'];

                    tmp=(stat.negclusterslabelmat==inegcluster);
                    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                    new = [{dir1},{file}, {ZONEid},{0} ];
                    infonew = [infonew;new];

                    file = [name,labelnode, GRname, num2str(nperm),' negCluster tval', num2str(inegcluster),' MCp=',num2str(stat.negclusters(inegcluster).prob),'.mat'];
                    tmp=squeeze(tval(:)).*(stat.negclusterslabelmat==inegcluster);
                    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                    new = [{dir1},{file}, {ZONEid},{0} ];
                    infonew = [infonew;new];
                end
            else
                disp('No negative cluster found')
            end

        end
    end
    try %SAVE Zone individual file
        for ifile=1:numel(g1)
            if job.m_nodeunit==1 %channel
                file = ['G1', info{g1(ifile)+1,2},'c.mat'];

            elseif job.m_nodeunit==2%zone
                file = ['G1', info{g1(ifile)+1,2},'z.mat'];
            end
            meancorr =  squeeze(MATall(g1(ifile),:,:));
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{1} ];
            infonew = [infonew;new];
        end
        for ifile = 1:numel(g2)
            if job.m_nodeunit==1 %channel
                file = ['G2', info{g2(ifile)+1,2},'c.mat'];
            elseif job.m_nodeunit==2 %zone
                file = ['G2', info{g2(ifile)+1,2},'z.mat'];
            end

            meancorr =squeeze(MATall(g2(ifile),:,:));
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{2} ];
            infonew = [infonew;new];
        end
    catch
    end


    try %optionnal cohend
        file = [name,labelnode,GRname, 'Cohend_TPositive','.mat'];
        tmp=squeeze(cohend(:)).*squeeze(tval(:)>0);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        file = [name,labelnode,GRname,'Cohend_TNegative','.mat'];
        tmp=squeeze(cohend(:)).*squeeze(tval(:)<0);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        file = [name,labelnode,GRname,'Cohend','.mat'];

        tmp = squeeze(cohend(:));
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

    catch
        disp(['Error save: ', fullfile(dir1,[file])]);
        disp('Error computeCohen_d')
        disp('Cohen_d function must be install in the set path see:  Ruggero G. Bettinardi (2022). computeCohen_d(x1, x2, varargin) (https://www.mathworks.com/matlabcentral/fileexchange/62957-computecohen_d-x1-x2-varargin), MATLAB Central File Exchange. Retrieved November 16, 2022. ')
        p = mfilename('fullpath');
        [pathname,filename,ext]=fileparts(p);
        disp(['Add the fonction in matlab path: addpath(''',pathname, '\External\computeCohen_d'')'])
        disp('And verify if Fisher transform is not apply twice, data must be real')
    end


    try
        file = [name,labelnode,GRname,'UnpairedTtest tval','.mat'];
        meancorr = halfmat2mat(tval,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        file = [name,labelnode,GRname,'UnpairedTtest tval_p', num2str(alpha_threshold),'unc','.mat'];
        tmp=squeeze(tval(:)).*squeeze(double(pval(:)<alpha_threshold));
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        file = [name,labelnode,GRname,'UnpairedTtest -tval','.mat'];
        tmp=-squeeze(tval(:));
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    catch
        disp(['File : ', file ' could not be create'])
    end

    try
        file = [name,labelnode,GRname,'UnpairedTtest pval','.mat'];
        tmp =squeeze(pval(:));
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    catch
        disp(['File : ', file ' could not be create'])
    end
    try
        if ~strcmp(fullfile(info{isubject,1}, ZONEid),fullfile(dir1,  ZONEid))
            copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid));
        end
    catch
    end
    %  dir1 = job.e_statmatrixPath{1};
    try
        if ismac
            % Code to run on Mac platform problem with xlswrite
            writetxtfile(fullfile(dir1,[name,labelnode,GRname,'UnpairedTtest.txt']),infonew);
            disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,GRname,'UnpairedTtest.txt'])]);
        else
            try
                xlswrite(fullfile(dir1,[name,labelnode,GRname,'UnpairedTtest.xlsx']),infonew);
                disp(['Result .xlsx file saved: ' fullfile(dir1,[name,labelnode,GRname,'UnpairedTtest.xlsx'])]);
            catch
                writetxtfile(fullfile(dir1,[name,labelnode,GRname,'UnpairedTtest.txt']),infonew);
                disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,GRname,'UnpairedTtest.txt'])]);
            end
        end
    catch
        disp(['Error could not save .xlsx file: ' fullfile(dir1,[name,labelnode,GRname 'UnpairedTtest.xlsx'])]);
    end

elseif isfield(job.c_statmatrix,'b_PermutationTest')
    ESTAD = 4;
    INDEP = 1;
    NPERM = str2num(job.c_statmatrix.b_PermutationTest.e_npermutation);
    g1 = find(sum(groupeall==job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR,2));
    g2 = find(sum(groupeall==job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR2,2));
    GRname = ['GR ',num2str(job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR),' vs ',num2str(job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR2),' '];
    % idG1 = find(  sum(groupeall==job.c_statmatrix.b_TtestOneSamplematrix.e_TtestOneSampleGR,2));

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
        %keep home made calculation TestPermut2Grupos give tval in absolute
        %value
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
                    %if Toij(i,j)==0
                    [FSupSup,FSupDeriv,FSupTime,pij,tij] = TestPermut2Grupos(ESTAD,INDEP,tmpcb1,tmppb1,NPERM);
                    %apply symetric
                    ncb1(i,j) = numel(tmpcb1);
                    npb1(i,j) = numel(tmppb1);
                    FUniv(i,j) = pij;
                    %   Toij(i,j) = tij;
                    ncb1(j,i) = numel(tmpcb1);
                    npb1(j,i) = numel(tmppb1);
                    FUniv(j,i) = pij;
                    %  Toij(j,i) = tij;

                    % end
                end
            end
        else
            [FSupSup,FSupDeriv,FSupTime,FUniv,Toij] = TestPermut2Grupos(ESTAD,INDEP,cb1,pb1,NPERM);
            ncb1 = ones(size(cb1,2),size(cb1,3))*numel(size(cb1,1));
            npb1 = ones(size(pb1,2),size(pb1,3))*numel(size(pb1,1));
        end
    end
    MeanG1 = squeeze(nanmean(cb1));
    MeanG2 = squeeze(nanmean(pb1));


    dir1 = job.e_statmatrixPath{1};
    if ~isdir(dir1)
        mkdir(dir1);
    end

    %WRITE IN A NEW FILE
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    file = [name,'_',labelnode,GRname,num2str(NPERM),'permutation tstat','.mat'];
    meancorr = Toij;
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ',fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];

    file = [name,'_',labelnode,GRname,num2str(NPERM),'permutation 1-pval','.mat'];
    meancorr = 1-FUniv;
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ',fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];


    file = [name,'_',labelnode,GRname,num2str(NPERM),'permutation mean', num2str(job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR),'-',num2str(job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR2) ,'.mat'];
    meancorr = real(squeeze((nanmean(cb1,1)- nanmean(pb1,1))));
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ',fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];


    file = [name,'_',labelnode,GRname,num2str(NPERM),'permutation mean ', num2str(job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR2),'-',num2str(job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR) ,'.mat'];
    meancorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1)));
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ',fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];

    file = [name,'_',labelnode,GRname,num2str(NPERM),'permutation mean ', num2str(job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR),'-',num2str(job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR2),' p',num2str(alpha_threshold),'.mat'];
    meancorr = real(squeeze(nanmean(cb1,1))-squeeze(nanmean(pb1,1))).*double(FUniv<alpha_threshold);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ',fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];

    file = [name,'_',labelnode,GRname,num2str(NPERM),'permutation mean ',num2str(job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR2),'-',num2str(job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR) ,' p',num2str(alpha_threshold),'.mat'];
    meancorr = real(squeeze(nanmean(pb1,1))-squeeze(nanmean(cb1,1))).*double(FUniv<alpha_threshold);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ',fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];

    file = [name,'_',labelnode,GRname,num2str(NPERM),'permutation mean G',num2str(job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR),'.mat'];
    meancorr = real(squeeze(nanmean(cb1,1)));
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ',fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];

    file = [name,'_',labelnode,GRname,num2str(NPERM),'permutation mean G',num2str(job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR2),'.mat'];
    meancorr = real(squeeze(nanmean(pb1,1)));
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ',fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];


    file = [name,'_',labelnode,GRname,num2str(NPERM),'permutation N G',num2str(job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR),'.mat'];
    meancorr = ncb1;
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ',fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];


    file = [name,'_',labelnode,GRname,num2str(NPERM),'permutation N G',num2str(job.c_statmatrix.b_PermutationTest.e_TtestOneSampleGR2),'.mat'];
    meancorr = npb1;
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ',fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{1} ];
    infonew = [infonew;new];
    try
        if ~strcmp(fullfile(info{isubject,1}, ZONEid),fullfile(dir1,  ZONEid))
            copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid))
        end
    catch
    end
    if ismac
        % Code to run on Mac platform problem with xlswrite
        [filepath,name,ext] = fileparts(xlslistfile);
        writetxtfile(fullfile(dir1,[name,GRname,labelnode,num2str(NPERM),'PermutationTtest.txt']),infonew);
        disp(['Result .txt file saved: ' fullfile(dir1,[name,'_',labelnode,num2str(NPERM),'PermutationTtest.txt'])]);

    else
        [filepath,name,ext] = fileparts(xlslistfile);
        try
            xlswrite(fullfile(dir1,[name,GRname,labelnode,num2str(NPERM),'PermutationTtest.xlsx']),infonew);
            disp(['Result .xlsx file saved ' fullfile(dir1,[name,'_',labelnode,num2str(NPERM),'PermutationTtest.xlsx'])]);
        catch
            writetxtfile(fullfile(dir1,[name,GRname,labelnode,num2str(NPERM),'PermutationTtest.txt']),infonew);
            disp(['Result .xlsx file saved ' fullfile(dir1,[name,'_',labelnode,num2str(NPERM),'PermutationTtest.txt'])]);
        end
    end

elseif isfield(job.c_statmatrix,'b_exportNBSformat')%
    dir1 = job.e_statmatrixPath{1};
    if ~isdir(dir1)
        mkdir(dir1);
    end

    Mat = permute(MATall,[3,2,1]);
    disp(['Save file to NBS: ', fullfile(dir1,'DATA.mat')]);
    save(fullfile(dir1,'DATA.mat'),'Mat');
    design= double([ groupeall==1,groupeall==2]);
    save(fullfile(dir1,'design.mat'),'design');
    contrast= [ 1,-1];
    save(fullfile(dir1,'G1higherG2contrast.mat'),'contrast','-mat');

    listelectrode = DATA{1}.ZoneList;
    fid = fopen(fullfile(dir1,'nodeLabels.txt'),'w');
    for i=1:numel(listelectrode)
        fprintf(fid,'%s',listelectrode{i});
        if i<numel(listelectrode)
            fprintf(fid,'\r');
        end
    end
    fclose(fid);
    %OUTOUT POUR NBSToolbox not realisitic position... to be built if needed
    fid = fopen(fullfile(dir1,'COG.txt'),'w');
    for i=1:numel(listelectrode)
        fprintf(fid,'%2.3f %2.3f %2.3f\r',0,0,0 );
    end
    fclose(fid);
elseif isfield(job.c_statmatrix,'b_PearsonCorr_Mat')
    dir1 = job.e_statmatrixPath{1};
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}, {job.c_statmatrix.b_PearsonCorr_Mat.b_Covariable_Mat}];
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
                        PearsonCoefSig(i,j)=rho*double(pval< alpha_threshold);
                    end
                end
            end
            iduse = find(~isnan(score))
            try %SAVE Zone individual file
                for ifile=1:numel(iduse)
                    if job.m_nodeunit==1 %channel
                        file = [ info{iduse(ifile)+1,2},'c.mat'];
                    elseif job.m_nodeunit==2%zone
                        file = [ info{iduse(ifile)+1,2},'z.mat'];
                    end
                    meancorr =  squeeze(MATall(iduse(ifile),:,:));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                    new = [{dir1},{file}, {ZONEid},{1},  {score(iduse(ifile))}];
                    infonew = [infonew;new];
                end
            catch
            end


            %WRITE IN A NEW FILE

            file = [name,'_',Pearsony,'PEARSON','.mat'];
            meancorr = PearsonCoef;
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{1},{nan} ];
            infonew = [infonew;new];


            file = [name,'_',Pearsony,'PEARSONp',num2str(alpha_threshold),'.mat'];
            meancorr = PearsonCoefSig;
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{1},{nan} ];
            infonew = [infonew;new];
        end
    end
    try
        copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid));
    catch
    end
    if ismac
        % Code to run on Mac platform problem with xlswrite
        [filepath,name,ext] = fileparts(xlslistfile);
        writetxtfile(fullfile(dir1,['PearsonCorrelation.txt']),infonew);
        disp(['Result .txt file saved: ' fullfile(dir1,['PearsonCorrelation.txt'])])
    else
        [filepath,name,ext] = fileparts(xlslistfile);
        try
            xlswrite(fullfile(dir1,['PearsonCorrelation.xlsx']),infonew);
            disp(['Result .xlsx file saved ' fullfile(dir1,['PearsonCorrelationPermutationTtest.xlsx'])]);
        catch
            writetxtfile(fullfile(dir1,['PearsonCorrelation.txt']),infonew);
            disp(['Result .xlsx file saved ' fullfile(dir1,['PearsonCorrelationPermutationTtest.xlsx'])]);
        end
    end
elseif isfield(job.c_statmatrix,'b_GLM_Mat')

    dir1 = job.e_statmatrixPath{1};

    covariableall=[];
    covariablestring = job.c_statmatrix. b_GLM_Mat.b_Covariable_Mat;
    [token,remain] =strtok(covariablestring,',');
    covariableall =  [covariableall,{strtrim(token)}];
    while ~isempty(remain)
        [token,remain] =strtok(remain,',');
        covariableall =  [covariableall,{strtrim(token)}];
    end

    for icov = 1:numel(covariableall)
        notfoundstophere = 1;
        Pearsony = covariableall{icov};
        ycol = 0;
        for icol=1:size(info,2)
            if ~isnan(info{1,icol})
                try
                    if strcmp(strtrim(upper(deblank(info{1,icol}))), strtrim(upper(Pearsony)));
                        ycol = icol;
                    end
                catch
                end
            end
        end
        if ycol
            id= 1;
            for i=2:size(info,1)
                if isnumeric(info{i,ycol })
                    score(id,icov) = info{i,ycol };
                else
                    score(id,icov) = nan;
                end
                id = id+1;
            end

        else
            disp(['Look for:' covariableall])
            disp([Pearsony,' column not found, no regression could be compute'])
            out='Stat tests';
            return
        end
    end




    mat_n_used = zeros(size(MATall,2));
    if strcmp(strtrim(upper(job.c_statmatrix.b_GLM_Mat.e_GLMGR)),'ALL');
        subjecttoapply = 1:numel(groupeall);
    else
        gruse  = str2num(job.c_statmatrix.b_GLM_Mat.e_GLMGR);
        subjecttoapply = find(sum(groupeall==gruse,2));
    end
    if isfield(job.c_statmatrix.b_GLM_Mat.c_statpermutation,'b_Nopermutation')
        for icov = 1:numel(covariableall)
            eval(['bCOV',num2str(icov),' = zeros(size(MATall,2),size(MATall,2));']);
            eval(['bCOV',num2str(icov),'sig = zeros(size(MATall,2),size(MATall,2));']);
        end
        for i=1:size(MATall,2)
            for j=1:size(MATall,2)
                iduse = find(sum(~isnan(score(subjecttoapply,:)),2)==size(score,2)& ~isnan(MATall(subjecttoapply,i,j)));
                mat_n_used(i,j) = numel(iduse);
                X = score(subjecttoapply(iduse),:);
                y = MATall(subjecttoapply(iduse),i,j);

                if ~isempty(iduse)
                    %R2 statistic, the F-statistic and its p-value, and an estimate of the error variance.
                    try
                        [b,bint,r,rint,stats] = regress(y,X);
                        [bp,bintp,rp,rintp,statsp] = regress(y,X(:,[1,3:end]));
                    catch
                        b(:) = nan;
                        stats(:) =  nan;

                    end
                    try
                        R2inc(i,j) = stats(1)-statsp(1);
                        %f2(i,j) = R2inc(i,j)./(1-R2inc(i,j)); %SPSS version Semi partial F
                        f2(i,j) = R2inc(i,j)./(1-stats(1)); %COHEN partial F
                    catch
                    end
                    try
                        for icov = 1:numel(covariableall)
                            eval(['bCOV',num2str(icov),'(',num2str(i),',',num2str(j),')=',num2str(b(icov)),';']);
                        end
                        eval(['r2','(',num2str(i),',',num2str(j),')=',num2str(stats(1)),';']);
                        eval(['F','(',num2str(i),',',num2str(j),')=',num2str(stats(2)),';']);
                        eval(['pval','(',num2str(i),',',num2str(j),')=',num2str(stats(3)),';']);
                        eval(['RMSerr','(',num2str(i),',',num2str(j),')=',num2str(sqrt(stats(4))),';']);
                    catch
                        1;
                    end
                    if stats(3)<alpha_threshold
                        for icov = 1:numel(covariableall)
                            eval(['bCOV',num2str(icov),'sig(',num2str(i),',',num2str(j),')=',num2str(b(icov)),';']);
                        end
                    end


                end
            end
        end


        1
        %check Frand
        if 1 %save residual covariable substracted
            infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}, covariableall];
            idcov2substract = 0;
            for id=1:numel(covariableall);
                if strcmp(covariableall{id},strtrim(job.c_statmatrix.b_GLM_Mat.b_substractidCovariable_Mat ));
                    idcov2substract = id; %indice de la covariable a soustraire, dans l'ordre d'entr�e
                end
            end
            if idcov2substract;
                residual = nan(size(MATall));
                residualsig = nan(size(MATall));
                for i=1:size(MATall,2);
                    for j=1:size(MATall,2);
                        iduse = find(sum(~isnan(score),2)==size(score,2)& ~isnan(MATall(:,i,j)));
                        X = score(iduse,:);
                        residual(iduse,i,j) = MATall(iduse,i,j) - eval(['bCOV',num2str(idcov2substract),'(i,j)'])*X(:,idcov2substract);
                        residualsig(iduse,i,j) = MATall(iduse,i,j) - eval(['bCOV',num2str(idcov2substract),'sig(i,j)'])*X(:,idcov2substract);

                    end
                end
                disp([upper(job.c_statmatrix.b_GLM_Mat.e_GLMGR), ' group where used for the regression n=',num2str(numel(subjecttoapply))])
                % for isubject = 1:numel(subjecttoapply)
                %     meancorr = squeeze(residual(subjecttoapply(isubject),:,:));
                %     file = ['res ' covariableall{idcov2substract}, info{subjecttoapply(isubject)+1,2}];
                %     save(fullfile(info{subjecttoapply(isubject)+1,1},file),'ZoneList','meancorr');
                %     new = [info{subjecttoapply(isubject)+1,1},{file}, {ZONEid},{1}, num2cell(score(subjecttoapply(isubject),:))  ];
                %     infonew = [infonew;new];
                %     disp( ['Save residual: ', fullfile(info{subjecttoapply(isubject)+1,1},file)])
                % end


                 % end
                for isubject=1:numel(subjecttoapply)
                    if job.m_nodeunit==1 %channel
                       % file = ['G1', info{g1(ifile)+1,2},'c.mat'];
                        file = ['res',covariableall{idcov2substract}, info{subjecttoapply(isubject)+1,2},'c.mat'];
                    elseif job.m_nodeunit==2%zone
                       % file = ['G1', info{g1(ifile)+1,2},'z.mat'];
                         file = ['res',covariableall{idcov2substract}, info{subjecttoapply(isubject)+1,2},'z.mat'];
                    end
                    meancorr = squeeze(residual(subjecttoapply(isubject),:,:));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                     new = [dir1,{file}, {ZONEid},{1}, num2cell(score(subjecttoapply(isubject),:))];
                    infonew = [infonew;new];
                end


            else
                infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}, covariableall];
                % for isubject = 1:numel(subjecttoapply)
                %     file = [info{subjecttoapply(isubject)+1,2}];
                %     new = [info{subjecttoapply(isubject)+1,1},{file}, {ZONEid},{1}, num2cell(score(subjecttoapply(isubject),:)) ];
                %     infonew = [infonew;new];
                % end
                for isubject=1:numel(subjecttoapply)
                    if job.m_nodeunit==1 %channel
                       % file = ['G1', info{g1(ifile)+1,2},'c.mat'];
                        file = [info{subjecttoapply(isubject)+1,2},'c.mat'];
                    elseif job.m_nodeunit==2%zone
                       % file = ['G1', info{g1(ifile)+1,2},'z.mat'];
                         file = [info{subjecttoapply(isubject)+1,2},'z.mat'];
                    end
                    meancorr =  squeeze(MATall(isubject,:,:));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                     new = [dir1,{file}, {ZONEid},{1}, num2cell(score(subjecttoapply(isubject),:))];
                    infonew = [infonew;new];
                end

                disp(['No subtracting residual ']);
            end
        end




        %WRITE IN A NEW FILE
        for icov = 1:numel(covariableall)
            file = ['Beta',name,'_',covariableall{icov},'.mat'];
            eval(['meancorr =','bCOV',num2str(icov),';']);
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            new = [{dir1},{file}, {ZONEid},{0}, num2cell(zeros(1,numel(covariableall)))];
            infonew = [infonew;new];
            file = ['Beta',name,'_',covariableall{icov},'p',num2str(alpha_threshold),'.mat'];
            eval(['meancorr =','bCOV',num2str(icov),'sig;']);
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            new = [{dir1},{file}, {ZONEid},{0}, num2cell(zeros(1,numel(covariableall)))];
            infonew = [infonew;new];
        end

        try
            file = ['R2_cov=',covariablestring];
            meancorr =r2;
            save(fullfile(dir1,file),'ZoneList','meancorr');
            new = [{dir1},{file},{ZONEid},{0},num2cell(zeros(1,numel(covariableall)))];
            infonew = [infonew;new];
            disp( ['Save R2: ', fullfile(info{isubject,1},file)])
        catch
        end

        file = ['F stat_cov=',covariablestring];
        meancorr =F;
        save(fullfile(dir1,file),'ZoneList','meancorr');
        new = [{dir1},{file},{ZONEid},{0},num2cell(zeros(1,numel(covariableall)))];
        infonew = [infonew;new];
        disp( ['Save Fstat: ', fullfile(info{isubject,1},file)])

        file = ['1-pval_cov=',covariablestring];
        meancorr =1-pval
        save(fullfile(dir1,file),'ZoneList','meancorr');
        new = [{dir1},{file},{ZONEid},{0},num2cell(zeros(1,numel(covariableall)))];
        infonew = [infonew;new];
        disp( ['Save 1-pval: ', fullfile(info{isubject,1},file)])

        file = ['RMSerr_cov=',covariablestring];
        meancorr =RMSerr;
        save(fullfile(dir1,file),'ZoneList','meancorr');
        new = [{dir1},{file}, {ZONEid},{0} ,num2cell(zeros(1,numel(covariableall)))];
        infonew = [infonew;new];
        disp( ['Save Root Mean Square Error: ',fullfile(info{isubject,1},file)])

        file = ['n',name];
        meancorr = mat_n_used;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        new = [{dir1},{file}, {ZONEid},{0},num2cell(zeros(1,numel(covariableall)))];
        infonew = [infonew;new];

        file = ['f2',name];
        meancorr = f2;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        new = [{dir1},{file}, {ZONEid},{0},num2cell(zeros(1,numel(covariableall)))];
        infonew = [infonew;new];


        file = ['f2BetaPosive',name];
        meancorr = f2.*double(bCOV2>0);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        new = [{dir1},{file}, {ZONEid},{0},num2cell(zeros(1,numel(covariableall)))];
        infonew = [infonew;new];


        file = ['f2BetaNegative',name];
        meancorr = f2.*double(bCOV2<0);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        new = [{dir1},{file}, {ZONEid},{0},num2cell(zeros(1,numel(covariableall)))];
        infonew = [infonew;new];


        file = ['R2inc',name];
        meancorr = R2inc;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        new = [{dir1},{file}, {ZONEid},{0},num2cell(zeros(1,numel(covariableall)))];
        infonew = [infonew;new];
        try
            copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid));
        catch
        end
        if ismac
            % Code to run on Mac platform problem with xlswrite
            [filepath,name,ext] = fileparts(xlslistfile);
            writetxtfile(fullfile(dir1,['GLM.txt']),infonew);
            disp(['Result .txt file saved: ' fullfile(dir1,['GLM.txt'])]);
        else
            [filepath,name,ext] = fileparts(xlslistfile);
            try
                xlswrite(fullfile(dir1,['GLM.xlsx']),infonew);
                disp(['Result .xlsx file saved ' fullfile(dir1,['GLM.xlsx'])]);
            catch
                writetxtfile(fullfile(dir1,['GLM.txt']),infonew);
                disp(['Result .txt file saved: ' fullfile(dir1,['GLM.txt'])]);
            end
        end
    elseif isfield(job.c_statmatrix.b_GLM_Mat.c_statpermutation,'b_permutation')
        'WARNING NO TESTED !'

        halfMAT = MATall(:,idhalf);
        for icov = 1:numel(covariableall)
            eval(['bCOV',num2str(icov),' = zeros(1,size(halfMAT,2));']);
            eval(['bCOV',num2str(icov),'sig = zeros(1,size(halfMAT,2));']);
        end
        for ilink = 1:size(halfMAT,2) % ici
            fprintf('%s ',num2str(ilink));
            if mod(ilink,30)==0
                fprintf('\n ');
            end
            iduse = find(sum(~isnan(score(subjecttoapply,:)),2)==size(score,2)& ~isnan(halfMAT(subjecttoapply,ilink)));

            X = score(subjecttoapply(iduse),:);
            y =  halfMAT(subjecttoapply(iduse),ilink);

            if ~isempty(iduse)
                %R2 statistic, the F-statistic and its p-value, and an estimate of the error variance.
                try
                    [b,bint,r,rint,stats] = regress(y,X);
                    [bp,bintp,rp,rintp,statsp] = regress(y,X(:,[1,3:end]));
                catch
                    b(:) = nan;
                    stats(:) =  nan;

                end
                try
                    R2inc(ilink) = stats(1)-statsp(1);
                    %f2(i,j) = R2inc(i,j)./(1-R2inc(i,j)); %SPSS version Semi partial F
                    f2(ilink) = R2inc(ilink)./(1-stats(1)); %COHEN partial F
                catch
                end
                try
                    for icov = 1:numel(covariableall)
                        eval(['bCOV',num2str(icov),'(',num2str(ilink),')=',num2str(b(icov)),';']);
                    end
                    eval(['r2','(',num2str(ilink),')=',num2str(stats(1)),';']);
                    eval(['F','(',num2str(ilink),')=',num2str(stats(2)),';']);
                    eval(['pval','(',num2str(ilink),')=',num2str(stats(3)),';']);
                    eval(['RMSerr','(',num2str(ilink),')=',num2str(sqrt(stats(4))),';']);
                catch
                    1;
                end
                if stats(3)<alpha_threshold
                    for icov = 1:numel(covariableall)
                        eval(['bCOV',num2str(icov),'sig(',num2str(ilink),')=',num2str(b(icov)),';']);
                    end
                end
            end

        end




        try
            nperm = str2num(job.c_statmatrix.b_GLM_Mat.c_statpermutation.b_permutation.e_npermutation);
            load(fullfile(dir1,'FPermutation.mat'))
            disp(['Previous computed permutation ', fullfile(dir1,'FPermutation.mat '), num2str(size(F_rand,1)), ' link x ',  num2str(size(F_rand,2)),' permutations',...
                'were open successfully if you wish to compute it again delete the file'])
        catch
            tic
            for ilink =1:size(halfMAT,2) %PERMUTATION F TEST
                for iperm=1:nperm
                    fprintf('%s ',num2str(ilink));
                    if mod(ilink,30)==0
                        fprintf('\n ');
                    end
                    iduse = find(sum(~isnan(score(subjecttoapply,:)),2)==size(score,2)& ~isnan(halfMAT(subjecttoapply,ilink)));
                    permfix =  randperm(numel(iduse));
                    X = score(subjecttoapply(iduse),:);
                    y =  halfMAT(subjecttoapply(iduse(permfix)),ilink);

                    if ~isempty(iduse)
                        %R2 statistic, the F-statistic and its p-value, and an estimate of the error variance.
                        try
                            [b,bint,r,rint,stats] = regress(y,X);
                            [bp,bintp,rp,rintp,statsp] = regress(y,X(:,[1,3:end]));
                        catch
                            b(:) = nan;
                            stats(:) =  nan;
                        end
                        try
                            eval(['F_rand','(',num2str(ilink),',',num2str(iperm),')=',num2str(stats(2)),';']);
                        catch

                        end
                    end
                end
            end
            toc
            save(fullfile(dir1,'FPermutation.mat'), 'F_rand')
        end

        [fecdf,xecdf] = ecdf(F_rand(:));
        clustercritical=xecdf(sum(fecdf<(1-alpha_threshold)));
        disp(['critical for clustering p=', num2str(alpha_threshold),' F=',num2str(clustercritical)])
        if isfield(job.c_statmatrix.b_GLM_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting,'b_MCT_ClusterBased')

            statobs = F';

            if nperm <= size(F_rand,2)
                statrand =F_rand(:,1:nperm);
            else
                nbrepeted = ceil( nperm/size(F_rand,2))
                F_rand =  repmat(F_rand,1,nbrepeted );
                statrand =F_rand(:,1:nperm);
                tmp = randperm(numel(statrand));
                statrand(tmp) = statrand;
            end
            if 0
                fileMean = 'C:\NIRx-actiCHAMP_Analysed\ELAN\Resting2021\Matrices\Stat\0mOneSample_ttest\_OHBM_MatriceCOH_Elan_BAYLEY_REVIEWmai2023cGR1 OneSampleTtest mean.mat';
                load(fileMean);
                MeanAmp = meancorr;
                statobs = F' .*double(MeanAmp(idhalf)>0);
                disp('Mask positive Pearson before cluster')
                statrand =statrand(:,:).*(double(MeanAmp(idhalf)>0)*ones(1,nperm));
            end
            neighbourdist = job.c_statmatrix.b_GLM_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_neighbourdist;
            minnbchan = job.c_statmatrix.b_GLM_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_minnbchan;
            disp(['critical for clustering p=', num2str(alpha_threshold),' T=',num2str(clustercritical)])
            [stat, matneig] = FindClusterBasedPermutationInMatrix(chanpos, neighbourdist,clustercritical, statobs, statrand,  minnbchan );
            save(fullfile(dir1,'GLM_stat.mat'),'stat')
            disp(['Save cluster', fullfile(dir1,'GLM_stat.mat') ])
            infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
            %                  figure
            % hist(stat.posdistribution)
            if isfield(stat,'posclusterslabelmat')
                for iposcluster = 1:numel(stat.posclusters)
                    file = [name,labelnode,num2str(nperm),' posCluster', num2str(iposcluster),' MCp=',num2str(stat.posclusters(iposcluster).prob),'.mat'];
                    tmp=(stat.posclusterslabelmat==iposcluster);
                    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                    new = [{dir1},{file}, {ZONEid},{0} ];
                    infonew = [infonew;new];


                    %WRITE beta and cluster in file
                    for icov = 1:numel(covariableall)
                        file = ['Beta','_',covariableall{icov},'posCluster', num2str(iposcluster),' MCp=',num2str(stat.posclusters(iposcluster).prob),'.mat'];
                        meancorr = zeros(size(MATall,2),size(MATall,2));
                        eval(['tmp =','bCOV',num2str(icov),';']);
                        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                        disp(['Save file: ', fullfile(dir1,[file])]);
                        save(fullfile(dir1,[file]),'ZoneList','meancorr');
                        new = [{dir1},{file}, {ZONEid},{0}];
                        infonew = [infonew;new];
                    end
                end

            else
                disp('No positive cluster found')
            end
            if isfield(stat,'negclusterslabelmat') %technicaly no neg cluster in f value...
                for inegcluster = 1:numel(stat.negclusters)
                    file = [name,labelnode,num2str(nperm),' negCluster', num2str(inegcluster),' MCp=',num2str(stat.negclusters(inegcluster).prob),'.mat'];
                    tmp=(stat.negclusterslabelmat==inegcluster);
                    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                    new = [{dir1},{file}, {ZONEid},{0} ];
                    infonew = [infonew;new];

                    %WRITE beta and cluster in file
                    for icov = 1:numel(covariableall)
                        file = ['Beta','_',covariableall{icov},'negCluster', num2str(inegcluster),' MCp=',num2str(stat.negclusters(inegcluster).prob),'.mat']
                        meancorr = zeros(size(MATall,2),size(MATall,2));
                        eval(['meancorr(idhalf) =','bCOV',num2str(icov),';']);
                        tmp=(stat.negclusterslabelmat==inegcluster);
                        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                        disp(['Save file: ', fullfile(dir1,[file])]);
                        save(fullfile(dir1,[file]),'ZoneList','meancorr');
                        new = [{dir1},{file}, {ZONEid},{0}];
                        infonew = [infonew;new];
                    end

                end
            else
                disp('No negative cluster found')
            end

        end
        [filepath,name,ext] = fileparts(xlslistfile);
        try
            xlswrite(fullfile(dir1,['CLUSTER.xlsx']),infonew);
            disp(['Result .xlsx file saved ' fullfile(dir1,['CLUSTER.xlsx'])]);
        catch
            writetxtfile(fullfile(dir1,['CLUSTER.xlsx']),infonew);
            disp(['Result .xlsx file saved ' fullfile(dir1,['CLUSTER.xlsx'])]);
        end

    end
elseif isfield(job.c_statmatrix,'b_LM_Mat')
    disp('Not tested GENERIC USE LM')

    dir1 = job.e_statmatrixPath{1};
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    datatable= readtable(xlslistfile)
    halfMAT = MATall(:,idhalf);
    LMformula = deblank(job.c_statmatrix.b_LM_Mat.b_LME_formula);
    fprintf('%s',['RUN LM : ', LMformula, ' link:'])
    try %load precomputed to speed up not keep as standard way of working
        datatable.MAT= halfMAT(:,1);
        lm = fitlm( datatable,LMformula);
        for icov = 1:size(lm.Coefficients,1)
            file = ['TCOV',num2str(icov),'_', strtrim(lm.CoefficientNames{1,icov}),'.mat'];
            filename = fullfile(dir1,file);
            eval(['load(''', filename,''')']);
            halfMATtemp = meancorr(idhalf);
            eval(['TCOV',num2str(icov),' = halfMATtemp']);
        end
        disp('load precomputed tval')
    catch

        for ilink=1:numel(idhalf);
            fprintf('%d,',ilink)
            if mod(ilink,30)==0
                fprintf('\n ');
            end
            datatable.MAT= halfMAT(:,ilink);
            %  LMEformula = 'AGENIRS_m~MAT+(1|ID)'
            try
                lm = fitlm( datatable,LMformula);
            catch
                disp(['Verify LM formula ', LMformula, ' and datatable Variable'])
                disp(datatable.Properties.VariableNames)
            end
            eval(['r2','(',num2str(ilink),')=' ,'lm.Rsquared.Ordinary;'])
            eval(['N','(',num2str(ilink),')=' ,' lm.NumObservations; '])

            for icov = 1:size(lm.Coefficients,1)
                eval(['TCOV',num2str(icov),'(',num2str(ilink),')=lm.Coefficients{',num2str(icov),',3};']);
                eval(['ECOV',num2str(icov),'(',num2str(ilink),')=lm.Coefficients{',num2str(icov),',1};']);
                eval(['pCOV',num2str(icov),'(',num2str(ilink),')=lm.Coefficients{',num2str(icov),',4};']);
            end

        end

        file = ['Nobservation.mat']
        meancorr = zeros(size(MATall,2),size(MATall,2));
        eval(['tmp =','N;']);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        disp(['Save file: ', fullfile(dir1,[file])]);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        new = [{dir1},{file}, {ZONEid},{0}];
        infonew = [infonew;new];

        file = ['R2.mat']
        meancorr = zeros(size(MATall,2),size(MATall,2));
        eval(['tmp =','r2']);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        disp(['Save file: ', fullfile(dir1,[file])]);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        new = [{dir1},{file}, {ZONEid},{0}];
        infonew = [infonew;new];

        % WRITE Tval in file
        for icov = 1:size(lm.Coefficients,1)
            file = ['TCOV',num2str(icov),'_', lm.CoefficientNames{1,icov},'.mat']
            eval(['tmp =','TCOV',num2str(icov),';']);
            eval(['matsig =','pCOV',num2str(icov),';']);

            tmp = tmp.* double(matsig<alpha_threshold);
            meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
            disp(['Save file: ', fullfile(dir1,[file])]);
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            new = [{dir1},{file}, {ZONEid},{0}];
            infonew = [infonew;new];
        end

        for icov = 1:size(lm.Coefficients,1)
            file = ['ECOV',num2str(icov),'_', lm.CoefficientNames{1,icov},'.mat']
            meancorr = zeros(size(MATall,2),size(MATall,2));
            eval(['tmp =','ECOV',num2str(icov),';']);
            eval(['matsig =','pCOV',num2str(icov),';']);
            tmp = tmp.* double(matsig<alpha_threshold);
            meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
            disp(['Save file: ', fullfile(dir1,[file])]);
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            new = [{dir1},{file}, {ZONEid},{0}];
            infonew = [infonew;new];
        end
        try
            xlswrite(fullfile(dir1,['lm', LMformula, '.xlsx']),infonew);
            disp(['Result .xlsx file saved ' fullfile(dir1,['lm', LMformula, '.xlsx'])]);
        catch
            writetxtfile(fullfile(dir1,['lm', LMformula, '.xlsx']),infonew)
            disp(['Result .xlsx file saved ' fullfile(dir1,['lm', LMformula, '.xlsx'])]);
        end
    end
    if isfield(job.c_statmatrix.b_LM_Mat.c_statpermutation,'b_permutation')
        try
            nperm = str2num(job.c_statmatrix.b_LM_Mat.c_statpermutation.b_permutation.e_npermutation);
            for icov = 1:size(lm.Coefficients,1)
                eval(['load(''',dir1,'\TrandCOV',num2str(icov),'.mat'')']);
            end
            disp(['Previous computed permutation ',...
                'were open successfully => delete the file if you wish to compute it again '])
        catch
            tic
            for ilink =1:size(halfMAT,2) %PERMUTATION F TEST
                for iperm=1:nperm
                    fprintf('%s ',num2str(ilink));
                    if mod(ilink,30)==0
                        fprintf('\n ');
                    end

                    fprintf('%d,',ilink);
                    if mod(ilink,30)==0
                        fprintf('\n ');
                    end
                    permfix =  randperm(size(halfMAT,1));
                    datatable.MAT= halfMAT(permfix,ilink);
                    lme = fitlm( datatable,LMformula);
                    for icov = 1:size(lm.Coefficients,1)
                        eval(['TrandCOV',num2str(icov),'(',num2str(ilink),',',num2str(iperm),')=lm.Coefficients{',num2str(icov),',3};']);

                    end
                end
            end
            toc
            for icov = 1:size(lme.Coefficients,1)
                eval(['save(''',dir1,'\TrandCOV',num2str(icov),'.mat'',''TrandCOV',num2str(icov),''')']);
            end

        end

        infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
        if isfield(job.c_statmatrix.b_LM_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting,'b_MCT_ClusterBased')
            for icov = 1:size(lm.Coefficients,1)
                eval(['[fecdf,xecdf] = ecdf(TrandCOV',num2str(icov),'(:));']);
                clustercritical=xecdf(sum(fecdf<(1-alpha_threshold)));
                %clustercritical=xecdf(sum(fecdf<(1-0.01)))
                disp(['critical for clustering cov: ' , lm.CoefficientNames{1,icov},' p=', num2str(alpha_threshold),' T=',num2str(clustercritical)])
                try
                    nperm = str2num(job.c_statmatrix.b_LME_Mat.c_statpermutation.b_permutation.e_npermutation);
                    for icov = 1:size(lme.Coefficients,1)
                        eval(['load(''',dir1,'\TrandCOV',num2str(icov),'.mat'')']);
                    end
                    disp(['Previous computed permutation ',...
                        'were open successfully => delete the file if you wish to compute it again '])
                catch
                end
                %statobs = TCOV2';
                eval(['statobs = TCOV',num2str(icov),';']);
                eval(['TrandCOV = TrandCOV',num2str(icov),';']);
                if nperm <= size(TrandCOV,2)
                    statrand =TrandCOV(:,1:nperm);
                else
                    nbrepeted = ceil( nperm/size(TrandCOV,2))
                    F_rand =  repmat(TrandCOV,1,nbrepeted );
                    statrand =TrandCOV(:,1:nperm);
                    tmp = randperm(numel(statrand));
                    statrand(tmp) = statrand;
                end
                neighbourdist = job.c_statmatrix.b_LM_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_neighbourdist;
                minnbchan = job.c_statmatrix.b_LM_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_minnbchan;
                disp(['critical for clustering p=', num2str(alpha_threshold),' T=',num2str(clustercritical)])
                disp(['Neighbor distance ', neighbourdist])
                disp(['Minimal commun neighbor ', minnbchan])

                [stat, matneig] = FindClusterBasedPermutationInMatrix(chanpos,neighbourdist,clustercritical, statobs, statrand,  minnbchan );
                %                     [statsig, matneig] = FindClusterBasedPermutationInMatrix(chanpos, 'sig',clustercritical, statobs, -statrand,  minnbchan );
                %
                %          figure;hist(statrand(:))


                save(fullfile(dir1,'FiedTrip_LMCluster_stat.mat'),'stat')
                disp(['Save fieldtrip cluster stat ', fullfile(dir1,'FiedTrip_LMCluster_stat',neighbourdist,'_',minnbchan,'.mat') ])


                if isfield(stat,'posclusterslabelmat')
                    for iposcluster = 1:numel(stat.posclusters)
                        file = [name,labelnode,num2str(nperm),'COV', num2str(icov),' posCluster', num2str(iposcluster),' MCp=',num2str(stat.posclusters(iposcluster).prob),'.mat'];
                        tmp=(stat.posclusterslabelmat==iposcluster);
                        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                        save(fullfile(dir1,[file]),'ZoneList','meancorr');
                        disp(['Save: ', fullfile(dir1,[file])]);
                        new = [{dir1},{file}, {ZONEid},{0} ];
                        infonew = [infonew;new];

                        %WRITE beta and cluster in file
                        file = ['Tval','COV', num2str(icov),'_','posCluster', num2str(iposcluster),' MCp=',num2str(stat.posclusters(iposcluster).prob),'.mat'];
                        meancorr = zeros(size(MATall,2),size(MATall,2));
                        eval(['tmp =','TCOV',num2str(icov),';']);
                        matsig=(stat.posclusterslabelmat==iposcluster);
                        tmp = tmp.*matsig;
                        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                        disp(['Save file: ', fullfile(dir1,[file])]);
                        save(fullfile(dir1,[file]),'ZoneList','meancorr');
                        new = [{dir1},{file}, {ZONEid},{0}];
                        infonew = [infonew;new];
                    end

                else
                    disp('No positive cluster found')
                end
                if isfield(stat,'negclusterslabelmat') %technicaly no neg cluster in f value...
                    for inegcluster = 1:numel(stat.negclusters)
                        file = [name,labelnode,num2str(nperm),'COV', num2str(icov),' negCluster', num2str(inegcluster),' MCp=',num2str(stat.negclusters(inegcluster).prob),'.mat'];
                        tmp=(stat.negclusterslabelmat==inegcluster);
                        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                        save(fullfile(dir1,[file]),'ZoneList','meancorr');
                        disp(['Save: ', fullfile(dir1,[file])]);
                        new = [{dir1},{file}, {ZONEid},{0} ];
                        infonew = [infonew;new];

                        %WRITE beta and cluster in file
                        file = ['Tval','COV', num2str(icov),'negCluster', num2str(inegcluster),' MCp=',num2str(stat.negclusters(inegcluster).prob),'.mat']

                        eval(['tmp =','TCOV',num2str(icov),';']);
                        matsig=(stat.negclusterslabelmat==inegcluster);
                        tmp =  tmp.* double(matsig);
                        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                        disp(['Save file: ', fullfile(dir1,[file])]);
                        save(fullfile(dir1,[file]),'ZoneList','meancorr');
                        new = [{dir1},{file}, {ZONEid},{0}];
                        infonew = [infonew;new];


                    end
                else
                    disp('No negative cluster found')
                end
            end
        end
        [filepath,name,ext] = fileparts(xlslistfile);
        try
            xlswrite(fullfile(dir1,['CLUSTER',neighbourdist,'_',minnbchan,'_',num2str(alpha_threshold),' .xlsx']),infonew);
            disp(['Result .xlsx file saved ' fullfile(dir1,['CLUSTER',neighbourdist,'_',minnbchan,'_',num2str(alpha_threshold),'.xlsx'])]);
        catch
            writetxtfile(fullfile(dir1,['CLUSTER',neighbourdist,'_',minnbchan,'_',num2str(alpha_threshold),'.xlsx']),infonew);
            disp(['Result .xlsx file saved ' fullfile(dir1,['CLUSTER',neighbourdist,'_',minnbchan,'_',num2str(alpha_threshold),'.xlsx'])]);
        end
    end
    %


elseif isfield(job.c_statmatrix,'b_LME_Mat')
    disp('Warning function LME in developpement');
    dir1 = job.e_statmatrixPath{1};
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    datatable= readtable(xlslistfile);
    halfMAT = MATall(:,idhalf);
    LMEformula = job.c_statmatrix.b_LME_Mat.b_LME_formula;
    try %load precomputed to speed up not keep as standard way of working
        datatable.MAT= halfMAT(:,1);
        lme = fitlme( datatable,LMEformula);
        for icov = 1:size(lme.Coefficients,1)
            file = ['TCOV',num2str(icov),'_', strtrim(lme.Coefficients{icov,1}),'.mat'];
            filename = fullfile(dir1,file);
            eval(['load(''', filename,''')']);
            halfMATtemp = meancorr(idhalf);
            eval(['TCOV',num2str(icov),' = halfMATtemp'])
        end
        disp('load precomputed tval')
    catch



        fprintf('%s',['RUN LME : ', LMEformula, ' link:'])
        for ilink=1:numel(idhalf);
            fprintf('%d,',ilink)
            if mod(ilink,30)==0
                fprintf('\n ');
            end
            datatable.MAT= halfMAT(:,ilink);  
            try
                lme = fitlme( datatable,LMEformula);
                eval(['R2ord(',num2str(ilink),')=lme.Rsquared.Ordinary;']);
                eval(['R2adj(',num2str(ilink),')=lme.Rsquared.Adjusted;']);
                for icov = 1:size(lme.Coefficients,1)
                    eval(['TCOV',num2str(icov),'(',num2str(ilink),')=lme.Coefficients{',num2str(icov),',4};']);
                    eval(['ECOV',num2str(icov),'(',num2str(ilink),')=lme.Coefficients{',num2str(icov),',2};']);
                    eval(['pCOV',num2str(icov),'(',num2str(ilink),')=lme.Coefficients{',num2str(icov),',6};']);
                end
                [beta,betanames,stats] =  fixedEffects(lme);
                [B,Bnames,statsrandom] = randomEffects(lme);
            % Group	Grouping variable associated with the random effect
            % Level	Level within the grouping variable corresponding to the random effect
            % Name	Name of the random-effect coefficient
            % Estimate	Best linear unbiased predictor (BLUP) of random effect
            % SEPred	Standard error of the estimate (BLUP minus random effect)
            % tStat	t-statistic for a test that the random effect is zero
            % DF	Estimated degrees of freedom for the t-statistic
            % pValue	p-value for the t-statistic
            % Lower	Lower limit of a 95% confidence interval for the random effect
            % Upper
            for irand = 1:size(statsrandom,1)
                eval(['TRand(',num2str(irand),',',num2str(ilink),')=statsrandom{',num2str(irand),',6};']);
                eval(['ERand(',num2str(irand),',',num2str(ilink),')=statsrandom{',num2str(irand),',4};']);
                eval(['pRand(',num2str(irand),',',num2str(ilink),')=statsrandom{',num2str(irand),',8};']);

            end

            catch
                try %v�rifier si l'erreur est du au valeurs manquante de la diagonale ou la formule 
                    datatable.MAT = ones(size(halfMAT,1),1);
                    lme = fitlme( datatable,LMEformula);  
                    [B,Bnames,statsrandom] = randomEffects(lme);
                catch
                      disp(['Verify LME formula ', LMEformula, ' and datatable Variable'])
                      disp(datatable.Properties.VariableNames) 
                end
    
                eval(['R2ord(',num2str(ilink),')=nan;']);
                eval(['R2adj(',num2str(ilink),')=nan;']);
                for icov = 1:size(lme.Coefficients,1)
                    eval(['TCOV',num2str(icov),'(',num2str(ilink),')=nan;']);
                    eval(['ECOV',num2str(icov),'(',num2str(ilink),')=nan;']);
                    eval(['pCOV',num2str(icov),'(',num2str(ilink),')=nan;']);
                end              
                 for irand = 1:size(statsrandom,1)
                eval(['TRand(',num2str(irand),',',num2str(ilink),')=nan;']);
                eval(['ERand(',num2str(irand),',',num2str(ilink),')=nan;']);
                eval(['pRand(',num2str(irand),',',num2str(ilink),')=nan;']);
                 end      
            end           
        end
   
        % WRITE result in file

        file = ['R2Ordinary.mat'];
        meancorr = zeros(size(MATall,2),size(MATall,2));
        eval(['tmp =','R2ord;']);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        disp(['Save file: ', fullfile(dir1,[file])]);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        new = [{dir1},{file}, {ZONEid},{0}];
        infonew = [infonew;new];

        file = ['R2Ajusted.mat'];

        eval(['tmp =','R2adj;']);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        disp(['Save file: ', fullfile(dir1,[file])]);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        new = [{dir1},{file}, {ZONEid},{0}];
        infonew = [infonew;new];

        for icov = 1:size(lme.Coefficients,1)
            file = ['TCOV',num2str(icov),'_', lme.Coefficients{icov,1},'.mat'];
            meancorr = zeros(size(MATall,2),size(MATall,2));
            eval(['tmp =','TCOV',num2str(icov),';']);
            eval(['matsig =','pCOV',num2str(icov),';']);
            tmp=  tmp.* double(matsig<alpha_threshold);
            meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
            disp(['Save file: ', fullfile(dir1,[file])]);
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            new = [{dir1},{file}, {ZONEid},{0}];
            infonew = [infonew;new];
        end



        for icov = 1:size(lme.Coefficients,1)
            file = ['ECOV',num2str(icov),'_', lme.Coefficients{icov,1},'.mat'];

            eval(['tmp =','ECOV',num2str(icov),';']);
            eval(['matsig =','pCOV',num2str(icov),';']);
            tmp =  tmp.* double(matsig<alpha_threshold);
            meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
            disp(['Save file: ', fullfile(dir1,[file])]);
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            new = [{dir1},{file}, {ZONEid},{0}];
            infonew = [infonew;new];
        end


        %  for irand = 1:size(statsrandom,1)
        %       %save random effect Tval
        %     file = ['Trand ', statsrandom{irand,1}, statsrandom{irand,2},statsrandom{irand,3},'.mat'];
        %     meancorr = zeros(size(MATall,2),size(MATall,2));
        %     eval(['meancorr(idhalf) =','TRand(',num2str(irand),',:);']);
        %     meancorr = meancorr +flipud(rot90(meancorr ));
        %     matsig = zeros(size(MATall,2),size(MATall,2));
        %     eval(['matsig =','pRand(',num2str(irand),',:);']);
        %     meancorr =  meancorr.* double(matsig<alpha_threshold);
        %     disp(['Save file: ', fullfile(dir1,[file])]);
        %     save(fullfile(dir1,[file]),'ZoneList','meancorr');
        %     new = [{dir1},{file}, {ZONEid},{0}];
        %     infonew = [infonew;new];
        %
        %
        %  %save random effect results Estimate
        %
        %     file = ['Erand ', statsrandom{irand,1}, statsrandom{irand,2},statsrandom{irand,3},'.mat'];
        %     meancorr = zeros(size(MATall,2),size(MATall,2));
        %     eval(['meancorr(idhalf) =','ERand(',num2str(irand),',:);']);
        %     meancorr = meancorr +flipud(rot90(meancorr ));
        %     matsig = zeros(size(MATall,2),size(MATall,2));
        %     eval(['matsig(idhalf) =','pRand(',num2str(irand),',:);']);
        %     matsig = matsig +flipud(rot90(matsig));
        %     meancorr =  meancorr.* double(matsig<alpha_threshold);
        %     disp(['Save file: ', fullfile(dir1,[file])]);
        %     save(fullfile(dir1,[file]),'ZoneList','meancorr');
        %     new = [{dir1},{file}, {ZONEid},{0}];
        %     infonew = [infonew;new];
        % end


        try

            xlswrite(fullfile(dir1,['lme', 'model', '.xlsx']),infonew);
            disp(['Result .xlsx file saved ' fullfile(dir1,['lme', LMEformula, '.xlsx'])]);
        catch
            writetxtfile(fullfile(dir1,['lme', 'model', '.xlsx']),infonew)
            disp(['Result .xlsx file saved ' fullfile(dir1,['lme', LMEformula, '.xlsx'])]);
        end
    end %use precomputed value to speed up test
    if isfield(job.c_statmatrix.b_LME_Mat.c_statpermutation,'b_permutation')
        try
            nperm = str2num(job.c_statmatrix.b_LME_Mat.c_statpermutation.b_permutation.e_npermutation);
            for icov = 1:size(lme.Coefficients,1)
                eval(['load(''',dir1,'\TrandCOV',num2str(icov),'.mat'')']);
            end
            disp(['Previous computed permutation ',...
                'were open successfully => delete the file if you wish to compute it again '])
        catch
            tic
            for ilink =933:1500 %1:size(halfMAT,2) %PERMUTATION TEST 216
                for iperm=1:nperm
                    fprintf('%s ',num2str(ilink));
                    if mod(ilink,30)==0
                        fprintf('\n ');
                    end

                    fprintf('%d,',ilink);
                    if mod(ilink,30)==0
                        fprintf('\n ');
                    end
                    permfix =  randperm(size(halfMAT,1));
                    datatable.MAT= halfMAT(permfix,ilink);
                    lme = fitlme( datatable,LMEformula);
                    for icov = 1:size(lme.Coefficients,1)
                        eval(['TrandCOV',num2str(icov),'(',num2str(ilink),',',num2str(iperm),')=lme.Coefficients{',num2str(icov),',4};']);

                    end
                end
            end
            toc
            for icov = 1:size(lme.Coefficients,1)
                eval(['save(''',dir1,'\TrandCOV',num2str(icov),'.mat'',''TrandCOV',num2str(icov),''')']);
            end

        end
        infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
        for icov = 1:size(lme.Coefficients,1)
            eval(['[fecdf,xecdf] = ecdf(TrandCOV', num2str(icov),'(:));']);
            clustercritical=xecdf(sum(fecdf<(1-alpha_threshold)));
            disp(['critical for clustering p=', num2str(alpha_threshold),' T=',num2str(clustercritical)])
            if isfield(job.c_statmatrix.b_LME_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting,'b_MCT_ClusterBased')

                eval(['statobs = TCOV', num2str(icov),';']);
                if nperm <= size(TrandCOV2,2)
                    eval(['statrand = TrandCOV', num2str(icov),'(:,1:nperm);']);
                else
                    eval(['nbrepeted = ceil( nperm/size(TrandCOV',num2str(icov),',2));']);
                    eval(['F_rand =  repmat(TrandCOV', num2str(icov),',1,nbrepeted );']);
                    eval(['statrand =TrandCOV', num2str(icov),'(:,1:nperm);']);
                    %                     nbrepeted = ceil( nperm/size(TrandCOV2,2))
                    %                     F_rand =  repmat(TrandCOV2,1,nbrepeted );
                    %                     statrand =TrandCOV2(:,1:nperm);
                    tmp = randperm(numel(statrand));
                    statrand(tmp) = statrand;
                end
                neighbourdist = job.c_statmatrix.b_LME_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_neighbourdist;
                minnbchan = job.c_statmatrix.b_LME_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_minnbchan;
                disp(['critical for clustering p=', num2str(alpha_threshold),' T=',num2str(clustercritical)])

                [stat, matneig] = FindClusterBasedPermutationInMatrix(chanpos, neighbourdist,clustercritical, statobs(:), statrand,  minnbchan );
                %  [stat, matneig] = FindClusterBasedPermutationInMatrix(chanpos, 'distance',clustercritical, statobs, statrand );
                save(fullfile(dir1,'FiedTrip_LMECluster_stat.mat'),'stat')
                disp(['Save fieldtrip cluster stat ', fullfile(dir1,'FiedTrip_LMECluster_stat.mat') ])

                %                  figure
                % hist(stat.posdistribution)
                if isfield(stat,'posclusterslabelmat')
                    for iposcluster = 1:numel(stat.posclusters)
                        file = [name,labelnode,num2str(nperm),'COV', num2str(icov), 'posCluster', num2str(iposcluster),' MCp=',num2str(stat.posclusters(iposcluster).prob),'.mat'];
                        tmp=(stat.posclusterslabelmat==iposcluster);
                        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                        save(fullfile(dir1,[file]),'ZoneList','meancorr');
                        disp(['Save: ', fullfile(dir1,[file])]);
                        new = [{dir1},{file}, {ZONEid},{0} ];
                        infonew = [infonew;new];

                        %WRITE beta and cluster in file
                        file = ['Tval','COV', num2str(icov),'_','posCluster', num2str(iposcluster),' MCp=',num2str(stat.posclusters(iposcluster).prob),'.mat'];
                        eval(['tmp =','TCOV',num2str(icov),';']);
                        matsig=(stat.posclusterslabelmat==iposcluster);

                        tmp =  tmp.* double(matsig);
                        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                        disp(['Save file: ', fullfile(dir1,[file])]);
                        save(fullfile(dir1,[file]),'ZoneList','meancorr');
                        new = [{dir1},{file}, {ZONEid},{0}];
                        infonew = [infonew;new];
                    end

                else
                    disp('No positive cluster found')
                end
                if isfield(stat,'negclusterslabelmat') %technicaly no neg cluster in f value...
                    for inegcluster = 1:numel(stat.negclusters)
                        file = [name,labelnode,num2str(nperm),'COV', num2str(icov),' negCluster', num2str(inegcluster),' MCp=',num2str(stat.negclusters(inegcluster).prob),'.mat'];
                        tmp=(stat.negclusterslabelmat==inegcluster);
                        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                        save(fullfile(dir1,[file]),'ZoneList','meancorr');
                        disp(['Save: ', fullfile(dir1,[file])]);
                        new = [{dir1},{file}, {ZONEid},{0} ];
                        infonew = [infonew;new];
                        %WRITE beta and cluster in file
                        file =  ['Tval','COV',num2str(icov),'negCluster', num2str(inegcluster),' MCp=',num2str(stat.negclusters(inegcluster).prob),'.mat']
                        meancorr = zeros(size(MATall,2),size(MATall,2));
                        eval(['tmp =','TCOV',num2str(icov),';']);

                        matsig=(stat.negclusterslabelmat==inegcluster);
                        tmp = tmp.* double(matsig);
                        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                        disp(['Save file: ', fullfile(dir1,[file])]);
                        save(fullfile(dir1,[file]),'ZoneList','meancorr');
                        new = [{dir1},{file}, {ZONEid},{0}];
                        infonew = [infonew;new];
                    end
                else
                    disp('No negative cluster found')
                end
            end
        end
        [filepath,name,ext] = fileparts(xlslistfile);
        try
            xlswrite(fullfile(dir1,['CLUSTER.xlsx']),infonew);
            disp(['Result .xlsx file saved ' fullfile(dir1,['CLUSTER.xlsx'])]);
        catch
            writetxtfile(fullfile(dir1,['CLUSTER.xlsx']),infonew);
            disp(['Result .xlsx file saved ' fullfile(dir1,['CLUSTER.xlsx'])]);
        end
    end
    %

elseif isfield(job.c_statmatrix,'b_PairedTtest')
    AllC = [];
    id =1;
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    %Use a specific groupev

    idG1 = find( sum( groupeall==job.c_statmatrix.b_PairedTtest.e_TtestOneSampleGR,2));
    idG2 = find(  sum(  groupeall==job.c_statmatrix.b_PairedTtest.e_TtestOneSampleGR2,2));
    GRname = [' GR',num2str(job.c_statmatrix.b_PairedTtest.e_TtestOneSampleGR),'vs',num2str(job.c_statmatrix.b_PairedTtest.e_TtestOneSampleGR2),' '];
    iduse = [idG1,idG2];

    halfMAT = MATall(:,idhalf);
    MATallG1 =  halfMAT( idG1,:);
    MATallG2 =  halfMAT( idG2,:);

    for ifile=1:numel(idG1)
        file = ['G1', info{idG1(ifile)+1,2}];
        meancorr = squeeze(MATall(idG1(ifile),:,:));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{1} ];
        infonew = [infonew;new];
    end
    for ifile=1:numel(idG2)
        file = ['G2', info{idG2(ifile)+1,2}];
        meancorr = squeeze(MATall(idG2(ifile),:,:));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{2} ];
        infonew = [infonew;new];
    end

    try
        matdiff = MATallG1 -  MATallG2;
    catch
        disp(['Error: Group ', num2str(job.c_statmatrix.b_PairedTtest.e_TtestOneSampleGR),...
            ' has ', num2str(size( MATallG1,1)), ' subjects that must be paired with Group ',...
            num2str(job.c_statmatrix.b_PairedTtest.e_TtestOneSampleGR2), ' which has ',...
            num2str(size( MATallG2,1)),' subjects'])
        out='Stat tests';
        return
    end
    %savein result folder
    dir1 = job.e_statmatrixPath{1};
    for ifile=1:size(matdiff,1)
        file = ['Paired',sprintf('%03.0f',(ifile)),' ',    list_subject{idG1(ifile) } '-', list_subject{idG2(ifile)}];
        tmp =squeeze(MATallG1(ifile,:) -  MATallG2(ifile,:));
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{3} ];
        infonew = [infonew;new];
    end
    disp(['Perform Groupe ', num2str(job.c_statmatrix.b_PairedTtest.e_TtestOneSampleGR) ' - Groupe ', num2str(job.c_statmatrix.b_PairedTtest.e_TtestOneSampleGR2)])
    meanall = squeeze(nanmean(  matdiff,1));
    refval  = job.c_statmatrix.b_PairedTtest.e_TtestOneSample_meanvalue;
    disp(['Paired ttest vs hypothetical mean ', num2str( job.c_statmatrix.b_PairedTtest.e_TtestOneSample_meanvalue)])

    tval = squeeze((nanmean( matdiff(:,:),1)-refval)./(nanstd( matdiff(:,:),1)./sqrt(sum(~isnan( matdiff(:,:)),1))  ));
    dfall = squeeze(sum(~isnan( matdiff(:,:)),1))-1;
    for i=1:size(tval,1)
        for j=1:size(tval,2)

            try
                % Compute the correct p-value for the test, and confidence intervals
                % if requested.
                if job.c_statmatrix.b_PairedTtest.m_TtestOneSample_matrix == 1; % two-tailed test
                    pval(i,j) = 2 * tcdf(-abs(-tval(i,j)), dfall(i,j));
                    %     if nargout > 2
                    %         crit = tinv((1 - alpha / 2), df) .* ser;
                    %         ci = cat(dim, xmean - crit, xmean + crit);
                    %     end
                elseif job.c_statmatrix.b_PairedTtest.m_TtestOneSample_matrix == 3; % right one-tailed test
                    pval(i,j) = tcdf(-tval(i,j), dfall(i,j));
                    %     if nargout > 2
                    %         crit = tinv(1 - alpha, df) .* ser;
                    %         ci = cat(dim, xmean - crit, Inf(size(p)));
                    %     end
                elseif job.c_statmatrixb_PairedTtest.m_TtestOneSample_matrix == 2; % left one-tailed test
                    pval(i,j) = tcdf(tval(i,j), dfall(i,j));
                    %     if nargout > 2
                    %         crit = tinv(1 - alpha, df) .* ser;
                    %         ci = cat(dim, -Inf(size(p)), xmean + crit);
                    %     end
                end
                if i==j
                    pval(i,j)=1;
                    tval(i,j)=0;
                end
            catch
                pval(i,j) = nan;
            end
        end
    end
    dir1 = job.e_statmatrixPath{1};

    if ~isdir(dir1)
        mkdir(dir1);
    end
    try
        [FDR,Q] = mafdr(pval(:));
        Q = reshape(Q,size(pval));
    catch
        Q = pval()
        Q(:) = 1
    end

    tmp = squeeze(meanall(:));
    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
    totaltrialgood = mean(dfall(:));
    file = [name,labelnode,GRname,'PAIREDTtest mean','.mat'];
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];

    file = [name,labelnode,GRname,'PAIREDTtest meanp',num2str(alpha_threshold),'.mat'];
    tmp = meanall.*double(pval<alpha_threshold);
    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];

    file = [name,labelnode,GRname,'PAIREDTtest meanFDR',num2str(alpha_threshold),'.mat'];
    tmp =  meanall.*double(Q<alpha_threshold);
    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];


    % ZoneList = MAT.ZoneList;
    totaltrialgood = mean(dfall(:));
    tmp = tval;
    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
    file = [name,labelnode,GRname,'PAIREDTtest tval','.mat'];
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];

    file = [name,labelnode,GRname,'PAIREDTtest tvalp',num2str(alpha_threshold),'.mat'];


    tmp =  tval.*double(pval<alpha_threshold);;
    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];

    file = [name,labelnode,GRname,'PAIREDTtest tvalFDR',num2str(alpha_threshold),'.mat'];
    tmp =   tval.*double(Q<alpha_threshold);
    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];


    file = [name,labelnode,GRname,'PAIREDTtest pval','.mat'];
    tmp =  pval;
    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];


    if isfield(job.c_statmatrix.b_PairedTtest.c_statpermutation,'b_permutation') %permutation here
        iduse = [idG1;idG2];
        nperm = str2num(job.c_statmatrix.b_PairedTtest.c_statpermutation.b_permutation.e_npermutation);
        for iperm=1:nperm
            permfix(iperm,:) =  randperm(numel(iduse));
        end
        disp(['Running ', num2str(nperm),' permutations on ', num2str(size(halfMAT,2)),' link :'])
        %loopslow to be improve
        for ilink = 1:size(halfMAT,2)
            fprintf('%s ',num2str(ilink));
            if mod(ilink,30)==0
                fprintf('\n ');
            end
            try
                idgood = find(~isnan(halfMAT(idG1,ilink))&~isnan(halfMAT(idG2,ilink)));
                [h,p,ci,stats] = ttest(halfMAT(idG1(idgood),ilink), halfMAT(idG2(idgood),ilink));
                try; cohend(:,ilink) = computeCohen_d(halfMAT(idG1(idgood),ilink), halfMAT(idG2(idgood),ilink), 'paired');catch; end;
                pval(:,ilink) = p;
                tval(:,ilink) = stats.tstat;
                for idperm=1:nperm
                    GR = iduse(permfix(idperm,:));
                    idgood = find(~isnan(halfMAT(idG1,ilink))&~isnan(halfMAT(idG2,ilink)));
                    [h,p,ci,stats] = ttest(halfMAT(GR(idG1(idgood)),ilink), halfMAT( GR(idG2(idgood)),ilink));
                    tvalperm(:,idperm,ilink) =stats.tstat;
                end
            catch
                fprintf('error link ');
            end
        end
        t_dist_all =  tvalperm ;
        t_dist_all(find(isnan(tvalperm)))=[];
        [fecdf,xecdf] = ecdf( t_dist_all);
        clustercritical=xecdf(sum(fecdf<(1-alpha_threshold)));
        if isfield(job.c_statmatrix.b_PairedTtest.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting,'b_MCT_ClusterBased')
            statrand = permute(tvalperm(1,:,:), [3 2 1]);
            statobs = permute(tval(1,:),[2,1]);
            tmp = randperm(numel(statrand));
            statrand =reshape(statrand(tmp),size(statrand,1),size(statrand,2));
            neighbourdist = job.c_statmatrix.b_PairedTtest.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_neighbourdist;
            minnbchan = job.c_statmatrix.b_PairedTtest.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_minnbchan;

            disp(['critical for clustering p=', num2str(alpha_threshold),' T=',num2str(clustercritical)])

            [stat, matneig] = FindClusterBasedPermutationInMatrix(chanpos, neighbourdist,clustercritical, statobs, statrand,minnbchan );
            disp(['Save cluster', fullfile(dir1,'Pairedttest_stat.mat') ])
            if isfield(stat,'posclusterslabelmat')
                for iposcluster = 1:numel(stat.posclusters)
                    file = [name,labelnode,num2str(nperm),' posCluster', num2str(iposcluster),' MCp=',num2str(stat.posclusters(iposcluster).prob),'.mat'];
                    tmp=(stat.posclusterslabelmat==iposcluster);
                    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                    new = [{dir1},{file}, {ZONEid},{0} ];
                    infonew = [infonew;new];
                end
            else
                disp('No positive cluster found')
            end
            if isfield(stat,'negclusterslabelmat') %technicaly no neg cluster in f value...
                for inegcluster = 1:numel(stat.negclusters)
                    file = [name,labelnode,num2str(nperm),' negCluster', num2str(iposcluster),' MCp=',num2str(stat.negclusters(inegcluster).prob),'.mat'];
                    tmp=(stat.negclusterslabelmat==inegcluster);
                    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                    new = [{dir1},{file}, {ZONEid},{0} ];
                    infonew = [infonew;new];
                end
            else
                disp('No negative cluster found')
            end

        end

        try %optionnal cohend
            file = [name,labelnode,'Cohend',num2str(alpha_threshold),'.mat'];
            tmp=squeeze(cohend(:));
            meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{1} ];
            infonew = [infonew;new];

        catch
            disp('Error computeCohen_d')
            disp('Cohen_d function must be install in the set path see:  Ruggero G. Bettinardi (2022). computeCohen_d(x1, x2, varargin) (https://www.mathworks.com/matlabcentral/fileexchange/62957-computecohen_d-x1-x2-varargin), MATLAB Central File Exchange. Retrieved November 16, 2022. ')
            p = mfilename('fullpath');
            [pathname,filename,ext]=fileparts(p);
            disp(['Add the fonction in matlab path: addpath(''',pathname, '\External\computeCohen_d'')'])
            disp('And verify if Fisher transform is not apply twice, data must be real')

        end
    end
    try
        if ~strcmp(fullfile(info{isubject,1}, ZONEid),fullfile(dir1,  ZONEid))
            copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid));
        end
    catch
        disp(['Error copy file:', fullfile(info{isubject,1}, ZONEid)])
    end
    %  dir1 = job.e_statmatrixPath{1};
    try
        if ismac
            % Code to run on Mac platform problem with xlswrite
            writetxtfile(fullfile(dir1,[name,labelnode,GRname,'PAIREDTtest.txt']),infonew);
            disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'PAIREDTtest.txt'])]);
        else
            try
                xlswrite(fullfile(dir1,[name,labelnode,GRname,'PAIREDTtest.xlsx']),infonew);
                disp(['Result .xlsx file saved: ' fullfile(dir1,[name,labelnode,'PAIREDTtest.xlsx'])]);
            catch
                writetxtfile(fullfile(dir1,[name,labelnode,GRname,'PAIREDTtest.txt']),infonew);
                disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,GRname,'PAIREDTtest.txt'])]);
            end
        end
    catch
        disp(['Error could not save .xlsx file: ' fullfile(dir1,[name,labelnode,GRname,'PAIREDTtest.xlsx'])]);
    end

elseif isfield(job.c_statmatrix,'b_zscore_Mat')
    AllC = [];
    id =1;
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    %Use a specific groupe
    idG1 = find(  groupeall==job.c_statmatrix.b_zscore_Mat.e_TtestOneSampleGR);
    idG2 = find(  groupeall==job.c_statmatrix.b_zscore_Mat.e_GRzscore);
    MATallG1 = MATall( idG1,:,:);
    MATallG2 = MATall( idG2,:,:);
    meanG1 = squeeze(nanmean(  MATallG1 ,1));
    stdG1 = squeeze(nanstd(  MATallG1 ,1));
    %savein result folder
    dir1 = job.e_statmatrixPath{1} ;

    for ifile=1:size(MATallG2)
        file = [name,'Zscore_',sprintf('%03.0f',(ifile)),list_subject{idG2(ifile)}  ,'.mat'];
        meancorr =  (squeeze(MATallG2(ifile,:,:))- meanG1)./ stdG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{2} ];
        infonew = [infonew;new];
    end
    if ~strcmp(fullfile(info{isubject,1}, ZONEid),fullfile(dir1,  ZONEid))
        copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid));
    end
    %  dir1 = job.e_statmatrixPath{1};
    try
        if ismac
            % Code to run on Mac platform problem with xlswrite
            writetxtfile(fullfile(dir1,[name,labelnode,'zscore.txt']),infonew);
            disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'zscore.txt'])]);
        else
            try
                xlswrite(fullfile(dir1,[name,labelnode,'zscore.xlsx']),infonew);
                disp(['Result .xlsx file saved: ' fullfile(dir1,[name,labelnode,'zscore.xlsx'])]);
            catch
                writetxtfile(fullfile(dir1,[name,labelnode,'zscore.txt']),infonew);
                disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'zscore.txt'])]);
            end
        end
    catch
        disp(['Error could not save .xlsx file: ' fullfile(dir1,[name,labelnode,'zscore.xlsx'])]);
    end
elseif isfield(job.c_statmatrix,'b_anova1_Mat')
    AllC = [];
    id =1;
    dir1 = job.e_statmatrixPath{1} ;
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    GROUPELIST = str2num(job.c_statmatrix.b_anova1_Mat.e_Anova1GR);
    %Use a specific groupe
    iduse = find(  sum(groupeall==GROUPELIST,2));
    if isfield(job.c_statmatrix.b_anova1_Mat.c_statpermutation,'b_Nopermutation')
        for i= 1:size(MATall,2)
            for j = 1:size(MATall,3)
                if i~=j

                    if i==1&j==2 %initialized once matrice results
                        isgood = find(~isnan(MATall( iduse,i,j)));
                        [p,tbl,stats] = anova1((MATall( iduse(isgood),i,j)), groupeall(iduse(isgood)),'off');
                        panova(i,j) = p;

                        [c,m,h,gnames] = multcompare(stats,'display','off');
                        diff = zeros(size(c,1),size(MATall,2),size(MATall,3) );
                        pdiff = zeros(size(c,1),size(MATall,2),size(MATall,3) );
                        for icomp=1:size(c,1)
                            labelmultcompare{icomp} = [gnames{c(icomp,1)}, '-' gnames{c(icomp,2)}];
                            labelmultcompareneg{icomp}  = [gnames{c(icomp,2)}, '-' gnames{c(icomp,1)}];
                        end
                    end
                    try
                        isgood = find(~isnan(MATall( iduse,i,j)));
                        [p,tbl,stats] = anova1((MATall( iduse(isgood),i,j)), groupeall(iduse(isgood)),'off');
                        try
                            sizeeffect = mes1way((MATall( iduse(isgood),i,j)),'eta2','group',groupeall(iduse(isgood)));
                            eta2(i,j) = sizeeffect.eta2;
                        catch
                        end
                        [c,m,h,gnames] = multcompare(stats,'display','off');
                        for icomp=1:size(c,1)
                            eval(['diff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=',num2str(c(icomp,4)),';']);
                            eval(['pdiff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=',num2str(c(icomp,6)),';']);
                        end
                        panova(i,j) = p;
                        nanova(i,j) = numel(iduse(isgood));
                        if p<alpha_threshold
                            filename = fullfile(dir1,['STATS_Anova ','GR ',num2str(GROUPELIST), ' ', ZoneList{i},' to ' ZoneList{j},' p= ',num2str(p),'.mat']);
                            save(filename,'stats');
                        end
                    catch %if pcould not be
                        for icomp=1:size(c,1)
                            eval(['diff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=nan;']);
                            eval(['pdiff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=','nan',';']);
                        end
                        panova(i,j) = 1;
                        nanova(i,j) = numel(iduse(isgood));
                    end
                else
                    panova(i,j) = nan;
                end
            end
        end

        [FRD,Q] = mafdr(panova(:)','LAMBDA',[0.0001:0.01:0.95]);
        Q = reshape(Q,size(panova));
        %SAVE FDR
        file = [name,labelnode,'ANOVAFDR p',num2str(alpha_threshold),'.mat'];
        meancorr = double( Q<alpha_threshold);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        for igr = 1:numel(GROUPELIST)
            iduse = find(  sum(groupeall==GROUPELIST(igr),2));
            file = [name,labelnode,'meanGR',num2str(GROUPELIST(igr)),'.mat'];
            meanG1 = squeeze(nanmean(MATall(iduse,:,:) ,1));
            meancorr = meanG1;
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
            infonew = [infonew;new];
        end

        for igr = 1:numel(GROUPELIST)
            iduse = find(  sum(groupeall==GROUPELIST(igr),2));
            file = [name,labelnode,'meanGR',num2str(GROUPELIST(igr)),'p',num2str(alpha_threshold),'.mat'];
            meanG1 = squeeze(nanmean(MATall(iduse,:,:) ,1)).*double(panova<alpha_threshold);
            meancorr = meanG1;
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
            infonew = [infonew;new];
        end


        for icomp = 1:size(diff,1)
            file = [name,labelnode,'diffGR',labelmultcompare{icomp},'.mat'];
            meanG1 = squeeze(diff(icomp,:,:));
            meancorr = meanG1;
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
            infonew = [infonew;new];
            file = [name,labelnode,'diffGR',labelmultcompareneg{icomp},'.mat'];
            meanG1 = squeeze(diff(icomp,:,:));
            meancorr = -meanG1;
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
            infonew = [infonew;new];

            file = [name,labelnode,'diffGR',labelmultcompare{icomp},'p',num2str(alpha_threshold),'.mat'];
            meanG1 = squeeze(diff(icomp,:,:)).*double(panova<alpha_threshold).*double(squeeze(pdiff(icomp,:,:)) <alpha_threshold);
            meancorr = meanG1;
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
            infonew = [infonew;new];
            file = [name,labelnode,'diffGR',labelmultcompareneg{icomp},'p',num2str(alpha_threshold),'.mat'];
            meanG1 = squeeze(diff(icomp,:,:)).*double(panova<alpha_threshold).*double(squeeze(pdiff(icomp,:,:)) <alpha_threshold);
            meancorr = -meanG1;
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
            infonew = [infonew;new];
        end
        file = [name,labelnode,'GR ',num2str(GROUPELIST),' pval',num2str(alpha_threshold),'.mat'];
        meancorr = double(panova);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        file = [name,labelnode,'GR ',num2str(GROUPELIST),' 1-pval',num2str(alpha_threshold),'.mat'];
        meancorr = double(1-panova);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];


        file = [name,labelnode,'GR ',num2str(GROUPELIST),' n valid measure ','.mat'];
        meancorr = double(nanova);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];


        try
            file = [name,labelnode,'GR ',num2str(GROUPELIST),' eta2 sizeeffect ','.mat'];
            meancorr = eta2;
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{0} ];
            infonew = [infonew;new];
        catch
            disp('Error size effect eta2')
            disp('mes1way function must be install in the set path see:  Harald Hentschke (2022). hhentschke/measures-of-effect-size-toolbox (https://github.com/hhentschke/measures-of-effect-size-toolbox), GitHub. Retrieved November 16, 2022. ')
            p = mfilename('fullpath');
            [pathname,filename,ext]=fileparts(p);
            disp(['Add the fonction in matlab path: addpath(''',pathname, '\External\EffectSizeToolbox'')'])
            disp('And verify if Fisher transform is not apply twice, data must be real')
        end
    elseif isfield(job.c_statmatrix.b_anova1_Mat.c_statpermutation,'b_permutation')     %not optimise in a loop for now... need to be improve when we have time
        nperm= str2num(job.c_statmatrix.b_anova1_Mat.c_statpermutation.b_permutation.e_npermutation);
        halfMAT = MATall(:,idhalf);
        f_dist = nan(nperm,numel(idhalf));
        fanova= nan(1,numel(idhalf));
        panova= nan(1,numel(idhalf));
        nanova= nan(1,numel(idhalf));
        for iperm=1:nperm
            permfix(iperm,:) =  randperm(numel(iduse));
        end
        disp(['Running ', num2str(nperm),' permutations on ', num2str(size(halfMAT,2)),' link :'])
        for ilink = 1:size(halfMAT,2)
            try
                fprintf('%s ',num2str(ilink));
                if mod(ilink,30)==0
                    fprintf('\n ');
                end
                isgood = find(~isnan(halfMAT( iduse,ilink)));
                [p,tbl,stats] = anova1((halfMAT( iduse(isgood),ilink)), groupeall(iduse(isgood)),'off');
                fanova(ilink) =tbl{2,5};
                panova(ilink) = p;
                nanova(ilink) = numel(iduse(isgood));
                for idperm=1:nperm
                    Yval = halfMAT(iduse(isgood),ilink);
                    [p,tbl,stats] = anova1(Yval, groupeall(permfix(idperm,isgood)),'off');
                    f_dist(idperm,ilink) = tbl{2,5};
                end
            catch
                disp(['error link', num2str(ilink)])
            end
        end
        [fecdf,xecdf] = ecdf(f_dist(:));
        disp( [num2str(nperm) ' Permutation p<',num2str(alpha_threshold), ' give a marginal F value=',num2str(xecdf(sum(fecdf<(1-alpha_threshold))))])

        %check whole distribution
        for ilink =1:size(halfMAT,2)
            if sum(xecdf<=fanova(ilink))==0
                pperm(ilink)= 1;
            else
                pperm(ilink) = 1-fecdf(sum(xecdf<=fanova(ilink)));
            end
        end
        clustercritical=xecdf(sum(fecdf<(1-alpha_threshold)));
        if isfield(job.c_statmatrix.b_anova1_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting,'b_MCT_ClusterBased')
            statrand = permute(f_dist, [2 1]);
            statobs = permute(fanova,[2,1]);
            tmp = randperm(numel(statrand));
            statrand =reshape(statrand(tmp),size(statrand,1),size(statrand,2));
            neighbourdist = job.c_statmatrix.b_anova1_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_neighbourdist;
            minnbchan = job.c_statmatrix.b_anova1_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_minnbchan;
            disp(['critical for clustering p=', num2str(alpha_threshold),' F=',num2str(clustercritical)])
            [stat, matneig] = FindClusterBasedPermutationInMatrix(chanpos, neighbourdist,clustercritical, statobs, statrand,minnbchan);
            disp(['Save cluster', fullfile(dir1,'Anova1_stat.mat') ])

            if isfield(stat,'posclusterslabelmat')
                for iposcluster = 1:numel(stat.posclusters)
                    file = [name,labelnode,num2str(nperm),' posCluster', num2str(iposcluster),' MCp=',num2str(stat.posclusters(iposcluster).prob),'.mat'];
                    tmp=(stat.posclusterslabelmat==iposcluster);
                    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                    new = [{dir1},{file}, {ZONEid},{0} ];
                    infonew = [infonew;new];
                end
            else
                disp('No positive cluster found')
            end
            if isfield(stat,'negclusterslabelmat') %technicaly no neg cluster in f value...
                for inegcluster = 1:numel(stat.negclusters)
                    file = [name,labelnode,num2str(nperm),' negCluster', num2str(iposcluster),' MCp=',num2str(stat.negclusters(inegcluster).prob),'.mat'];
                    tmp=(stat.negclusterslabelmat==inegcluster);
                    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                    new = [{dir1},{file}, {ZONEid},{0} ];
                    infonew = [infonew;new];
                end
            else
                disp('No negative cluster found')
            end

        end

        %equivalent tmax Lage-Castellanos, A., Mart�nez-Montes, E., Hern�ndez-Cabrera, J.A., Gal�n, L., 2010. False discovery rate and permutation test: An evaluation in ERP data analysis. Statistics in Medicine 29, 63�74. https://doi.org/10.1002/sim.3784
        supvalue =reshape(f_dist, nperm,  size(f_dist,2));
        dist_sup = max( supvalue');
        disp( [' Max Permutation over all link p<',num2str(alpha_threshold), ' give a marginal F value=',num2str(xecdf(sum(fecdf<=(1-alpha_threshold))))])

        [fecdf,xecdf] = ecdf( [dist_sup]);
        disp( [' Max Permutation over all link p<',num2str(alpha_threshold), ' give a marginal F value=',num2str(xecdf(sum(fecdf<=(1-alpha_threshold))))])
        for  ilink = 1:size(halfMAT,2)
            if sum(xecdf<=fanova(ilink))==0
                ppermmax(ilink)=1;
            else
                ppermmax(ilink) = 1-fecdf(sum(xecdf<=fanova(ilink)));
            end
        end
        try
            [FRD,Q] = mafdr(panova(:)','LAMBDA',[0.0001:0.01:0.95]);
            Q = reshape(Q,size(panova));
            %SAVE FDR
            file = [name,labelnode,'Anova_FDR Qval.mat'];
            tmp=double(Q);
            meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{0} ];
            infonew = [infonew;new];
        catch
        end

        %SAVE p anova
        file = [name,labelnode,'Anova_pval.mat'];
        tmp=double(panova);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        %SAVE p prem
        file = [name,labelnode,'Anova_pperm', num2str(nperm),'.mat'];
        tmp= double(pperm);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        file = [name,labelnode,'Anova_ppermmax', num2str(nperm),'.mat']; %problem possible avec le ppermmax plusieurs groupe lire +
        tmp= double(ppermmax);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        file = [name,labelnode,'Anova_pval 1-p.mat'];
        tmp= double(1-panova);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        file = [name,labelnode,'Anova_pperm', num2str(nperm),' 1-p.mat'];
        tmp =  double(1-pperm);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        file = [name,labelnode,'Anova_ppermmax', num2str(nperm),' 1-p.mat']; %problem possible avec le ppermmax plusieurs groupe lire +
        tmp =  double(1-ppermmax);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        file = [name,labelnode,'Anova_F','.mat'];
        tmp =  double(fanova);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    end

    try
    if ~strcmp(fullfile(info{isubject,1}, ZONEid),fullfile(dir1,  ZONEid))
        copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid))
    end
    catch
    end
    

    if ismac
        % Code to run on Mac platform problem with xlswrite
        writetxtfile(fullfile(dir1,[name,labelnode,'GR ',num2str(GROUPELIST),' Anova.txt']),infonew);
        disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'GR ',num2str(GROUPELIST),' Anova.txt'])]);
    else
        try
            xlswrite(fullfile(dir1,[name,labelnode,'GR ',num2str(GROUPELIST),' Anova.xlsx']),infonew);
            disp(['Result .xlsx file saved: ' fullfile(dir1,[name,labelnode,'GR ',num2str(GROUPELIST),' Anova.xlsx'])]);
        catch
            writetxtfile(fullfile(dir1,[name,labelnode,'GR ',num2str(GROUPELIST),' Anova.txt']),infonew);
            disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'GR ',num2str(GROUPELIST),' Anova.txt'])]);
        end
    end
elseif isfield(job.c_statmatrix,'b_ANCOVA_Mat') %version de matlab juste une covariable abandonn� remplac� par mancovan
    dir1 = job.e_statmatrixPath{1};
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    covariableall=[];
    covariablestring = job.c_statmatrix.b_ANCOVA_Mat.b_Covariable_Mat;
    [token,remain] =strtok(covariablestring,',');
    covariableall =  [covariableall,{strtrim(token)}];
    while ~isempty(remain)
        [token,remain] =strtok(remain,',');
        covariableall =  [covariableall,{strtrim(token)}];
    end

    for icov = 1:numel(covariableall)
        notfoundstophere = 1;
        Pearsony = covariableall{icov};
        ycol = 0;
        for icol=1:size(info,2)
            if ~isnan(info{1,icol})
                if strcmp(strtrim(upper(deblank(info{1,icol}))), strtrim(upper(Pearsony)));
                    ycol = icol;
                end
            end
        end
        if ycol
            id= 1;
            for i=2:size(info,1)
                if isnumeric(info{i,ycol })
                    score(id,icov) = info{i,ycol };
                    id = id+1;
                else
                    score(id,icov) = nan;
                    id = id+1;
                end
            end

        else
            disp([Pearsony,' column not found, covariable unavailable'])
            out='Stat tests';
            return
        end
    end
    AllC = [];
    id =1;
    dir1 = job.e_statmatrixPath{1} ;
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    GROUPELIST = str2num(job.c_statmatrix.b_ANCOVA_Mat.e_Anova1GR);
    %Use a specific groupe
    iduse = find(  sum(groupeall==GROUPELIST,2)& ~sum(isnan(score),2));

    for i= 1:size(MATall,2)
        for j = 1:size(MATall,3)
            if i~=j
                isgood = find(~isnan(MATall( iduse,i,j)));
                try
                    [h,atab,ctab,stats] = aoctool(score(iduse(isgood),:), (MATall( iduse(isgood),i,j)), groupeall(iduse(isgood)),alpha_threshold,'x','y','Groupe','off');

                    if 0
                        [h,atab,ctab,stats] = aoctoolpermutation(score(iduse(isgood)), (MATall( iduse(isgood),i,j)), groupeall(iduse(isgood)),500);
                    end
                    panovaGroupe(i,j)= atab{2,end};
                    panovaCoV(i,j)= atab{3,end};
                    panovaGroupebyCoV(i,j)= atab{4,end};
                    nanova(i,j) = numel(iduse(isgood));
                    p_groupe = atab{2,6};
                    p_CoV= atab{3,6};
                    if i==1&j==2 %initialized once matrice results
                        [c,m,h,gnames] = multcompare(stats,'display','off');

                        diff = zeros(size(c,1),size(MATall,2),size(MATall,3) );
                        pdiff = zeros(size(c,1),size(MATall,2),size(MATall,3) );
                        for icomp=1:size(c,1)
                            labelmultcompare{icomp} = [gnames{c(icomp,1)}, '-' gnames{c(icomp,2)}];
                            labelmultcompareneg{icomp}  = [gnames{c(icomp,2)}, '-' gnames{c(icomp,1)}];
                        end
                    end
                    [c,m,h,gnames] = multcompare(stats,'display','off');
                    for icomp=1:size(c,1)
                        eval(['diff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=',num2str(c(icomp,4)),';']);
                        eval(['pdiff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=',num2str(c(icomp,6)),';']);
                    end
                    if p_groupe < alpha_threshold %p groupe
                        filename = fullfile(dir1,['STATS_Ancova ','GR ',num2str(GROUPELIST),' x groupe' ZoneList{i},' to ' ZoneList{i},' p= ',num2str(p_groupe),'.mat']);
                        save(filename,'stats');
                        % multcompare(stats)
                    end
                    if p_CoV < alpha_threshold %p groupe
                        filename = fullfile(dir1,['STATS_Ancova ','GR ',num2str(GROUPELIST),' x COV ' ZoneList{i},' to ' ZoneList{i},' p= ',num2str(p_CoV),'.mat']);
                        save(filename,'stats');
                        %  multcompare(stats)

                    end

                catch
                    panovaGroupe(i,j) = 1;
                    panovaCoV(i,j) = 1;
                    panovaGroupebyCoV(i,j) = 1;
                    nanova(i,j) = numel(iduse(isgood));
                    try
                        for icomp=1:size(c,1)
                            eval(['diff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=','nan',';']);
                            eval(['pdiff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=','nan',';']);
                        end
                        disp(['Error stat for link ' num2str(i) 'to ',num2str(j) ])
                    catch
                        disp(['Error stat for link ' num2str(i) 'to ',num2str(j) ])
                    end
                end
            else
                panovaGroupe(i,j) = nan;
                panovaCoV(i,j) = nan;
                panovaGroupebyCoV(i,j) = nan;
            end

        end
    end
    [FRD,Q] = mafdr(panovaGroupe(:));
    Q = reshape(Q,size(panovaGroupe));
    %SAVE FDR
    file = [name,labelnode,'ANCOVAFDR groupe p',num2str(alpha_threshold),'.mat'];
    meancorr = double( Q<alpha_threshold);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];

    [FRD,Q] = mafdr( panovaCoV(:));
    Q = reshape(Q,size( panovaCoV));
    %figure;hist(panovaCoV(:))
    %SAVE FDR
    file = [name,labelnode,'ANCOVAFDR COV p',num2str(alpha_threshold),'.mat'];
    meancorr = double( Q<alpha_threshold);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];

    [FRD,Q] = mafdr(panovaGroupebyCoV(:));
    Q = reshape(Q,size(panovaGroupebyCoV));
    %SAVE FDR
    file = [name,labelnode,'ANCOVAFDR groupe p',num2str(alpha_threshold),'.mat'];
    meancorr = double( Q<alpha_threshold);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];


    for igr = 1:numel(GROUPELIST)
        iduse = find(  sum(groupeall==GROUPELIST(igr),2));
        file = [name,labelnode,'meanGR',num2str(GROUPELIST(igr)),'.mat'];
        meanG1 = squeeze(nanmean(MATall(iduse,:,:) ,1));
        meancorr = meanG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
        infonew = [infonew;new];
    end
    for icomp = 1:size(diff,1)
        file = [name,labelnode,'Ancova diffGR',labelmultcompare{icomp},'.mat'];
        meanG1 = squeeze(diff(icomp,:,:)).*double(panovaGroupe<=alpha_threshold);
        meancorr = meanG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
        infonew = [infonew;new];
        file = [name,labelnode,'Ancova diffGR',labelmultcompareneg{icomp},'.mat'];
        meanG1 = squeeze(diff(icomp,:,:)).*double(panovaGroupe<=alpha_threshold);
        meancorr = -meanG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
        infonew = [infonew;new];

        file = [name,labelnode,'Ancova diffGR',labelmultcompare{icomp},'p',num2str(alpha_threshold),'.mat'];
        meanG1 = squeeze(diff(icomp,:,:)).*double(panovaGroupe<alpha_threshold).*double(squeeze(pdiff(icomp,:,:)) <alpha_threshold);
        meancorr = meanG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
        infonew = [infonew;new];
        file = [name,labelnode,'Ancova diffGR',labelmultcompareneg{icomp},'p',num2str(alpha_threshold),'.mat'];
        meanG1 = squeeze(diff(icomp,:,:)).*double(panovaGroupe<alpha_threshold).*double(squeeze(pdiff(icomp,:,:)) <alpha_threshold);
        meancorr = -meanG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
        infonew = [infonew;new];
    end
    file = [name,labelnode,'Ancova GR ',num2str(GROUPELIST),', pGroupe',num2str(alpha_threshold),'.mat'];
    meancorr = double(panovaGroupe<alpha_threshold);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];


    file = [name,labelnode,'Ancova GR ',num2str(GROUPELIST), ', pCoV',num2str(alpha_threshold),'.mat'];
    meancorr = double(panovaCoV<alpha_threshold);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];

    file = [name,labelnode,'Ancova GR ',num2str(GROUPELIST), ', pGroupebyCoV',num2str(alpha_threshold),'.mat'];
    meancorr = double(panovaGroupebyCoV<alpha_threshold);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];
    if ~strcmp(fullfile(info{isubject,1}, ZONEid),fullfile(dir1,  ZONEid))
        copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid))
    end
    if ismac
        % Code to run on Mac platform problem with xlswrite
        writetxtfile(fullfile(dir1,[name,labelnode,'Ancova', 'GR ',num2str(GROUPELIST), '.txt']),infonew);
        disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'Ancova', 'GR ',num2str(GROUPELIST), '.txt'])]);
    else
        try
            xlswrite(fullfile(dir1,[name,labelnode,'Ancova', ' GR ',num2str(GROUPELIST), '.xlsx']),infonew);
            disp(['Result .xlsx file saved: ' fullfile(dir1,[name,labelnode,'Ancova', ' GR ',num2str(GROUPELIST), '.xlsx'])]);
        catch
            writetxtfile(fullfile(dir1,[name,labelnode,'Ancova', ' GR ',num2str(GROUPELIST), '.txt']),infonew);
            disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'Ancova', ' GR ',num2str(GROUPELIST), '.txt'])]);
        end
    end



elseif isfield(job.c_statmatrix,'b_anovarep_Mat')
    dir1 = job.e_statmatrixPath{1};
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    covariableall=[];
    covariablestring = job.c_statmatrix.b_anovarep_Mat.e_Anovarep_predictor;
    [token,remain] =strtok(covariablestring,',');
    covariableall =  [covariableall,{strtrim(token)}];
    while ~isempty(remain)
        [token,remain] =strtok(remain,',');
        covariableall =  [covariableall,{strtrim(token)}];
    end
 
    for icov = 1:numel(covariableall)
        notfoundstophere = 1;
        Pearsony = covariableall{icov};
        ycol = 0;
        for icol=1:size(info,2)
            if ~isnan(info{1,icol})
                if strcmp(strtrim(upper(deblank(info{1,icol}))), strtrim(upper(Pearsony)));
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
    GROUPELIST = str2num(job.c_statmatrix.b_anovarep_Mat.e_Anovarep_withinsubjectGR);
    if isempty(GROUPELIST)
        disp('Ensure to enter numerical value for groupe id to include')
    end
    for igr = 1:numel(GROUPELIST)
        iduse = find(  sum(groupeall==GROUPELIST(igr),2));
        file = [name,labelnode,'meanGR',num2str(igr),'.mat'];
        meanG1 = squeeze(nanmean(MATall(iduse,:,:) ,1));

        meancorr = meanG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
        infonew = [infonew;new];
    end


    for i= 1:size(MATall,2)
        disp(num2str(i));
        for j = 1:size(MATall,3)
           % if i~=j
            %    try
                    T1= table();
                    %Use a specific groupe
                    for igr = 1:numel(GROUPELIST)
                        iduse = find(  sum(groupeall==GROUPELIST(igr),2));
                        GR= MATall( iduse,i,j);
                        if isempty(T1)
                            T1= table(GR,'VariableNames',{['G',num2str(igr)]});
                        else
                            if size(T1,1)~=size(GR,1)
                                disp(['Error: repeated measure needs an equal number of measurement. First group has ', num2str(size(T1,1)), ' measurements while the other has ', num2str(size(GR,1)), ' measurements.' ]);
                                disp('Stats failed');
                                return
                            else
                                T1 = addvars(T1,GR,'NewVariableNames',{['G',num2str(igr)]});
                            end
                        end
                    end

 
                    % Add Predictor measure
                    idvariable = 1;
                    for idvariable = 1:numel(covariableall)
                        X = score(iduse,idvariable);
                        T1 = addvars(T1,num2str(X),'NewVariableNames',covariableall{idvariable});
                    end
                    %job.c_statmatrix.b_anovarep_Mat.e_anovarep_model
                    rm = fitrm(T1,job.c_statmatrix.b_anovarep_Mat.e_anovarep_model);
                   %  rm = fitrm(T1,'G1,G2 ~ Age + PREMA');
                     ranovatbl = anova(rm);
 
    
                    Fanova = table2array(ranovatbl(:,6));
                    Panova = table2array(ranovatbl(:,7));
                    tbl = multcompare(rm,covariableall{1});
                    for itbl = 1:size(tbl,1)
                        DiffGR(itbl,i,j)= tbl.Difference(itbl);
                        PvalGR(itbl,i,j)= tbl.pValue(itbl);
                        tmp1 = tbl{itbl,1};
                        tmp2 = tbl{itbl,2};
                        nameGR{itbl} = [covariableall{1},tmp1{1},'-',tmp2{1}];
                    end

                    % tbl = multcompare(rm,covariableall{2});
                    % for itbl = 1:size(tbl,1)
                    %     DiffAGE(itbl,i,j)= tbl.Difference(itbl);
                    %     PvalAGE(itbl,i,j)= tbl.pValue(itbl);
                    %     tmp1 = tbl{itbl,1};
                    %     tmp2 = tbl{itbl,2};
                    %     nameGR{itbl} = [covariableall{1},tmp1{1},'-',tmp2{1}];
                    % end


                    tbl = multcompare(rm,'Time');
                    for itbl = 1:size(tbl,1)
                        DiffTime(itbl,i,j)= tbl.Difference(itbl);
                        PvalTime(itbl,i,j)= tbl.pValue(itbl);
                        tmp1 = num2str(tbl{itbl,1});
                        tmp2 = num2str(tbl{itbl,2});
                        nameTime{itbl} = ['Time',tmp1,'-',tmp2];
                    end

                    for ipval = 1:numel(Panova)
                        eval(['panova',num2str(ipval),'(',num2str(i),',',num2str(j),')=  Panova(',num2str(ipval),');']);
                    end
                    for ifval = 1:numel(Fanova)
                        eval(['fanova',num2str(ifval),'(',num2str(i),',',num2str(j),')=  Fanova(',num2str(ifval),');']);
                    end

                % catch
                %     disp('Error fitrm')
                % 
                % end
            %end
        end
    end


    for ifval = 1:numel(Panova)
        file = [name,labelnode,' fanova', num2str(ifval) ,'.mat'];
        meancorr = eval(['fanova', num2str(ifval)]);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    end


    for ipval = 1:numel(Panova)
        file = [name,labelnode,'1-panova', num2str(ipval),'.mat'];
        meancorr = eval(['1-panova', num2str(ipval)]);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    end

    for igr= 1:size(DiffGR,1)
        file = [name,labelnode,nameGR{igr},'diff.mat'];
        meancorr = squeeze(DiffGR(igr,:,:)) ;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        file = [name,labelnode,nameGR{igr},'diffp,',num2str(alpha_threshold),'.mat'];
        meancorr = squeeze(DiffGR(igr,:,:)).*double(squeeze(PvalGR(igr,:,:))<alpha_threshold) ;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    end


    for igr= 1:size(DiffTime,1)
        file = [name,labelnode,nameTime{igr},'diff.mat'];
        meancorr = squeeze(DiffTime(igr,:,:)) ;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
        file = [name,labelnode,nameTime{igr},'diffp,',num2str(alpha_threshold),'.mat'];

        meancorr = squeeze(DiffTime(igr,:,:)).*double(squeeze(PvalTime(igr,:,:))<alpha_threshold) ;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    end

    try
    if ~strcmp(fullfile(info{isubject,1}, ZONEid),fullfile(dir1,  ZONEid))
        copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid))
    end
    catch
    end

    if ismac
        % Code to run on Mac platform problem with xlswrite
        writetxtfile(fullfile(dir1,[name,labelnode,'RepAnova.txt']),infonew);
        disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'RepAnova.txt'])]);
    else
        try
            xlswrite(fullfile(dir1,[name,labelnode,'RepAnova.xlsx']),infonew);
            disp(['Result .xlsx file saved: ' fullfile(dir1,[name,labelnode,'RepAnova.xlsx'])]);
        catch
            writetxtfile(fullfile(dir1,[name,labelnode,'RepAnova.txt']),infonew);
            disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'RepAnova.txt'])]);
        end
    end
elseif isfield(job.c_statmatrix,'b_NWayANOVA_Mat')

     GROUPE =  str2num(job.c_statmatrix.b_NWayANOVA_Mat.e_NWayANOVA_model_GR);          
     iduse = find(sum(groupeall==GROUPE,2));
     eval(['formule=',job.c_statmatrix.b_NWayANOVA_Mat.e_NWayANOVA_model,';'])   


       formule_to_eval= [];
       varsname_to_eval= [];
        for i=1:numel(formule)
            formule_to_eval = [formule_to_eval, 'T.', formule{i},'(iduse),'] 
            varsname_to_eval = [ varsname_to_eval , '''',formule{i},''';']
        end
        formule_to_eval(end) = []; 
        varsname_to_eval(end)=[];
       

    dir1 = job.e_statmatrixPath{1};
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];

T = readtable(xlslistfile);


%g1 =find(sum(groupeall==1,2)) ;
%g2 =find(sum(groupeall==2,2)) ;
 for i= 1:size(MATall,2)
        disp(num2str(i));
        for j = 1:size(MATall,3)

        y = squeeze(MATall(:,i,j));
      
      % T.PREMA(iduse)
       eval(['[Panova,tbl2,stats,terms] = anovan(y(iduse),{', formule_to_eval,'},''model'',''full'',''varnames'',{',varsname_to_eval ,'},''display'',''off'');'])
        %[h,p,ci,stats] = ttest2(y(g1), y(g2));
        %tval(i,j) = stats.tstat;
        
       Panova = 1-Panova;                  
                    for ipval = 1:numel(Panova)
                        eval(['panova',num2str(ipval),'(',num2str(i),',',num2str(j),')=  Panova(',num2str(ipval),');']);
                    end
        end
 end

   for ipval = 1:numel(Panova)
       condition =tbl2{1+ipval,1};
       k = findstr(condition,':');
       if ~isempty(k)
       condition(k)=[' '];
       end
       file = [name,labelnode,'1-panova', num2str(ipval),condition, '.mat'];
        meancorr = eval(['panova', num2str(ipval)]);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        % file = [name,labelnode,'tval.mat'];
        % meancorr = tval;
        % save(fullfile(dir1,[file]),'ZoneList','meancorr');
        % disp(['Save: ', fullfile(dir1,[file])]);
        % new = [{dir1},{file}, {ZONEid},{0} ];
        % infonew = [infonew;new];
        % 
        % 
        % file = [name,labelnode,'-tval.mat'];
        % meancorr = -tval;
        % save(fullfile(dir1,[file]),'ZoneList','meancorr');
        % disp(['Save: ', fullfile(dir1,[file])]);
        % new = [{dir1},{file}, {ZONEid},{0} ];
        % infonew = [infonew;new];

   end

    if ismac
        % Code to run on Mac platform problem with xlswrite
        writetxtfile(fullfile(dir1,[name,labelnode,'NWAYAnova.txt']),infonew);
        disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'NWAYAnova.txt'])]);
    else
        try
            xlswrite(fullfile(dir1,[name,labelnode,'NWAYAnova.xlsx']),infonew);
            disp(['Result .xlsx file saved: ' fullfile(dir1,[name,labelnode,'NWAYAnova.xlsx'])]);
        catch
            writetxtfile(fullfile(dir1,[name,labelnode,'NWAYAnova.txt']),infonew);
            disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'NWAYAnova.txt'])]);
        end
    end

elseif isfield(job.c_statmatrix,'b_kruskalwallis_Mat')
    AllC = [];
    id =1;
    dir1 = job.e_statmatrixPath{1} ;
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    GROUPELIST = str2num(job.c_statmatrix.b_kruskalwallis_Mat.e_Anova1GR);
    %Use a specific groupe
    iduse = find(  sum(groupeall==GROUPELIST,2));

    for i= 1:size(MATall,2)
        for j = 1:size(MATall,3)
            if i~=j
                if i==1&j==2 %initialized once matrice results
                    isgood = find(~isnan(MATall( iduse,i,j)));
                    [p,tbl,stats] = kruskalwallis((MATall( iduse(isgood),i,j)), groupeall(iduse(isgood)),'off');
                    [c,m,h,gnames] = multcompare(stats,'display','off');
                    diff = zeros(size(c,1),size(MATall,2),size(MATall,3) );
                    pdiff = zeros(size(c,1),size(MATall,2),size(MATall,3) );
                    for icomp=1:size(c,1)
                        labelmultcompare{icomp} = [gnames{c(icomp,1)}, '-' gnames{c(icomp,2)}];
                        labelmultcompareneg{icomp}  = [gnames{c(icomp,2)}, '-' gnames{c(icomp,1)}];
                    end
                end
                try
                    isgood = find(~isnan(MATall( iduse,i,j)));
                    [p,tbl,stats] = kruskalwallis((MATall( iduse(isgood),i,j)), groupeall(iduse(isgood)),'off');
                    [c,m,h,gnames] = multcompare(stats,'display','off');
                    for icomp=1:size(c,1)
                        eval(['diff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=',num2str(c(icomp,4)),';']);
                        eval(['pdiff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=',num2str(c(icomp,6)),';']);
                    end
                    panova(i,j) = p;
                    nanova(i,j) = numel(iduse(isgood));
                    if p<alpha_threshold
                        if job.m_nodeunit ==1
                            filename = fullfile(dir1,['STATS_KruskalWallis ','GR ',num2str(GROUPELIST), ' ', ZoneList{i},' to ' ZoneList{j},' p= ',num2str(p),'.mat']);
                            save(filename,'stats');
                        else job.m_nodeunit ==2
                            filename = fullfile(dir1,['STATS_KruskalWallis ','GR ',num2str(GROUPELIST), ' ', ZoneLabel{i},' to ' ZoneLabel{j},' p= ',num2str(p),'.mat']);
                            save(filename,'stats')
                        end
                    end
                catch %if pcould not be
                    for icomp=1:size(c,1)
                        eval(['diff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=nan;']);
                        eval(['pdiff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=','nan',';']);
                    end
                    panova(i,j) = 1;
                    nanova(i,j) = numel(iduse(isgood));
                end
            else
                panova(i,j) = nan;
            end
        end
    end
    try
        [FRD,Q] = mafdr(panova(:)','LAMBDA',[0.0001:0.01:0.95]);
        Q = reshape(Q,size(panova));
        %SAVE FDR
        file = [name,labelnode,'KrustalWallisFDR Q','.mat'];
        meancorr = double(Q);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
        file = [name,labelnode,'KrustalWallisFDR Q 1-Q','.mat'];
        meancorr = double(1-Q);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    catch
    end

    if isfield(job.c_statmatrix.b_kruskalwallis_Mat.c_statpermutation,'b_permutation')
        % implementation max for each
        nperm= str2num(job.c_statmatrix.b_kruskalwallis_Mat.c_statpermutation.b_permutation.e_npermutation);
        halfMAT = MATall(:,idhalf);
        chi_dist_groupe = nan(nperm,numel(idhalf));
        for iperm=1:nperm
            permfix(iperm,:) =  randperm(numel(iduse));
        end
        disp(['Running ', num2str(nperm),' permutations on ', num2str(size(halfMAT,2)),' link :'])
        for ilink = 1:size(halfMAT,2)
            fprintf('%s ',num2str(ilink));
            if mod(ilink,30)==0
                fprintf('\n ');
            end
            isgood = find(~isnan(MATall( iduse,ilink)));
            [p,tbl,stats] = kruskalwallis((halfMAT( iduse(isgood),ilink)), groupeall(iduse(isgood)),'off');
            chi_groupe(ilink) = tbl{2,5};
            for idperm=1:nperm
                isgood = find(~isnan(halfMAT( iduse(permfix(idperm,:)),ilink)));
                [p,tbl,stats] = kruskalwallis((halfMAT( iduse(isgood),ilink)), groupeall(iduse(isgood)),'off');
                chi_dist_groupe(idperm,ilink) = tbl{2,5};
            end
        end

        [fecdf,xecdf] = ecdf(chi_dist_groupe(:));
        disp( [num2str(nperm) ' Permutation p<',num2str(alpha_threshold), ' give a marginal Chi value=',num2str(xecdf(sum(fecdf<(1-alpha_threshold))))])
        %figure;plot(xecdf,fecdf)
        clustercritical=xecdf(sum(fecdf<(1-alpha_threshold)));
        if isfield(job.c_statmatrix.b_kruskalwallis_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting,'b_MCT_ClusterBased')
            statrand = permute(chi_dist_groupe,[2,1]);
            statobs = permute(chi_groupe,[2,1]);
            tmp = randperm(numel(statrand));
            statrand =reshape(statrand(tmp),size(statrand,1),size(statrand,2));
            neighbourdist = job.c_statmatrix.b_kruskalwallis_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_neighbourdist;
            minnbchan = job.c_statmatrix.b_kruskalwallis_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_minnbchan;
            disp(['critical for clustering p=', num2str(alpha_threshold),' F=',num2str(clustercritical)])
            [stat, matneig] = FindClusterBasedPermutationInMatrix(chanpos, neighbourdist,clustercritical, statobs, statrand,minnbchan  );
            save(fullfile(dir1,['KruskalWallis_Cluster stat.mat']),'stat')
            disp(['Save cluster', fullfile(dir1,['KruskalWallis_Cluster stat.mat']) ])

            if isfield(stat,'posclusterslabelmat')
                for iposcluster = 1:numel(stat.posclusters)
                    file = [name,num2str(nperm),' posCluster', num2str(iposcluster),' MCp=',num2str(stat.posclusters(iposcluster).prob),'.mat'];
                    tmp=(stat.posclusterslabelmat==iposcluster);
                    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                    new = [{dir1},{file}, {ZONEid},{0} ];
                    infonew = [infonew;new];
                end
            else
                disp('No positive cluster found')
            end
            if isfield(stat,'negclusterslabelmat') %technicaly no neg cluster in f value...
                for inegcluster = 1:numel(stat.negclusters)
                    file = [name, num2str(nperm),' negCluster', num2str(iposcluster),' MCp=',num2str(stat.negclusters(inegcluster).prob),'.mat'];
                    tmp=(stat.negclusterslabelmat==inegcluster);
                    meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                    save(fullfile(dir1,[file]),'ZoneList','meancorr');
                    disp(['Save: ', fullfile(dir1,[file])]);
                    new = [{dir1},{file}, {ZONEid},{0} ];
                    infonew = [infonew;new];
                end
            else
                disp('No negative cluster found')
            end
        end

        %check p univarite permutation distribution
        for ilink=1:numel(chi_groupe)
            if sum(xecdf<=chi_groupe(ilink))==0
                pperm(ilink)= 1;
            else
                pperm(ilink) = 1-fecdf(sum(xecdf<=chi_groupe(ilink)));
            end
        end

        %equivalent tmax Lage-Castellanos, A., Mart�nez-Montes, E., Hern�ndez-Cabrera, J.A., Gal�n, L., 2010. False discovery rate and permutation test: An evaluation in ERP data analysis. Statistics in Medicine 29, 63�74. https://doi.org/10.1002/sim.3784
        supvalue =reshape(chi_dist_groupe, nperm,  size(chi_dist_groupe,2));
        dist_sup = max( supvalue');
        %          figure;
        %            hist(dist_sup(:),100)
        [fecdf,xecdf] = ecdf( [dist_sup]);
        %'Exal Prob chi2 permutation max';
        for ilink=1:numel(chi_groupe)
            if sum(xecdf<=chi_groupe(ilink))==0
                ppermmax(ilink)=1;
            else
                ppermmax(ilink) = 1-fecdf(sum(xecdf<=chi_groupe(ilink)));
            end
        end


        file = [name,labelnode,'KrustalWallis',num2str(nperm), 'Permp.mat'];
        tmp=double(pperm);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        file = [name,labelnode,'KrustalWallis',num2str(nperm), 'Permmaxp.mat'];
        tmp=double(ppermmax);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,file),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        file = [name,labelnode,'KrustalWallis',num2str(nperm), 'Perm 1-p.mat'];
        tmp= double(1-pperm);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        file = [name,labelnode,'KrustalWallis',num2str(nperm), 'Permmax 1-p.mat'];
        tmp= double( 1-ppermmax);
        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];


    end  % Fin permutation






    for igr = 1:numel(GROUPELIST)
        iduse = find(  sum(groupeall==GROUPELIST(igr),2));
        file = [name,labelnode,'meanGR',num2str(GROUPELIST(igr)),'.mat'];
        meanG1 = squeeze(nanmean(MATall(iduse,:,:) ,1));
        meancorr = meanG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
        infonew = [infonew;new];
    end

    for igr = 1:numel(GROUPELIST)
        iduse = find(  sum(groupeall==GROUPELIST(igr),2));
        file = [name,labelnode,'meanGR',num2str(GROUPELIST(igr)),'KrustalWallisp',num2str(alpha_threshold),'.mat'];
        meanG1 = squeeze(nanmean(MATall(iduse,:,:) ,1)).*double(panova<alpha_threshold);
        matcorr = meanG1;
        meancorr = meanG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
        infonew = [infonew;new];
    end

    try
        for icomp = 1:size(diff,1)
            file = [name,labelnode,'diff_rank_GR',labelmultcompare{icomp},'.mat'];
            meanG1 = squeeze(diff(icomp,:,:));
            meancorr = meanG1;
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
            infonew = [infonew;new];
            file = [name,labelnode,'diff_rank_GR',labelmultcompareneg{icomp},'.mat'];
            meanG1 = squeeze(diff(icomp,:,:));
            meancorr = -meanG1;
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
            infonew = [infonew;new];

            file = [name,labelnode,'diff_rank_GR',labelmultcompare{icomp},'p',num2str(alpha_threshold),'.mat'];
            meanG1 = squeeze(diff(icomp,:,:)).*double(panova<alpha_threshold).*double(squeeze(pdiff(icomp,:,:)) <alpha_threshold);
            meancorr = meanG1;
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
            infonew = [infonew;new];
            file = [name,labelnode,'diff_rank_GR',labelmultcompareneg{icomp},'p',num2str(alpha_threshold),'.mat'];
            meanG1 = squeeze(diff(icomp,:,:)).*double(panova<alpha_threshold).*double(squeeze(pdiff(icomp,:,:)) <alpha_threshold);
            meancorr = -meanG1;
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
            infonew = [infonew;new];
        end
    catch
    end
    file = [name,labelnode,'GR KruskalWallis',num2str(GROUPELIST),' pval',num2str(alpha_threshold),'.mat'];
    meancorr = double(panova);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];

    file = [name,labelnode,'GR KruskalWallis',num2str(GROUPELIST),' n valid measure ','.mat'];
    meancorr = double(nanova);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];
    if job.m_nodeunit==1
        try
            if ~strcmp(fullfile(info{isubject,1}, ZONEid),fullfile(dir1,  ZONEid))
                copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid))
                disp(['Copy file: ',  fullfile(info{isubject,1}, ZONEid) ' to ',fullfile(dir1,  ZONEid) ])
            end
        catch
            disp(['Could not copy file: ',  fullfile(info{isubject,1}, ZONEid)])
        end
    end
    if ismac
        % Code to run on Mac platform problem with xlswrite
        writetxtfile(fullfile(dir1,[name,labelnode,'GR ',num2str(GROUPELIST),' KruskalWallis.txt']),infonew);
        disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'GR ',num2str(GROUPELIST),' KruskalWallis.txt'])]);
    else
        try
            xlswrite(fullfile(dir1,[name,labelnode,'GR ',num2str(GROUPELIST),' KruskalWallis.xlsx']),infonew);
            disp(['Result .xlsx file saved: ' fullfile(dir1,[name,labelnode,'GR ',num2str(GROUPELIST),' KruskalWallis.xlsx'])]);
        catch
            writetxtfile(fullfile(dir1,[name,labelnode,'GR ',num2str(GROUPELIST),' KruskalWallis.txt']),infonew);
            disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'GR ',num2str(GROUPELIST),' KruskalWallis.txt'])]);
        end
    end
elseif isfield(job.c_statmatrix,'b_manova1_Mat')
    disp('not tested');
    return
    dir1 = job.e_statmatrixPath{1};
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    covariableall=[];
    covariablestring = job.c_statmatrix.b_manova1_Mat.b_Covariable_Mat;
    [token,remain] =strtok(covariablestring,',');
    covariableall =  [covariableall,{strtrim(token)}];
    while ~isempty(remain)
        [token,remain] =strtok(remain,',');
        covariableall =  [covariableall,{strtrim(token)}];
    end

    for icov = 1:numel(covariableall)
        notfoundstophere = 1;
        Pearsony = covariableall{icov};
        ycol = 0;
        for icol=1:size(info,2)
            if ~isnan(info{1,icol})
                if strcmp(strtrim(upper(deblank(info{1,icol}))), strtrim(upper(Pearsony)));
                    ycol = icol;
                end
            end
        end
        if ycol
            id= 1;
            for i=2:size(info,1)
                if isnumeric(info{i,ycol })
                    score(id,icov) = info{i,ycol };
                    id = id+1;
                else
                    score(id,icov) = nan;
                    id = id+1;
                end
            end

        else
            disp([Pearsony,' column not found, covariable unavailable'])
            out='Stat tests';
            return
        end
    end
    AllC = [];
    id =1;
    dir1 = job.e_statmatrixPath{1} ;
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    GROUPELIST = str2num(job.c_statmatrix.b_manova1_Mat.e_Anova1GR);
    %Use a specific groupe
    iduse = find(  sum(groupeall==GROUPELIST,2)& ~sum(isnan(score),2));

    for i= 1:size(MATall,2)
        for j = 1:size(MATall,3)
            if i~=j
                isgood = find(~isnan(MATall( iduse,i,j)));
                try
                    X = [(MATall( iduse(isgood),i,j)), score(iduse(isgood),:)];
                    [d,p,stats] = manova1(X, groupeall(iduse(isgood)));


                    panovaGroupe(i,j)= p;
                    panovaCoV(i,j)=1;% atab{3,end};
                    panovaGroupebyCoV(i,j)= 1;%atab{4,end};
                    nanova(i,j) = numel(iduse(isgood));
                    p_groupe = 1; %atab{2,6};
                    p_CoV= 1 ;%atab{3,6};
                    if i==1&j==2 %initialized once matrice results
                        [c,m,h,gnames] = multcompare(stats,'display','off');

                        diff = zeros(size(c,1),size(MATall,2),size(MATall,3) );
                        pdiff = zeros(size(c,1),size(MATall,2),size(MATall,3) );
                        for icomp=1:size(c,1)
                            labelmultcompare{icomp} = [gnames{c(icomp,1)}, '-' gnames{c(icomp,2)}];
                            labelmultcompareneg{icomp}  = [gnames{c(icomp,2)}, '-' gnames{c(icomp,1)}];
                        end
                    end
                    [c,m,h,gnames] = multcompare(stats,'display','off');
                    for icomp=1:size(c,1)
                        eval(['diff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=',num2str(c(icomp,4)),';']);
                        eval(['pdiff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=',num2str(c(icomp,6)),';']);
                    end
                    if p_groupe < alpha_threshold %p groupe
                        filename = fullfile(dir1,['STATS_Ancova ','GR ',num2str(GROUPELIST),' x groupe' ZoneList{i},' to ' ZoneList{i},' p= ',num2str(p_groupe),'.mat']);
                        save(filename,'stats');
                        % multcompare(stats)
                    end
                    if p_CoV < alpha_threshold %p groupe
                        filename = fullfile(dir1,['STATS_Ancova ','GR ',num2str(GROUPELIST),' x COV ' ZoneList{i},' to ' ZoneList{i},' p= ',num2str(p_CoV),'.mat']);
                        save(filename,'stats');
                        %  multcompare(stats)

                    end

                catch
                    panovaGroupe(i,j) = 1;
                    panovaCoV(i,j) = 1;
                    panovaGroupebyCoV(i,j) = 1;
                    nanova(i,j) = numel(iduse(isgood));
                    try
                        for icomp=1:size(c,1)
                            eval(['diff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=','nan',';']);
                            eval(['pdiff(',num2str(icomp),',',num2str(i),',',num2str(j),')','=','nan',';']);
                        end
                        disp(['Error stat for link ' num2str(i) 'to ',num2str(j) ])
                    catch
                        disp(['Error stat for link ' num2str(i) 'to ',num2str(j) ])
                    end
                end
            else
                panovaGroupe(i,j) = nan;
                panovaCoV(i,j) = nan;
                panovaGroupebyCoV(i,j) = nan;
            end

        end
    end
    try
        [FRD,Q] = mafdr(panovaGroupe(:));
        Q = reshape(Q,size(panovaGroupe));
        %SAVE FDR
        file = [name,labelnode,'ANCOVAFDR groupe p',num2str(alpha_threshold),'.mat'];
        meancorr = double( Q<alpha_threshold);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        [FRD,Q] = mafdr( panovaCoV(:));
        Q = reshape(Q,size( panovaCoV));
        %figure;hist(panovaCoV(:))
        %SAVE FDR
        file = [name,labelnode,'ANCOVAFDR COV p',num2str(alpha_threshold),'.mat'];
        meancorr = double( Q<alpha_threshold);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];

        [FRD,Q] = mafdr(panovaGroupebyCoV(:));
        Q = reshape(Q,size(panovaGroupebyCoV));
        %SAVE FDR
        file = [name,labelnode,'ANCOVAFDR groupe p',num2str(alpha_threshold),'.mat'];
        meancorr = double( Q<alpha_threshold);
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    catch
        disp('No FDR')
    end

    for igr = 1:numel(GROUPELIST)
        iduse = find(  sum(groupeall==GROUPELIST(igr),2));
        file = [name,labelnode,'meanGR',num2str(GROUPELIST(igr)),'.mat'];
        meanG1 = squeeze(nanmean(MATall(iduse,:,:) ,1));
        meancorr = meanG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
        infonew = [infonew;new];
    end
    for icomp = 1:size(diff,1)
        file = [name,labelnode,'Ancova diffGR',labelmultcompare{icomp},'.mat'];
        meanG1 = squeeze(diff(icomp,:,:)).*double(panovaGroupe<=alpha_threshold);
        meancorr = meanG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
        infonew = [infonew;new];
        file = [name,labelnode,'Ancova diffGR',labelmultcompareneg{icomp},'.mat'];
        meanG1 = squeeze(diff(icomp,:,:)).*double(panovaGroupe<=alpha_threshold);
        meancorr = -meanG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
        infonew = [infonew;new];

        file = [name,labelnode,'Ancova diffGR',labelmultcompare{icomp},'p',num2str(alpha_threshold),'.mat'];
        meanG1 = squeeze(diff(icomp,:,:)).*double(panovaGroupe<alpha_threshold).*double(squeeze(pdiff(icomp,:,:)) <alpha_threshold);
        meancorr = meanG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
        infonew = [infonew;new];
        file = [name,labelnode,'Ancova diffGR',labelmultcompareneg{icomp},'p',num2str(alpha_threshold),'.mat'];
        meanG1 = squeeze(diff(icomp,:,:)).*double(panovaGroupe<alpha_threshold).*double(squeeze(pdiff(icomp,:,:)) <alpha_threshold);
        meancorr = -meanG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
        infonew = [infonew;new];
    end
    file = [name,labelnode,'Ancova GR ',num2str(GROUPELIST),', pGroupe',num2str(alpha_threshold),'.mat'];
    meancorr = double(panovaGroupe<alpha_threshold);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];

    file = [name,labelnode,'Ancova GR ',num2str(GROUPELIST), ', pCoV',num2str(alpha_threshold),'.mat'];
    meancorr = double(panovaCoV<alpha_threshold);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];

    file = [name,labelnode,'Ancova GR ',num2str(GROUPELIST), ', pGroupebyCoV',num2str(alpha_threshold),'.mat'];
    meancorr = double(panovaGroupebyCoV<alpha_threshold);
    save(fullfile(dir1,[file]),'ZoneList','meancorr');
    disp(['Save: ', fullfile(dir1,[file])]);
    new = [{dir1},{file}, {ZONEid},{0} ];
    infonew = [infonew;new];
    if ~strcmp(fullfile(info{isubject,1}, ZONEid),fullfile(dir1,  ZONEid))
        copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid))
    end
    if ismac
        % Code to run on Mac platform problem with xlswrite
        writetxtfile(fullfile(dir1,[name,labelnode,'Ancova', 'GR ',num2str(GROUPELIST), '.txt']),infonew);
        disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'Ancova', 'GR ',num2str(GROUPELIST), '.txt'])]);
    else
        try
            xlswrite(fullfile(dir1,[name,labelnode,'Ancova', ' GR ',num2str(GROUPELIST), '.xlsx']),infonew);
            disp(['Result .xlsx file saved: ' fullfile(dir1,[name,labelnode,'Ancova', ' GR ',num2str(GROUPELIST), '.xlsx'])]);
        catch
            writetxtfile(fullfile(dir1,[name,labelnode,'Ancova', ' GR ',num2str(GROUPELIST), '.txt']),infonew);
            disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'Ancova', ' GR ',num2str(GROUPELIST), '.txt'])]);
        end
    end
elseif isfield(job.c_statmatrix,'b_fitMANCOVAN_Mat')
    dir1 = job.e_statmatrixPath{1};
    infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
    covariableall=[];
    covariablestring = job.c_statmatrix.b_fitMANCOVAN_Mat.b_fitMANCOVAN_Covariable;
    [token,remain] =strtok(covariablestring,',');
    covariableall =  [covariableall,{strtrim(token)}];
    while ~isempty(remain)
        [token,remain] =strtok(remain,',');
        covariableall =  [covariableall,{strtrim(token)}];
    end
    for icov = 1:numel(covariableall)
        Pearsony = covariableall{icov};
        ycol = 0;
        for icol=1:size(info,2)
            if ~isnan(info{1,icol})
                if strcmp(strtrim(upper(deblank(info{1,icol}))), strtrim(upper(Pearsony)));
                    ycol = icol;
                end
            end
        end
        if ycol
            id= 1;
            for i=2:size(info,1)
                try
                    score(id,icov) = info{i,ycol };
                    id = id+1;
                catch
                    score(id,icov) = nan;
                    id = id+1;
                end
            end

        else
            disp([Pearsony,' column not found, no regression could be compute'])
            out='Stat tests';
            return
        end
    end
    GROUPELIST = str2num(job.c_statmatrix.b_fitMANCOVAN_Mat.e_fitMANCOVAN_GR);
    for igr = 1:numel(GROUPELIST)
        iduse = find(  sum(groupeall==GROUPELIST(igr),2));
        file = [name,labelnode,'meanGR',num2str(igr),'.mat'];
        meanG1 = squeeze(nanmean(MATall(iduse,:,:) ,1));
        meancorr = meanG1;
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{GROUPELIST(igr)} ];
        infonew = [infonew;new];
    end


    % if isfield(job.c_statmatrix.b_fitMANCOVAN_Mat.c_statpermutation,'b_Nopermutation')

    LabelSTAT = ['Groupe', covariableall];
    for i= 1:size(MATall,2)
        for j = 1:size(MATall,3)
            if i~=j
                try
                    iduse = find( sum(groupeall==GROUPELIST,2));
                    Yval = MATall( iduse,i,j);
                    GR = groupeall(iduse);
                    CovariateX = [];
                    for idvariable = 1:numel(covariableall)
                        X = score(iduse,idvariable);
                        CovariateX =  [CovariateX, X];
                    end

                    idok = find(~(isnan(Yval)|sum(isnan(CovariateX),2)));
                    %F MANCOVAN Same results as spss
                    [ T, p,  Fanova, Panova, stats ] = mancovan( Yval( idok), GR( idok),CovariateX( idok,:));
                    Fanoval(:,i,j) = Fanova;
                    Panoval(:,i,j) = Panova;
                    Banoval(:,i,j) = stats.B(2:end);

                catch
                    disp('Error MANCOVAN')
                    disp('MANCOVAN function must be install in the set path see:  William Gruner (2022). MANCOVAN (https://www.mathworks.com/matlabcentral/fileexchange/27014-mancovan), MATLAB Central File Exchange. Retrieved October 25, 2022. ')
                    p = mfilename('fullpath');
                    [pathname,filename,ext]=fileparts(p);
                    disp(['Add the fonction in matlab path: addpath(''',pathname, '\External\mancovan_496'')'])
                    disp('And verify if Fisher transform is not apply twice, data must be real')
                    %return
                end
            end
        end
    end


    for icond = 1:numel(LabelSTAT)
        file = [name,labelnode,' fanova', LabelSTAT{icond},'.mat'];
        meancorr =squeeze(Fanoval(icond,:,:));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    end
    for icond = 1:numel(LabelSTAT)
        file = [name,labelnode,'panova 1-p', LabelSTAT{icond},'.mat'];
        meancorr =  squeeze(1-Panoval(icond,:,:));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    end
    for icond = 1:numel(LabelSTAT)
        file = [name,labelnode,'panova', LabelSTAT{icond},'.mat'];
        meancorr = squeeze(Panoval(icond,:,:));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    end
    for icond  = 1:numel(LabelSTAT)
        file = [name,labelnode,'Beta', LabelSTAT{icond },'.mat'];
        meancorr = squeeze(Banoval(icond ,:,:));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    end
    for icond = 1:numel(LabelSTAT)
        file = [name,labelnode,'Beta', LabelSTAT{icond },'p',num2str(alpha_threshold),'.mat'];
        meancorr =  squeeze(Banoval(icond ,:,:).*double(Panoval(icond ,:,:)<alpha_threshold));
        save(fullfile(dir1,[file]),'ZoneList','meancorr');
        disp(['Save: ', fullfile(dir1,[file])]);
        new = [{dir1},{file}, {ZONEid},{0} ];
        infonew = [infonew;new];
    end

    if isfield(job.c_statmatrix.b_fitMANCOVAN_Mat.c_statpermutation,'b_permutation')
        LabelSTAT = ['Groupe', covariableall];
        halfMAT = MATall(:,idhalf);
        nperm= str2num(job.c_statmatrix.b_fitMANCOVAN_Mat.c_statpermutation.b_permutation.e_npermutation);
        f_dist = nan(numel(LabelSTAT), nperm,numel(idhalf));
        Fanoval = nan(numel(LabelSTAT),numel(idhalf));
        Panoval= nan(numel(LabelSTAT),numel(idhalf));
        Banoval= nan(numel(LabelSTAT),numel(idhalf));
        ppermmax= str2num(job.c_statmatrix.b_fitMANCOVAN_Mat.c_statpermutation.b_permutation.e_npermutation);
        %% SECTION TITLE
        % DESCRIPTIVE TEXT
        pperm = nan(numel(LabelSTAT),numel(idhalf));
        iduse = find( sum(groupeall==GROUPELIST,2));
        for iperm=1:nperm
            permfix(iperm,:) =  randperm(numel(iduse));
        end
        disp(['Running ', num2str(nperm),' permutations on ', num2str(size(halfMAT,2)),' link :'])


        for ilink = 1:size(halfMAT,2)
            fprintf('%s ',num2str(ilink));
            if mod(ilink,30)==0
                fprintf('\n ');
            end
            try
                Yval = halfMAT(iduse,ilink);
                GR = groupeall(iduse);
                CovariateX = [];
                for idvariable = 1:numel(covariableall)
                    X = score(iduse,idvariable);
                    CovariateX =  [CovariateX, X];
                end
                idok = find(~(isnan(Yval)|sum(isnan(CovariateX),2)));
                %F MANCOVAN Same results as spss
                try
                    [ T, p,  Fanova, Panova, stats ] = mancovan( Yval( idok), GR( idok),CovariateX( idok,:));
                catch
                    disp('Error MANCOVAN')
                    disp('MANCOVAN function must be install in the set path see:  William Gruner (2022). MANCOVAN (https://www.mathworks.com/matlabcentral/fileexchange/27014-mancovan), MATLAB Central File Exchange. Retrieved October 25, 2022. ')
                    p = mfilename('fullpath');
                    [pathname,filename,ext]=fileparts(p);
                    disp(['Add the fonction in matlab path: addpath(''',pathname, '\External\mancovan_496'')'])
                    disp('And verify if Fisher transform is not apply twice, data must be real')
                    Fanova = 0;
                    Panova =1;
                end
                Fanoval(:,ilink) = Fanova;
                Panoval(:,ilink) = Panova;
                try
                    Banoval(:,ilink) =   stats.B(2:end);
                catch
                    Banoval(:,ilink)=nan;
                end
                for idperm=1:nperm
                    Yval = halfMAT( iduse,ilink);
                    GR = groupeall(iduse(permfix(idperm,:)));
                    CovariateX = [];
                    for idvariable = 1:numel(covariableall)
                        X = score(iduse(permfix(idperm,:)),idvariable);
                        CovariateX =  [CovariateX, X];
                    end
                    idok = find(~(isnan(Yval)|sum(isnan(CovariateX),2)));
                    [ T, p,  Fanova, Panova, stats ] = mancovan( Yval( idok), GR( idok),CovariateX( idok,:));
                    f_dist(:, idperm,ilink) = Fanova(:,1);
                end
            catch
                fprintf('error link ');
            end

        end

        f_dist_all =  f_dist ;
        f_dist_all(find(isnan(f_dist)))=[];
        [fecdf,xecdf] = ecdf( f_dist_all);
        clustercritical=xecdf(sum(fecdf<(1-alpha_threshold)));
        %         figure;
        %         plot(xecdf,fecdf)
        %         xlabel('xecfg=valeur de F')
        %         ylabel('fecdf=F(x) cumulative distribution function')
        if isfield(job.c_statmatrix.b_fitMANCOVAN_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting,'b_MCT_ClusterBased')
            for icomp = 1:size(Fanoval,1)
                statrand = permute(f_dist(icomp,:,:), [3 2 1]);
                statobs = permute(Fanoval(icomp,:),[2,1]);
                tmp = randperm(numel(statrand));
                statrand =reshape(statrand(tmp),size(statrand,1),size(statrand,2));
                %                 neighbourdist = 3;
                %                 disp(['critical for clustering p=', num2str(alpha_threshold),' F=',num2str(clustercritical)])
                %                 stat3 = FindClusterBasedPermutationInMatrix(chanpos, neighbourdist,clustercritical, statobs, statrand );
                %
                neighbourdist = job.c_statmatrix.b_fitMANCOVAN_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_neighbourdist;
                minnbchan = job.c_statmatrix.b_fitMANCOVAN_Mat.c_statpermutation.b_permutation.c_statMultipleComparaisonTesting.b_MCT_ClusterBased.e_minnbchan;

                disp(['critical for clustering p=', num2str(alpha_threshold),' F=',num2str(clustercritical)])
                [stat, matneig] = FindClusterBasedPermutationInMatrix(chanpos, neighbourdist,clustercritical, statobs, statrand,minnbchan  );
                save(fullfile(dir1,[LabelSTAT{icomp},'stat.mat']),'stat')
                disp(['Save cluster', fullfile(dir1,[LabelSTAT{icomp},'stat.mat']) ])

                if isfield(stat,'posclusterslabelmat')
                    for iposcluster = 1:numel(stat.posclusters)
                        file = [name,labelnode, LabelSTAT{icomp},num2str(nperm),' posCluster', num2str(iposcluster),' MCp=',num2str(stat.posclusters(iposcluster).prob),'.mat'];
                        tmp=(stat.posclusterslabelmat==iposcluster);
                        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                        save(fullfile(dir1,[file]),'ZoneList','meancorr');
                        disp(['Save: ', fullfile(dir1,[file])]);
                        new = [{dir1},{file}, {ZONEid},{0} ];
                        infonew = [infonew;new];
                    end
                else
                    disp('No positive cluster found')
                end
                if isfield(stat,'negclusterslabelmat') %technicaly no neg cluster in f value...
                    for inegcluster = 1:numel(stat.negclusters)
                        file = [name,labelnode, LabelSTAT{icomp},num2str(nperm),' negCluster', num2str(iposcluster),' MCp=',num2str(stat.negclusters(inegcluster).prob),'.mat'];
                        tmp=(stat.negclusterslabelmat==inegcluster);
                        meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
                        save(fullfile(dir1,[file]),'ZoneList','meancorr');
                        disp(['Save: ', fullfile(dir1,[file])]);
                        new = [{dir1},{file}, {ZONEid},{0} ];
                        infonew = [infonew;new];
                    end
                else
                    disp('No negative cluster found')
                end
            end
        end


        %check whole distribution
        for icomp = 1:size(Fanoval,1)
            for ilink = 1:size(halfMAT,2)
                if sum(xecdf<=Fanoval(icomp,ilink))==0
                    pperm(icomp,ilink )= nan;
                else
                    pperm(icomp,ilink ) = 1-fecdf(sum(xecdf<=Fanoval(icomp,ilink)));
                end
            end
        end
        disp( [num2str(nperm) ' Permutation p<',num2str(alpha_threshold), ' give a marginal F value=',num2str(xecdf(sum(fecdf<(1-alpha_threshold))))])

        if 1 % check maxdistribution
            for  icond =1:size(Fanoval,1)
                dist_sup = squeeze(max(squeeze(f_dist(icond,:,:))'));
                dist_sup = squeeze(max(squeeze(f_dist(icond,:,:))'));
                [fecdf,xecdf] = ecdf( [dist_sup]);
                disp( [ LabelSTAT{icond}, ' ' ,num2str(nperm) ' Max Permutation over all link p<',num2str(alpha_threshold), ' give a marginal F value=',num2str(xecdf(sum(fecdf<=(1-alpha_threshold))))])
            end
        end


        for icond = 1:numel(LabelSTAT)
            file = [name,labelnode,' fanova', num2str(icond), LabelSTAT{icond},'.mat'];
            tmp=Fanoval(icond,:);
            meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{0} ];
            infonew = [infonew;new];
        end

        for icond = 1:numel(LabelSTAT)
            file = [name,labelnode,'panova', ' 1-p',LabelSTAT{icond},'.mat'];
            tmp=1-Panoval(icond,:);
            meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{0} ];
            infonew = [infonew;new];
        end

        for icond= 1:numel(LabelSTAT)
            file = [name,labelnode,'panova',num2str(nperm),'perm',' 1-p', LabelSTAT{icond},'.mat'];
            tmp=squeeze(1-pperm(icond,:));
            meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{0} ];
            infonew = [infonew;new];
        end


        for icond= 1:numel(LabelSTAT)
            file = [name,labelnode,'panova', LabelSTAT{icond},'.mat'];
            tmp = squeeze(1-Panoval(icond,:))
            meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{0} ];
            infonew = [infonew;new];
        end

        for icond= 1:numel(LabelSTAT)
            file = [name,labelnode,'panova',num2str(nperm),'perm', LabelSTAT{icond},'.mat'];
            tmp=squeeze(pperm( icond,:));
            meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{0} ];
            infonew = [infonew;new];
        end

        for icond= 1:numel(LabelSTAT)
            file = [name,labelnode,'beta',LabelSTAT{ icond},'.mat'];
            tmp = squeeze(Banoval( icond,:));
            meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{0} ];
            infonew = [infonew;new];
        end

        for icond= 1:numel(LabelSTAT)
            file = [name,labelnode,'beta',LabelSTAT{ icond},'pperm',num2str(alpha_threshold),'.mat'];
            tmp = squeeze(Banoval( icond,:).*double(pperm( icond,:)<alpha_threshold));
            meancorr = halfmat2mat(tmp,idhalf, size(MATall,2));
            save(fullfile(dir1,[file]),'ZoneList','meancorr');
            disp(['Save: ', fullfile(dir1,[file])]);
            new = [{dir1},{file}, {ZONEid},{0} ];
            infonew = [infonew;new];
        end



    end
    if ~strcmp(fullfile(info{isubject,1}, ZONEid),fullfile(dir1,  ZONEid))
        copyfile(fullfile(info{isubject,1}, ZONEid),  fullfile(dir1,  ZONEid))
    end
    if ismac
        % Code to run on Mac platform problem with xlswrite
        writetxtfile(fullfile(dir1,[name,labelnode,'MANCOVAN.txt']),infonew);
        disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'MANCOVAN.txt'])]);
    else
        try
            xlswrite(fullfile(dir1,[name,labelnode,'MANCOVAN.xlsx']),infonew);
            disp(['Result .xlsx file saved: ' fullfile(dir1,[name,labelnode,'MANCOVAN.xlsx'])]);
        catch
            writetxtfile(fullfile(dir1,[name,labelnode,'MANCOVAN.txt']),infonew);
            disp(['Result .txt file saved: ' fullfile(dir1,[name,labelnode,'MANCOVAN.txt'])]);
        end
    end
end
disp('Visualize the results using: GUI_LookMatrices')
out='Stat tests';

matlabbatch{1}.spm.tools.nirsHSJ.M_Connectivity.E_statmatrix = job;
save(fullfile(dir1,'Batch_stat.mat'),'matlabbatch')

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

function mat = halfmat2mat(halfmat,idhalf, n)
%fonction pour remettre la demi matrice en colonne en matrice carr�
%sym�trique en faisant attention de ne pas doubl� la diagonale.
mat = zeros(n,n);
mat(idhalf)=squeeze(halfmat);
mat = mat+flipud(rot90(mat));
mat =  mat.*(ones(n,n)-eye(n)*0.5 );  %ajustement pour ne pas doubl� la diagonale dans l'addition .
