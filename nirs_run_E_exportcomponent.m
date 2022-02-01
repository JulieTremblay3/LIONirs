function out = nirs_run_E_exportcomponent(job)
%Look in the list of selected factor to export the value
  
     [pathstr, name, ext] = fileparts(job.NIRSmat{1});
    load(job.NIRSmat{1});
    load(fullfile(pathstr,'SelectedFactors.mat'));
    label = job.i_exportcomponent_label;
    switch job.m_exportcomponent_type
        case 0 
            type='GLM';
        case 1 
            type='PARAFAC';
        case 2 
            type='PCA';
    end
    if isfield(job.NewDirCopyNIRSTRUE,'CreateNIRSCopy')
        NewNIRSdir = job.NewDirCopyNIRSTRUE.CreateNIRSCopy.NewNIRSdir;
        disp(['Create directory for condition ',NewNIRSdir])
        NewDirCopyNIRS = 1;
        pathout = [pathstr,filesep, NewNIRSdir];
        if ~isdir(pathout)
            mkdir(pathout);
        end
    else
        NewDirCopyNIRS = 0;
    end
    all = [];
    alllabel = [];
   
    NC = NIRS.Cf.H.C.N;

       findch = 0;
       if findch
            curveall = [];
            onset = [4.966413
            54.937741
            104.909069
            154.880396
            204.851724
            254.823052
            304.79438
            354.765708
            404.737036
            454.708364
            504.679692
            554.65102
            604.622348
            654.593676
            704.565004
            754.536332]
        fs = NIRS.Cf.dev.fs
        ALLPARAFAC = zeros(1010,17);
        idonset = 1;
        tref =   (1/fs:1/fs:1/fs*1010)-5
       end
    for icomp = 1:numel(PARCOMP)
         if ~isempty(findstr(PARCOMP(icomp).label,label))
            if strcmp(PARCOMP(icomp).type,type)
             listgood=PARCOMP(icomp).listgood;
             topo=PARCOMP(icomp).topo;
             tmp = NaN(NC,1);
              tmp(listgood)=topo;
              all = [all,tmp ];
              alllabel = [alllabel,{PARCOMP(icomp).label}];
              if findch
              id = find(listgood==31); %CH 31
              if ~isempty(id)
                  
            fs = NIRS.Cf.dev.fs
                 indt = PARCOMP(1)
                 t = PARCOMP(icomp).indt*1/fs;
                 
                 
                  idstart = find( tref <= (PARCOMP(icomp).indt(1)*1/fs - onset(1)))
                 
                  curveall = PARCOMP(icomp).Xm(:,id,1)*2
                  nb = numel(curveall)-1
                  ALLPARAFAC(idstart:idstart+nb,idonset)= curveall 
                  
                  idonset = idonset +1
              end
              end
            end
        end
    end
    if findch
        figure
        plot(tref,mean(ALLPARAFAC,2)) 
        figure
        plot(tref,ALLPARAFAC)
    end
    idkeep = [1:size(all,2)];
    Amean = nanmean(all(:,idkeep),2);
    Astd = nanstd(all(:,idkeep),1,2);
    Anb  = sum(~isnan(all(:,idkeep)),2)-1;
    tmp = [alllabel;num2cell(all)];
    
    all(:,idkeep)
    if ismac
        writetxt_asxlswrite(fullfile(pathout,['Component' type,label,'.xls']),tmp);
        disp(['Create: ', fullfile(pathout,['Component' type,label,'.xls']) ]);
    else
        try
            xlswrite(fullfile(pathout,['Component' type,label,'.xls']),tmp);
            disp(['Create: ', fullfile(pathout,['Component' type,label,'.xls']) ]);
        catch
             writetxtfile(fullfile(pathout,['Component' type,label,'.xls'],tmp);
             disp(['Create: ', fullfile(pathout,['Component' type,label,'.xls']) ]);
        end
    end
     
    A = Amean(1:end/2);
    save(fullfile(pathout,['D1matrix ','mean', label,type,'.mat']),'A' );
    A = Astd(1:end/2);
    save(fullfile(pathout,['D1matrix ','std', label,type,'.mat']),'A' );
    for ich = 1:size(A,1)
        tmp = all(ich ,idkeep);
        idbad = find(isnan(tmp)|isinf(tmp));
        if isempty(idbad)
            tmp(idbad)=[];
        end
        [h,p,ci,stats] = ttest(tmp);
        T(ich) = stats.tstat;
    end
    idbad = find(isnan(T)|isinf(T));
    T(idbad)=0;
    A = T;
    save(fullfile(pathout,['D1matrix ','T', label,type,'.mat']),'A' );
    figure;plot(T,'x')
    
%     A = A(1:end/2)
%     save(fullfile(pathout,['D1matrix ','T', label,type,'.mat']),'A' );

    out.NIRSmat = job.NIRSmat;