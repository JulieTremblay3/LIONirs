function out = nirs_run_E_extractcomponent(job)
% Batch detection of the components
if isfield(job.c_extractcomponent,'b_extractcomponent_PCA')
    NIRSmatlst = job.c_extractcomponent.b_extractcomponent_PCA.NIRSmat;
    try 
         pca(rand(10)); %is pca function is install ?
    catch
         disp('Uncomplete Extract PCA, Please Install Matlab Statistics and Machine Learning Toolbox')
         out.NIRSmat = job.c_extractcomponent.b_extractcomponent_PCA.NIRSmat;
         return
    end
    i_extract_pourcentagech= job.c_extractcomponent.b_extractcomponent_PCA.i_extract_pourcentagech;
    m_extractcomponentfigure= job.c_extractcomponent.b_extractcomponent_PCA.m_extractcomponentfigure;
    m_extractcomponent= job.c_extractcomponent.b_extractcomponent_PCA.m_extractcomponent;
    for filenb=1:size(NIRSmatlst,1) %Loop over all subjects
        %Load NIRS.mat information
        tic
        NIRS = [];
        NIRSmat= NIRSmatlst{filenb,1};
        load(NIRSmat);
        
        
        [dirout,~,~] = fileparts(NIRSmat);
        lst = length(NIRS.Dt.fir.pp);
        rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
        NC = NIRS.Cf.H.C.N;
        fs = NIRS.Cf.dev.fs;
        indconsecutifthreshold = 1;     %en dur�e en sample
        nbchminimum = i_extract_pourcentagech/100; %0.05;             %en pourcentage
        fprintf('%s\n','File processed');
        %Open file 1 and check where are the noise marked.
        for f = 1:size(rDtp,1)
            
            d = fopen_NIR(rDtp{f,1},NC)';
            [dir1,fil1,ext1] = fileparts(rDtp{f});
            vmrk_path = fullfile(dir1,[fil1,'.vmrk']);
            %         handles.file_vmrk = handles.NIRS.Dt.fir.pp(end).p{idfile}; %used
            %         in save noise
            mrk_type_arr = cellstr('bad_step');
            mrks = [];
            ind = [];
            noise =  logical(zeros(size(d)));
            time  = 1/fs:1/fs:1/fs*size(d,1);
            [ind_dur_ch] = read_vmrk_find(vmrk_path,mrk_type_arr);
            if ~isempty(ind_dur_ch)
                maxpoint  = ind_dur_ch(:,1)+ind_dur_ch(:,2);
                badind = find(maxpoint>size(noise,1));
                if ~isempty(badind)
                    disp(['Warning file ' vmrk_path ' marker : ' num2str(badind') ' are out of range in the data file'])
                    ind_dur_ch(badind,2)=size(noise,2)- ind_dur_ch(badind,1);
                end
                for Idx = 1:size(noise,2)
                    mrks = find(ind_dur_ch(:,3)==Idx|ind_dur_ch(:,3)==0);
                    ind = ind_dur_ch(mrks,1);
                    indf = ind + ind_dur_ch(mrks,2) - 1;
                    if ~isempty(ind)
                        try
                            for i = 1:numel(ind)
                                noise(ind(i):indf(i),Idx) = 1;
                            end
                        catch
                           disp('Noise reading problem')
                        end
                    end
                end
            end
            
            
            % ici group channel with the same noise latency
            figureon = m_extractcomponentfigure;
            ind = find((sum(noise,2)./size(noise,2))>nbchminimum );
            ind = find((sum(noise,2)./size(noise,2)));
            inddiff = diff(ind);
            if isempty(ind)
                disp(['No noisy event found then no component are extracted in file ', fil1])
                break
            end
            idsep = find(inddiff>indconsecutifthreshold);
            if isempty(idsep)
                idstart  =[ind(1)];
                idstop = [ind(end)];
            else
                idstart  =[ind(1);ind(idsep(1:end)+1)];
                idstop =  [ind(idsep(1:end)-1);ind(end)];
            end
            eventbadstartstop = [idstart,idstop] ;
            
            for ievent=1:size(eventbadstartstop,1)
                
                %listchannel to find component
                indt = [eventbadstartstop(ievent,1):eventbadstartstop(ievent,2)];
                if numel(indt)>1
                    if 0 %figureon
                        figure;subplot(3,2,[1,2]); plot(sum(noise,2)/size(noise,2));hold on
                        for i=1:size(eventbadstartstop,1)
                            plot([eventbadstartstop(i,1),eventbadstartstop(i,1)],[0, 1],'b')
                            plot([eventbadstartstop(i,2),eventbadstartstop(i,2)],[0, 1],'r')
                        end
                        plot(indt(round(numel(indt)/2)),1,'x');
                        title(fil1)
                    end
                    
                    plotLst= find(sum(noise(indt,:),1)>0);
                    listok=  find(NIRS.Cf.H.C.ok(plotLst,f)); %bon canaux ?
                    plotLst = plotLst(listok);
                    intensnorm = d(indt,:);
                    % take both wavelengt channel
                    MeasListActplotLst = zeros(size(d,2)./2,1); %combine both lambda
                    id830 = find((plotLst-numel(MeasListActplotLst)) <=0);
                    if ~isempty(id830)
                        MeasListActplotLst(plotLst(id830 ),:)=1;
                    end
                    id690 =  find((plotLst-numel(MeasListActplotLst)) >0);
                    if ~isempty(id690)
                        MeasListActplotLst(plotLst(id690)-numel(MeasListActplotLst),:)=1;
                    end
                    listgood = find(MeasListActplotLst);
                    listgoodPCA = [listgood;listgood+NC/2];
                    %Detrent DATA segment for centrering
                    X = 1:1:size(intensnorm,1);
                    Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
                    Mb2 =  intensnorm(1,:)'; %offset
                    A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
                    spar = intensnorm - A;
                    % figure;plot(intensnorm )
                    %    figure;plot( spar)
                    spar = cat(3,spar(:,1:end/2),spar(:,end/2+1:end)); 
                    spartmp = spar(:,listgood,:);
                    listgood = find(MeasListActplotLst);
                    %PCA first component !
                    if 1
                        idwhile = 1;
                        nbPCA= 2; %maximal number of iteration
                        while idwhile==1  %find(std(data(:,listgood,1)./mean(data(:,listgood,1))>0.1))
                            
                            d1 = cat(2,spartmp(:,:,1),spartmp(:,:,2));
                            c = d1'*d1;
                            
                            [coeff,score,latent,tsquared,explained] = pca(c);
                            
                             
                            iscumlowerthantop = cumsum(explained)<job.c_extractcomponent.b_extractcomponent_PCA.i_extractnoiseupto_nbPCA;
                            ishigherthanmin = explained>job.c_extractcomponent.b_extractcomponent_PCA.i_extractnoise_nbPCA;
                            atleastone = zeros(size(explained));
                            atleastone(1) = 1;
                            lstSV = find((iscumlowerthantop&ishigherthanmin)|atleastone);
                            labelexplainvar =     sprintf('%02.0fn=%d',sum(explained(lstSV)),numel(lstSV));
                            if isempty(lstSV)
                                break
                            end
                            [v,s,foo]=svd(c);
                            svs = diag(s);
                            u = d1*v*inv(s);
                            Xm =  u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)';
                            
                            try
                                load(fullfile(dirout,'SelectedFactors.mat'))
                                newfile = 0;
                            catch
                                %donot exist create the stucture
                                PARCOMP.file= f;
                                PARCOMP.filestr =  fil1;
                                PARCOMP.module  =numel(NIRS.Dt.fir.pp);
                                PARCOMP.modulestr = NIRS.Dt.fir.pp(end).pre;
                                PARCOMP.listgood =  listgoodPCA;
                                PARCOMP.indt = indt; %indice de temps.
                                PARCOMP.data =d1;
                                PARCOMP.u =  u;
                                PARCOMP.s = s;
                                PARCOMP.v = v;
                                PARCOMP.ComponentToKeep =1;
                                lstSV =1;
                                PARCOMP.Xm = PARCOMP.u(:,lstSV)*PARCOMP.s(lstSV,lstSV)*PARCOMP.v(:,lstSV)';
                                labelid  =  job.c_extractcomponent.b_extractcomponent_PCA.m_extractcomponent;
                                PARCOMP.label= [labelid,labelexplainvar,' ' ,sprintf('%03.0f',size(PARCOMP,2)),' ',fil1];
                                PARCOMP.type = 'PCA';
                                PARCOMP.topo = s(lstSV,lstSV)*v(:,lstSV)';
                                newfile = 1;
                            end
                            if newfile == 0
                                id = numel(PARCOMP);
                                PARCOMP(id+1).file= f;
                                PARCOMP(id+1).filestr =  fil1;
                                PARCOMP(id+1).module  =numel(NIRS.Dt.fir.pp);
                                PARCOMP(id+1).modulestr =  NIRS.Dt.fir.pp(end).pre;
                                PARCOMP(id+1).listgood =  listgoodPCA;
                                PARCOMP(id+1).indt = indt; %indice de temps.
                                PARCOMP(id+1).data = d1;
                                PARCOMP(id+1).u =  u;
                                PARCOMP(id+1).s = s;
                                PARCOMP(id+1).v = v;
                                PARCOMP(id+1).ComponentToKeep =1;
                                lstSV =1;
                                PARCOMP(id+1).Xm = PARCOMP(id+1).u(:,lstSV)*PARCOMP(id+1).s(lstSV,lstSV)*PARCOMP(id+1).v(:,lstSV)';
                                labelid  = job.c_extractcomponent.b_extractcomponent_PCA.m_extractcomponent;;
                                PARCOMP(id+1).label= [labelid,labelexplainvar,' ' , sprintf('%03.0f',size(PARCOMP,2)),' ',fil1];
                                PARCOMP(id+1).type = 'PCA';
                                PARCOMP(id+1).topo = s(lstSV,lstSV)*v(:,lstSV)';
                                %  plot(PARCOMP(id+1).u(:,lstSV),'r');
                                %  title('RED pca, Blue PARAFAC');
                            end
                            %figure;plot(d(indt,listgoodPCA));
                            idwhile=0;  %no loop
                            
                            if 0
                                if nbPCA==1 %reference zscore on the first loop
                                    zmean = mean(d1(:));
                                    zstd = std(d1(:));
                                end
                                tmp = (d1(:)-zmean)./zstd;
                                
                                spartmp = spartmp - reshape(Xm,size(Xm,1),size(Xm,2)/2,2);
                                
                                zscorecorr = (spartmp(:)-zmean)/zstd;
                                if figureon
                                    figure;
                                    subplot(3,2,1)
                                    plot(time(indt), reshape(d1,size(d1)));
                                    title([num2str(indt(1)),'to',num2str(indt(end))])
                                    
                                    subplot(3,2,2)
                                    plot(time(indt),reshape(tmp,size(d1)))
                                    title(['zscore > 3'] )
                                    subplot(3,2,3);plot(time(indt),spartmp(:,:,1))
                                    title(['FIRST ' num2str(explained(1)),'SECOND ',num2str(explained(2))])
                                    subplot(3,2,4);plot(time(indt),reshape(zscorecorr,size(d1)))
                                end
                                if sum(abs(zscorecorr(:))>3)>10 %if still a lot a variation go up too 3 decompositon
                                    idwhile = 1; %continue add other parafac not clean yet
                                    spartmp = spartmp - reshape(Xm,size(Xm,1),size(Xm,2)/2,2);
                                    if figureon
                                        subplot(3,2,5);plot(time(indt),spartmp(:,:,1))
                                        subplot(3,2,6);plot(time(indt),Xm)
                                        title('loop PCA')
                                    end
                                else
                                    idwhile = 0; %stop loop
                                    if figureon
                                        title('STOP HERE')
                                    end
                                end
                                nbPCA =nbPCA + 1;
                                if nbPCA >1
                                    idwhile=0;  %stop loop
                                end
                            end
                            %FIND the component more equal and less
                            save(fullfile(dirout,'SelectedFactors.mat'),'PARCOMP');
                        end
                    end
                end
            end
            
        end
    end
elseif isfield(job.c_extractcomponent,'b_extractnoise_PARAFAC')
    NIRSmatlst = job.c_extractcomponent.b_extractnoise_PARAFAC.NIRSmat;
    
    i_extract_pourcentagech= job.c_extractcomponent.b_extractnoise_PARAFAC.i_extract_pourcentagech;
    m_extractcomponentfigure= job.c_extractcomponent.b_extractnoise_PARAFAC.m_extractcomponentfigure;
    for filenb=1:size(NIRSmatlst,1) %Loop over all subjects
        %Load NIRS.mat information
        tic
        NIRS = [];
        NIRSmat= NIRSmatlst{filenb,1};
        load(NIRSmat);
        [dirout,~,~] = fileparts(NIRSmat);
        lst = length(NIRS.Dt.fir.pp);
        rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
        NC = NIRS.Cf.H.C.N;
        fs = NIRS.Cf.dev.fs;
        indconsecutifthreshold = 1;     %en dur�e en sample
        nbchminimum = i_extract_pourcentagech/100; %0.05;             %en pourcentage
        fprintf('%s\n','File processed');
        %Open file 1 and check where are the noise marked.
        for f = 1:size(rDtp,1)
            
            d = fopen_NIR(rDtp{f,1},NC)';
            [dir1,fil1,ext1] = fileparts(rDtp{f});
            vmrk_path = fullfile(dir1,[fil1,'.vmrk']);
            %         handles.file_vmrk = handles.NIRS.Dt.fir.pp(end).p{idfile}; %used
            %         in save noise
            mrk_type_arr = cellstr('bad_step');
            mrks = [];
            ind = [];
            noise =  logical(zeros(size(d)));
            time  = 1/fs:1/fs:1/fs*size(d,1);
            [ind_dur_ch] = read_vmrk_find(vmrk_path,mrk_type_arr);
            if ~isempty(ind_dur_ch)
                maxpoint  = ind_dur_ch(:,1)+ind_dur_ch(:,2);
                badind = find(maxpoint>size(noise,1));
                if ~isempty(badind)
                    disp(['Warning file ' vmrk_path ' marker : ' num2str(badind') ' are out of range in the data file'])
                    ind_dur_ch(badind,2)=size(noise,1)- ind_dur_ch(badind,1);
                end
                for Idx = 1:size(noise,2)
                    mrks = find(ind_dur_ch(:,3)==Idx |ind_dur_ch(:,3)==0 );
                    ind = ind_dur_ch(mrks,1);
                    indf = ind + ind_dur_ch(mrks,2) - 1;
                    if ~isempty(ind)
                        try
                            for i = 1:numel(ind)
                                noise(ind(i):indf(i),Idx) = 1;
                            end
                        catch
                            disp('Noise reading problem')
                        end
                    end
                end
            end
            % figure;imagesc(noise)
            
            % ici group channel with the same noise latency
            figureon = m_extractcomponentfigure;
            
            % figure;imagesc(noise')
            % figure;plot(sum(noise,2))
            ind = find((sum(noise,2)./size(noise,2))>nbchminimum );
            ind = find((sum(noise,2)./size(noise,2)) );
            inddiff = diff(ind);
            if isempty(ind)
                disp(['No noisy event found then no component are extracted in file ', fil1])
                break
            end
            try
            idsep = find(inddiff>indconsecutifthreshold);
            if isempty(idsep)
                idstart  =[ind(1)];
                idstop = [ind(end)];
            else
                idstart  =[ind(1);ind(idsep(1:end)+1)];
                idstop =  [ind(idsep(1:end)-1);ind(end)];
            end
            catch
                1
            end
            eventbadstartstop = [idstart,idstop] ;
            for ievent=1:size(eventbadstartstop,1)
                %listchannel to find component
                indt = [eventbadstartstop(ievent,1):eventbadstartstop(ievent,2)];
                if numel(indt)>1
                    if 0 %figureon
                        figure;subplot(3,2,[1,2]); plot(sum(noise,2)/size(noise,2));hold on
                        for i=1:size(eventbadstartstop,1)
                            plot([eventbadstartstop(i,1),eventbadstartstop(i,1)],[0, 1],'b')
                            plot([eventbadstartstop(i,2),eventbadstartstop(i,2)],[0, 1],'r')
                        end
                        plot(indt(round(numel(indt)/2)),1,'x');
                        title(fil1)
                    end
                    
                    plotLst= find(sum(noise(indt,:),1)>0);
                    listok=  find(NIRS.Cf.H.C.ok(plotLst,f)); %bon canaux ?
                    plotLst = plotLst(listok);
                    intensnorm = d(indt,:);
                    % take both wavelengt channel
                    MeasListActplotLst = zeros(size(d,2)./2,1); %combine both lambda
                    id830 = find((plotLst-numel(MeasListActplotLst)) <=0);
                    if ~isempty(id830)
                        MeasListActplotLst(plotLst(id830 ),:)=1;
                    end
                    id690 =  find((plotLst-numel(MeasListActplotLst)) >0);
                    if ~isempty(id690)
                        MeasListActplotLst(plotLst(id690)-numel(MeasListActplotLst),:)=1;
                    end
                    listgood = find(MeasListActplotLst);
                    listgoodPCA = [listgood;listgood+NC/2];
                    %Detrent DATA segment for centrering
                    X = 1:1:size(intensnorm,1);
                    Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
                    Mb2 =  intensnorm(1,:)'; %offset
                    A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
                    spar = intensnorm - A;
                    spar = cat(3,spar(:,1:end/2),spar(:,end/2+1:end));
                    spartmp = spar(:,listgood,:);
                    listgood = find(MeasListActplotLst);
                    if ~isempty(listgood)
                        idwhile = 1;
                        nbPARAFAC = 1; %maximal number of iteration
                        while idwhile==1  %find(std(data(:,listgood,1)./mean(data(:,listgood,1))>0.1))
                            if 0
                                Nc = 1; %force le nombre de component a 2
                                opt(1) = 1e-6;
                                opt(2) = 0;    %initialisation 10 all methods
                                opt(3) = 0;     %plotting
                                opt(5) = 0;     %how often to show fit.
                                opt(6) = 0;     %Max num iterations
                                const = [0 0 0]; % constraints 0- nothing; 1-orthog; 2- nonneg; 3- unimod
                                Oldload{2} = []; % rand(21,3);
                                fixMode = [0 0 0]; %force 2 WAVELENGHT EQUAL
                                weights = []; %Mean for example put one at good time point and zero at noisy one
                                %if too many variance apply multiple PARAFAC
                                [Factors,it,err,corcondia] = parafac(spartmp,Nc,opt,const,Oldload,fixMode,weights);
                                %                             if corcondia<0.9
                                %                                 1
                                %                             end
                                ComponentToKeep = 1;
                            else
                                for itry=1 %si pas vraiment de gain 1 seul essais pour sauver du temps
                                    for icom = 1:job.c_extractcomponent.b_extractnoise_PARAFAC.i_extractnoise_nbPARAFAC
                                        Nc = icom; %force le nombre de component a 2
                                        opt(1) = 1e-6;
                                        opt(2) = 0;    %initialisation 10 all methods
                                        opt(3) = 0;     %plotting
                                        opt(5) = 0;     %how often to show fit.
                                        opt(6) = 0;     %Max num iterations
                                        const = [0 0 0]; % constraints 0- nothing; 1-orthog; 2- nonneg; 3- unimod
                                        Oldload{2} = []; % rand(21,3);
                                        fixMode = [0 0 0]; %force 2 WAVELENGHT EQUAL
                                        weights = []; %Mean for example put one at good time point and zero at noisy one
                                        %if too many variance apply multiple PARAFAC
                                        try
                                            [Factors,it,err,corcondia] = parafac(spartmp,Nc,opt,const,Oldload,fixMode,weights);
                                        catch
                                            1
                                        end
                                        Factorsall{itry,icom} = Factors;
                                        errall(itry,icom) = err;
                                        corcondiaall(itry,icom)=corcondia;
                                    end
                                end
                                [val, idnc]= sort(sum(corcondiaall,1),'descend');
                                % normalisederror to be comparable au concordia
                                errorscale = ( max(errall(:)) - errall)/ max(errall(:))*100;
                                [val, idnc]= sort(sum(errorscale,1)+sum(corcondiaall,1),'descend');
                                if m_extractcomponentfigure
                                    figure;
                                    subplot(3,2,1);hold on
                                    plot(1:numel(corcondiaall),corcondiaall,'color','k')
                                    plot(idnc(1),corcondiaall(idnc(1)),'x','linewidth',4,'color','r')
                                    title('Concordia')
                                    ylabel('Concordia')
                                    xlabel('Nb components')
                                    subplot(3,2,2);hold on
                                    plot(1:numel(corcondiaall),errall,'color','k')
                                    plot(idnc(1),errall(idnc(1)),'x','linewidth',4,'color','r')
                                    title('Error')
                                    ylabel('Error')
                                    xlabel('Nb components')
                                end
                                %try multiple nc
                                
                                
                                
                                
                                %evaluate factor for best concordia and smallest
                                %error
                                for itry=1:size(Factorsall,1)
                                    Factors=Factorsall{itry,idnc(1)};
                                    A = Factors{1};
                                    B = Factors{2};
                                    C = Factors{3};
                                    distC(itry,:) = abs( C(1,:) - C(2,:));
                                    sumA(itry,:) = abs(sum(A));
                                end
                                %Choose among try the smallest time course and the largest sum distance
                                
                                rejeterhightimecourse =(mean(sumA(:)))<= sumA';
                                rejectwavelengthdistance = mean(abs(distC(:))) >= abs(distC)';
                                oneortheother = rejeterhightimecourse|rejectwavelengthdistance;
                                tcmp= indt*1/fs;
                                clear   rejeterhightimecourse  rejectwavelengthdistance
                                ComponentToKeep = find(oneortheother(:,end));% 1:(idnc-1)
                                if m_extractcomponentfigure
                                    subplot(3,2,3);hold on
                                    plot(distC,'x','color','k')
                                    plot(ComponentToKeep,distC(ComponentToKeep),'x','linewidth',4,'color','r')
                                    plot([1,numel(oneortheother)],[mean(abs(distC(:))),mean(abs(distC(:)))] )
                                    title('Distance between wavelength')
                                    xlabel('iteration')
                                    
                                    xlim([0 ,size(distC,2)+1])
                                    subplot(3,2,4);hold on
                                    plot(1:numel(oneortheother),(sumA),'x','color','k')
                                    plot(ComponentToKeep,sumA(ComponentToKeep),'x','linewidth',4,'color','r')
                                    plot([1,numel(oneortheother)],[(mean(sumA(:))),(mean(sumA(:)))] )
                                    title('Abs sum time course ')
                                    xlabel('iteration')
                                    xlim([0 ,size(distC,2)+1])
                                end
                                if m_extractcomponentfigure
                                    subplot(3,2,5);hold on
                                    plot(C,'color','k')
                                    plot(C(:,ComponentToKeep),'linewidth',4,'color','r')
                                    title('Wavelength')
                                    xlabel('Wavelenth')
                                    xlim([0 ,size(distC,2)+1])
                                    subplot(3,2,6);hold on
                                    plot(tcmp,A,'color','k')
                                    plot(tcmp,A(:,ComponentToKeep),'linewidth',4,'color','r')
                                    title(['Time course file', num2str(f)])
                                    xlabel('Time (s)')
                                end
                                clear rejeterhightimecourse rejectwavelengthdistance distC sumA
                            end
                            
                            %                             if 1 %minimal distance between wavelenght et maximum event amplitude
                            %                                 [val,ComponentToKeep]=max([1-abs(C(1,:)-C(2,:)) + sum(A)]);
                            %                             end
                            try
                            Ac = A(:,ComponentToKeep); Bc = B(:,ComponentToKeep); Cc = C(:,ComponentToKeep);
                            [Xm]=nmodel(({Ac,Bc,Cc}));
                            data = cat(3,d(indt,1:end/2),d(indt,end/2+1:end));
                            data(:,listgood,:) = data(:,listgood,:)-Xm;
                            catch
                                1
                            end
                            if 0 %figureon==1
                                figure
                                subplot(3,2,3);plot(Factors{1});subplot(3,2,4);plot(Factors{2});subplot(3,2,5);plot(Factors{3});
                                subplot(3,2,6);plot(Ac);hold on
                            end
                            %             %Optimise the component to keep as lowest distance between wavelength.
                            %figure;plot(std(data(:,listgood,1)./mean(data(:,listgood,1))))
                            %decomposed until zscore
                           
                            if 0 %sum(std(data(:,listgood,1)),1)./mean(data(:,listgood,1))>0.1 & mean(data(:,listgood,1))> max(data(:))*0.05 %test when not much variation except for low intensity channel
                                idwhile = 1; %continue add other parafac not clean yet
                                spartmp = spartmp - Xm;
                                %figure;plot(spartmp(:,:,1))
                            else
                                idwhile = 0; %stop loop
                            end
                            nbPARAFAC =nbPARAFAC+ 1;
                            if nbPARAFAC>1 %only once
                                idwhile=0;  %stop loop
                            end
                            %check the get PARAFAC
                            try
                                load(fullfile(dirout,'SelectedFactors.mat'));
                                newfile = 0;
                            catch
                                %donot exist create the stucture
                                PARCOMP.file= f;
                                PARCOMP.filestr =  fil1;
                                PARCOMP.module  =numel(NIRS.Dt.fir.pp);
                                PARCOMP.modulestr = NIRS.Dt.fir.pp(end).pre;
                                PARCOMP.listgood =  listgood;
                                PARCOMP.indt = indt; %indice de temps.
                                PARCOMP.data =data(:,listgood,:);
                                PARCOMP.Xm = Xm;
                                PARCOMP.FacA = Factors{1};
                                PARCOMP.FacB = Factors{2};
                                PARCOMP.FacC = Factors{3};
                                PARCOMP.ComponentToKeep = ComponentToKeep;
                                labelid  = job.c_extractcomponent.b_extractnoise_PARAFAC.i_extractnoise_labelPARAFAC;%   'MVTPARAFAC';
                                PARCOMP.label= [labelid, sprintf('%03.0f',size(PARCOMP,2)),' ',fil1];
                                PARCOMP.type = 'PARAFAC';
                                B = Factors{2};
                                PARCOMP.topo =   B(:,ComponentToKeep);
                                newfile = 1;
                            end
                            if newfile == 0
                                id = numel(PARCOMP);
                                PARCOMP(id+1).file= f;
                                PARCOMP(id+1).filestr =  fil1;
                                PARCOMP(id+1).module   =numel(NIRS.Dt.fir.pp);
                                PARCOMP(id+1).modulestr  = NIRS.Dt.fir.pp(end).pre;
                                PARCOMP(id+1).listgood =  listgood;
                                PARCOMP(id+1).indt = indt; %indice de temps.
                                PARCOMP(id+1).data = data(:,listgood,:);
                                PARCOMP(id+1).Xm = Xm;
                                PARCOMP(id+1).FacA = Factors{1};
                                PARCOMP(id+1).FacB = Factors{2};
                                PARCOMP(id+1).FacC = Factors{3};
                                PARCOMP(id+1).ComponentToKeep = ComponentToKeep;
                                labelid   = job.c_extractcomponent.b_extractnoise_PARAFAC.i_extractnoise_labelPARAFAC;%
                                PARCOMP(id+1).label= [labelid, sprintf('%03.0f',size(PARCOMP,2)),' ',fil1];
                                PARCOMP(id+1).type = job.c_extractcomponent.b_extractnoise_PARAFAC.i_extractnoise_labelPARAFAC; %'PARAFAC';
                                B = Factors{2};
                                PARCOMP(id+1).topo =  B(:,ComponentToKeep);
                                
                            end
                            
                            
                            save(fullfile(dirout,'SelectedFactors.mat'),'PARCOMP');
                        end
                    end
                    1
                end
            end
        end
    end
    %IDENTIFY GLM COMPONENT
elseif isfield(job.c_extractcomponent,'b_extractcomponent_phys')
    % This applies short-separation regression as described in:
    % (1) Saager and Berger 2005 https://www.ncbi.nlm.nih.gov/pubmed/16211814
    % (2) Scholkmann et al 2014 https://www.ncbi.nlm.nih.gov/pubmed/24622337
    NIRSmat = job.c_extractcomponent.b_extractcomponent_phys.NIRSmat{1};
    load(NIRSmat);
    NC = NIRS.Cf.H.C.N;
    i_extract_pourcentagech= 5; %job.c_extractcomponent.b_extractnoise_PARAFAC.i_extract_pourcentagech;
    [dirout,~,~] = fileparts(NIRSmat);
    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N;
    fs = NIRS.Cf.dev.fs;
    indconsecutifthreshold = 1;     %en dur�e en sample
    nbchminimum = i_extract_pourcentagech/100; %0.05;             %en pourcentage
    fprintf('%s\n','File processed');
    %Open file 1 and check where are the noise marked.
    for f = 1:size(rDtp,1)
        d = fopen_NIR(rDtp{f,1},NC)';
        [dir1,fil1,ext1] = fileparts(rDtp{f});
        vmrk_path = fullfile(dir1,[fil1,'.vmrk']);
        %         handles.file_vmrk = handles.NIRS.Dt.fir.pp(end).p{idfile}; %used
        %         in save noise
        mrk_type_arr = cellstr('bad_step');
        mrks = [];
        ind = [];
        noise =  logical(zeros(size(d)));
        [ind_dur_ch] = read_vmrk_find(vmrk_path,mrk_type_arr);
        if ~isempty(ind_dur_ch)
            maxpoint  = ind_dur_ch(:,1)+ind_dur_ch(:,2);
            badind = find(maxpoint>size(noise,1));
            if ~isempty(badind)
                disp(['Warning file ' vmrk_path ' marker : ' num2str(badind') ' are out of range in the data file'])
                ind_dur_ch(badind,2)=size(noise,2)- ind_dur_ch(badind,1);
            end
            for Idx = 1:size(noise,2)
                mrks = find(ind_dur_ch(:,3)==Idx);
                ind = ind_dur_ch(mrks,1);
                indf = ind + ind_dur_ch(mrks,2) - 1;
                if ~isempty(ind)
                    try
                        for i = 1:numel(ind)
                            noise(ind(i):indf(i),Idx) = 1;
                        end
                    catch
                        disp('Noise reading problem')
                    end
                end
            end
        end
        
        
        % ici group channel with the same noise latency
        
        ind = find((sum(noise,2)./size(noise,2))>nbchminimum );
        inddiff = diff(ind);
        % %              if isempty(ind)    %%%LAURA A MIS EN COMMENTAIRE VOIR  PROPOSITION JUSTE PLUS BAS
        % %                  disp(['No noisy event found then no component are extracted in file ', fil1])
        % %                  break
        % %              end
        %             idsep = find(inddiff>indconsecutifthreshold);
        %             if isempty(idsep)
        %                 idstart  =[ind(1)];
        %                 idstop = [ind(end)];
        %             else
        %                 idstart  =[ind(1);ind(idsep(1:end)+1)];
        %                 idstop =  [ind(idsep(1:end)-1);ind(end)];
        %             end
        %
        %             % add event start
        %           eventbadstartstop = [idstart,idstop] ;
        %           temp =   permute(eventbadstartstop, [2,1])
        %           temp= [1;temp(:); size(noise,1)];
        %           eventint=  reshape(temp,2,numel(temp)/2)
        %           eventgoodstartstop=permute(eventint, [2,1])
        
        if isempty(ind)
            disp(['No noisy event found in file ', fil1])
            temp= [1; size(noise,1)];
            eventint=  reshape(temp,2,numel(temp)/2)
            eventgoodstartstop=permute(eventint, [2,1])
        else
            idsep = find(inddiff>indconsecutifthreshold);
            if isempty(idsep)
                idstart  =[ind(1)];
                idstop = [ind(end)];
            else
                idstart  =[ind(1);ind(idsep(1:end)+1)];
                idstop =  [ind(idsep(1:end)-1);ind(end)];
            end
            
            % add event start
            eventbadstartstop = [idstart,idstop] ;
            temp =   permute(eventbadstartstop, [2,1]);
            temp= [1;temp(:); size(noise,1)];
            eventint=  reshape(temp,2,numel(temp)/2);
            eventgoodstartstop=permute(eventint, [2,1]);
        end
        %%%%%%%jusqu'ici proposition Laura
        
        if 0
            figure;subplot(3,2,[1,2]); plot(sum(noise,2)/size(noise,2));hold on
            for i=1:size(eventbadstartstop,1)
                plot([eventbadstartstop(i,1),eventbadstartstop(i,1)],[0, 1],'b')
                plot([eventbadstartstop(i,2),eventbadstartstop(i,2)],[0, 1],'r')
            end
            
            figure;subplot(3,2,[1,2]); plot(sum(noise,2)/size(noise,2));hold on
            for i=1:size(eventbadstartstop,1)
                plot([eventgoodstartstop(i,1),eventgoodstartstop(i,1)],[0, 1],'b')
                plot([eventgoodstartstop(i,2),eventgoodstartstop(i,2)],[0, 1],'r')
            end
            title('Event form blue to red')
            
        end
        
        
        load(job.c_extractcomponent.b_extractcomponent_phys.f_extractcomponent_physzone{1}   ,'-mat');
        
        Regressorzone = 1;
        ListRegressorZone= [];
        ListChannelZone= [];%(first column regresor channel to average, second row column ch to apply)
        for izone =1:numel(zone.label)
            tmp = upper(zone.label{izone});
            if numel(tmp)>9
                if strcmp(tmp(1:9),'REGRESSOR')
                    ListRegressorZone = [ListRegressorZone,izone];
                end
            end
        end
        
        for iRegressor = 1:numel(ListRegressorZone)
            tmp = upper(zone.label{ListRegressorZone(iRegressor)});
            zoneidentification = tmp(10:end);
            for izone = 1:numel(zone.label)
                tmpzone = upper(zone.label{izone});
                if strcmp( strtrim(zoneidentification), strtrim(tmpzone));
                    ListChannelZone = [ListChannelZone,izone];
                end
            end
        end
        if numel(ListChannelZone) ~= numel(ListRegressorZone)
            msgbox('Regressor Zone and Zone must be in equal number in the zone list')
            return
        end
        
        %Do the regression for each event
        for ievent=1:size( eventgoodstartstop,1)
            if eventgoodstartstop(ievent,2) -eventgoodstartstop(ievent,1) > 10 %don't consider intervall less then a sample
                tstart = eventgoodstartstop(ievent,1);
                tstop = eventgoodstartstop(ievent,2);
                tmpGLM.indt = [tstart,tstop];%Time indice
                tmpGLM.spar = d(tmpGLM.indt(1):tmpGLM.indt(end),:);
                tmpGLM.listgood   = zone.plotLst{1};
                idchgood = [];
                idorder = 1;
                %check good channel for regressor
                badch = NIRS.Cf.H.C.ok(:,f);
                if badch==0 %if all channels are bad, switch to next block! %%ADDED BY LAURA
                    continue %%ADDED BY LAURA
                end %%ADDED BY LAURA
                
                for izoneRegressor = 1:numel(ListRegressorZone)
                    %do the wavelength 1
                    chlistRegressor = zone.plotLst{ListRegressorZone(izoneRegressor)};
                    idbad = find(badch( chlistRegressor)==0); %remove exclude channel from regressor
                    if ~isempty(idbad)
                        chlistRegressor(idbad) = [];
                        if isempty(chlistRegressor)
                            disp(['No good channel in the regressor zone please verify your zone ', zone.label{ListRegressorZone(izoneRegressor)}])
                            continue %break %%%%MIS EN COMMENTAIRE PAR LAURA
                        end
                    end
                    chlistApply = zone.plotLst{ListChannelZone(izoneRegressor)};
                    if job.c_extractcomponent.b_extractcomponent_phys.m_extractcomponent_physzone == 0 %mean if many channel
                        XmeanSD = nanmean(tmpGLM.spar(:,chlistRegressor),2)-  nanmean(nanmean(tmpGLM.spar(:,chlistRegressor),2)); %center to zero ? necessary ?
                    elseif  job.c_extractcomponent.b_extractcomponent_phys.m_extractcomponent_physzone == 1  %PCA  if many channel
                        d1 = tmpGLM.spar(:,chlistRegressor);
                        c = d1'*d1;
                        [v,s,foo]=svd(c);
                        svs = diag(s);
                        u = d1*v*inv(s);
                        lstSV = 1;
                        Xm =  u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)';
                        XmeanSD = Xm(:,1)-mean(Xm(:,1));
                    end
                    
                    for ich = 1:numel(chlistApply)
                        Xlong = tmpGLM.spar(:,chlistApply(ich));
                        tmpGLM.beta(idorder) = dot(XmeanSD ,Xlong)/dot(XmeanSD ,XmeanSD );
                        idchgood = [idchgood, chlistApply(ich)];
                        Xcorr = Xlong - tmpGLM.beta(idorder).*XmeanSD;
                        %figure;plot(Xlong,'b','displayname', 'Initial');hold on;plot(XmeanSD,'r','displayname','Short distance');plot(Xcorr,'k','displayname','corrected')
                        tmpGLM.Xm(:,idorder) = tmpGLM.beta(idorder).*XmeanSD;
                        idorder = idorder + 1;
                    end
                    %do the wavelength 2
                    chlistRegressor =  chlistRegressor + size(d,2)/2;
                    chlistApply =  chlistApply + size(d,2)/2;;
                    if job.c_extractcomponent.b_extractcomponent_phys.m_extractcomponent_physzone == 0 %mean if many channel
                        XmeanSD = nanmean(tmpGLM.spar(:,chlistRegressor),2)-  nanmean(nanmean(tmpGLM.spar(:,chlistRegressor),2)); %center to zero ? necessary ?
                    elseif  job.c_extractcomponent.b_extractcomponent_phys.m_extractcomponent_physzone == 1  %PCA  if many channel
                        d1 = tmpGLM.spar(:,chlistRegressor);
                        c = d1'*d1;
                        [v,s,foo]=svd(c);
                        svs = diag(s);
                        u = d1*v*inv(s);
                        lstSV = 1
                        Xm =  u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)';
                        XmeanSD = Xm(:,1)-mean(Xm(:,1));
                    end
                    
                    for ich = 1:numel(chlistApply)
                        Xlong = tmpGLM.spar(:,chlistApply(ich)); -mean( tmpGLM.spar(:,chlistApply(ich)));
                        tmpGLM.beta(idorder) = dot(XmeanSD ,Xlong)/dot(XmeanSD ,XmeanSD );  %b = regress( Xlong , XmeanSD) identique
                        idchgood = [idchgood, chlistApply(ich)];
                        Xcorr = Xlong - tmpGLM.beta(idorder).*XmeanSD;
                        %  figure;plot(Xlong,'b','displayname', 'Initial');hold on;plot(XmeanSD,'r','displayname','Short distance');plot(Xcorr,'k','displayname','corrected')
                        tmpGLM.Xm(:,idorder) = tmpGLM.beta(idorder).*XmeanSD;
                        idorder = idorder + 1;
                    end
                    
                end
                
                tmpGLM.idreg = 1;
                tmpGLM.selected  = 1;
                label = 'SHORTGLM';
                tmpGLM.listgood =  idchgood ;
                
                try
                    load(fullfile(dirout,'SelectedFactors.mat'));
                    newfile = 0;
                catch
                    clear PARCOMP
                    %donot exist create the stucture
                    PARCOMP.file= f;
                    PARCOMP.filestr =  sprintf('Bloc%03.0f',ievent);
                    PARCOMP.module  = numel(NIRS.Dt.fir.pp);
                    PARCOMP.modulestr = NIRS.Dt.fir.pp(end).pre;
                    PARCOMP.listgood =  tmpGLM.listgood;
                    PARCOMP.beta = tmpGLM.beta;
                    PARCOMP.std= 0;
                    PARCOMP.AUX= 0;
                    PARCOMP.indt = tmpGLM.indt(1):tmpGLM.indt(end); %indice de temps.
                    PARCOMP.data = d(tmpGLM.indt(1):tmpGLM.indt(end),tmpGLM.listgood) ;%data(:,listgood,:);
                    PARCOMP.Xm = tmpGLM.Xm;
                    
                    PARCOMP.ComponentToKeep = tmpGLM.selected;
                    PARCOMP.idreg = tmpGLM.idreg;
                    PARCOMP.label= ['GLM',label , sprintf('%03.0f',size(PARCOMP,2))];
                    PARCOMP.type = 'GLM';
                    PARCOMP.topo =  tmpGLM.beta(tmpGLM.selected,:);
                    newfile = 1;
                end
                if newfile == 0
                    id = numel(PARCOMP);
                    PARCOMP(id+1).file= f;
                    PARCOMP(id+1).filestr =  sprintf('Bloc%03.0f',ievent);
                    PARCOMP(id+1).module  = numel(NIRS.Dt.fir.pp);
                    PARCOMP(id+1).modulestr = NIRS.Dt.fir.pp(end).pre;
                    PARCOMP(id+1).listgood =  tmpGLM.listgood;
                    PARCOMP(id+1).beta = tmpGLM.beta;
                    PARCOMP(id+1).std= 0;
                    PARCOMP(id+1).AUX= 0;
                    PARCOMP(id+1).indt = tmpGLM.indt(1):tmpGLM.indt(end); %indice de temps.
                    PARCOMP(id+1).data = d(tmpGLM.indt(1):tmpGLM.indt(end),tmpGLM.listgood) ;%data(:,listgood,:);
                    PARCOMP(id+1).Xm = tmpGLM.Xm;
                    PARCOMP(id+1).ComponentToKeep = tmpGLM.selected;
                    PARCOMP(id+1).idreg = tmpGLM.idreg;
                    PARCOMP(id+1).label= ['GLM',label , sprintf('%03.0f',size(PARCOMP,2))];
                    PARCOMP(id+1).type = 'GLM';
                    PARCOMP(id+1).topo =  tmpGLM.beta(tmpGLM.selected,:);
                    newfile = 1;
                end
                save(fullfile(dirout,'SelectedFactors.mat'),'PARCOMP');
                clear tmpGLM
            end
        end
    end
    
    
elseif isfield(job.c_extractcomponent,'b_extractcomponent_glm')
  for ixlsfile = 1:size(job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist,1)
       [currentpath,~,~] = fileparts(job.c_extractcomponent.b_extractcomponent_glm.NIRSmat{1});
      if isempty(job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{ixlsfile})         
          job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{ixlsfile} = fullfile(currentpath,'ExtractHRF.xlsx');
      end
    [~,~,ext] =fileparts(job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{ixlsfile});
    if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
        try
            [data, text, rawData] = xlsread(job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{ixlsfile});
            [dirxls,filexls,extxls] = fileparts(job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{ixlsfile});
            id.Regressor = [];
        catch
            try
            [data, text, rawData] = readtxtfile_asxlsread(job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{ixlsfile});
            [dirxls,filexls,extxls] = fileparts(job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{ixlsfile});
            id.Regressor = []; 
            catch
                disp(['Could not open: ', job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{ixlsfile}]);
            end
        end
        
    elseif strcmp(ext,'.txt')
        [data, text, rawData] = readtxtfile_asxlsread(job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{ixlsfile});
        [dirxls,filexls,extxls] = fileparts(job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{ixlsfile});
        id.Regressor = [];
    end
    for icol=1:size(rawData,2)
        if strcmp(upper(deblank(rawData{1,icol})),deblank(upper('NIRS.mat folder')))
            id.NIRSDtp = icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('File'))
            id.fileDtp =  icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('tStart'))
            id.startDtp =  icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('tStop'))
            id.stopDtp  =  icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('Label'))
            id.labelDtp =  icol;
        end
        for i=1:21 %up to a max of 20 regressor
            if strcmp(upper(deblank(rawData{1,icol})),upper(['X',num2str(i-1)]))
                id.Regressor = [id.Regressor,icol];
            end
        end
    end
    try
        NIRSDtp = rawData(2:end,id.NIRSDtp );
        fileDtp = rawData(2:end,2);
        chDtp = 'HBO' ; rawData(2:end,3); %do both to modify except ....
        % trigDtp = rawData(2:end,4); non utilis�
        startDtp = rawData(2:end,id.startDtp);
        stopDtp = rawData(2:end,id.stopDtp);
        labelDtp = rawData(2:end,id.labelDtp);
        AUXid = rawData(2:end,id.Regressor );
    catch
        msgbox(['Please verify the GLM extract xls file : ', job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{ixlsfile},...
            ' have the following column information : NIRS.mat folder, File, tStart, tStop, label and Xn regressors'])
    end
    
    try
       regress(ones(10,1),ones(10,1));
    catch
        disp('Uncomplete Extract GLM, Please Install Matlab Statistics and Machine Learning Toolbox or regress.m function')
        NIRSmat = fullfile(NIRSDtp{1},'NIRS.mat');
        out.NIRSmat = {NIRSmat};
        return
    end
    %Load data and aux for regression
    for ievent=1:size(NIRSDtp,1)
        %try
        warning off
        tmp=NIRSDtp{ievent};
        try
        if strcmp(tmp(end-7:end),'NIRS.mat')
            NIRSmat = NIRSDtp{ievent};
            NIRSDtp{ievent} = tmp(1:end-8);
        else
            NIRSmat = fullfile(NIRSDtp{ievent},'NIRS.mat');
        end
        catch
            NIRSmat = fullfile(NIRSDtp{ievent},'NIRS.mat');
        end
        %try
       
        disp(['load',NIRSmat])
        load(NIRSmat);
        
           
       
        Regressorlist=AUXid{ievent,:};
        disp(NIRSmat);
        NC = NIRS.Cf.H.C.N;
        fDtp = NIRS.Dt.fir.pp(end).p;
        d1 = fopen_NIR(fDtp{fileDtp{ievent}},NC)';
        [dir1,fil1,ext] = fileparts(fDtp{fileDtp{ievent}});
        vmrk_path = fullfile(dir1,[fil1,'.vmrk']);
        %         handles.file_vmrk = handles.NIRS.Dt.fir.pp(end).p{idfile}; %used
        %         in save noise
        mrk_type_arr = cellstr('bad_step');
        mrks = [];
        ind = [];
        noise =  logical(zeros(size(d1)));
        [ind_dur_ch] = read_vmrk_find(vmrk_path,mrk_type_arr);
        if ~isempty(ind_dur_ch)
            maxpoint  = ind_dur_ch(:,1)+ind_dur_ch(:,2);
            badind = find(maxpoint>size(noise,1));
            if ~isempty(badind)
                disp(['Warning file ' vmrk_path ' marker : ' num2str(badind') ' are out of range in the data file']);
                ind_dur_ch(badind,2)=size(noise,2)- ind_dur_ch(badind,1);
            end
            for Idx = 1:size(noise,2)
                mrks = find(ind_dur_ch(:,3)==Idx);
                ind = ind_dur_ch(mrks,1);
                indf = ind + ind_dur_ch(mrks,2) - 1;
                if ~isempty(ind)
                    try
                        for i = 1:numel(ind)
                            noise(ind(i):indf(i),Idx) = 1;
                        end
                    catch
                        disp('Noise reading problem')
                    end
                end
            end
        end
        if 1 %trcmp(upper(chDtp{ievent}),'HBO')
            tmpGLM.listgood = 1:NC; %do all
        elseif strcmp(upper(chDtp{ievent}),'HBR')
            tmpGLM.listgood = (1:(NC/2)) + NC/2;
        end
        
        tHRF = 1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:size(d1,1)*1/NIRS.Cf.dev.fs;
        fsNIRS = NIRS.Cf.dev.fs;
        tstart = find(tHRF<=startDtp{ievent});
        
        if isempty(tstart) 
            tstart = 1;
        end
        tstop = find(tHRF<= stopDtp{ievent});
        tmpGLM.indt = [tstart(end),tstop(end)];%Time indice
        tmpGLM.spar = d1(tmpGLM.indt(1):tmpGLM.indt(end),:); %fix weard edge offset by detrending first
%         tmpGLM.spar = tmpGLM.spar - detrend(tmpGLM.spar);
%         figure;plot(tmpGLM.spar)
% figure;plot( d1(tmpGLM.indt(1):tmpGLM.indt(end),:))
        %add labelisbad
        pourcentagenoise = sum(sum(noise(tmpGLM.indt(1):tmpGLM.indt(end),:)))./numel(d1(tmpGLM.indt(1):tmpGLM.indt(end),:));
        if pourcentagenoise > (job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.i_glmlist_autoexport_labelbad_threshold/100)
            labelisbad = 'bad';
        else
             labelisbad = 'ok';
        end
        pourcentagenoisebychHbO = sum(noise(tmpGLM.indt(1):tmpGLM.indt(end),1:end/2))./numel(tmpGLM.indt(1):tmpGLM.indt(end))*100;
        pourcentagenoisebychHbR = sum(noise(tmpGLM.indt(1):tmpGLM.indt(end),end/2+1:end))./numel(tmpGLM.indt(1):tmpGLM.indt(end))*100;

        pourcentagenoise = pourcentagenoise*100;
        if job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.i_glmlist_autoexport_nan_chrejected==1
            %(sum(noise(tmpGLM.indt(1):tmpGLM.indt(end),:))./numel(d1(tmpGLM.indt(1):tmpGLM.indt(end),1)))
            ch_remove_for_int =(sum(noise(tmpGLM.indt(1):tmpGLM.indt(end),:))./numel(d1(tmpGLM.indt(1):tmpGLM.indt(end),1)))< (job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.i_glmlist_autoexport_labelbad_threshold/100);
            %reshape(sum(noise(tmpGLM.indt(1):tmpGLM.indt(end),:)),NC/2 ,2) 
        end 
            iRegressor =  2; 
        if isfield(NIRS.Dt,'AUX')
            for iAUX = 1:numel(NIRS.Dt.AUX)
                nameAUX = NIRS.Dt.AUX(iAUX).pp(end).p{fileDtp{ievent}};
                if isfield(NIRS.Dt.AUX(iAUX).pp(end),'sync_timesec')
                    tstartf = NIRS.Dt.AUX(iAUX).pp(end).sync_timesec{fileDtp{ievent}};
                else
                    tstartf = 0;
                    disp('No segmentation have been made ensure that aux synchronisation are ok')
                end
                tstopf = tstartf+tHRF(end);
                [pathtmp,filetmp,exttmp]=fileparts(nameAUX);
                [data,infoBV,label,ind_dur_ch] = fopen_EEG(nameAUX, tstartf, tstopf);
                
                %             [pathtmp,filetmp,exttmp]=fileparts(nameAUX);
                %             [data,infoBV,label,ind_dur_ch] = fopen_EEG(nameAUX);
                for ich=1:numel(infoBV.name_ele)
                    tmpGLM.AUX.label{iRegressor} =[filetmp,' ',infoBV.name_ele{ich}];
                    fsAUX =1/(infoBV.SamplingInterval/1000000); %Frequence echantillonage Hz
                    tmp = data(:,ich);
                    [p,q] = rat(fsNIRS/fsAUX,0.0001);
                    
                    %tmpr=resample( tmp , p, q); %ancien Julie - a cause
                    %d'un jump au d�but fin, Laura l'a modifie par
                    %downsample - environ le meme fonctionnement - voir la
                    %photo Guide DownsampleEEgdata dans GITHUB student
                    
                    %not always integer factor then resample is be back as
                    %default option it crash with downsample %Julie
                   if ~isempty(tmp)
                    tmpr=resample( tmp , p, q); %copy resample function here to avoid using different version of the resample function which have change with years... 
                   end
                    
                    % we cut aux to the data initial size
                    if numel(tmpr)<numel(tHRF)
                        nplus = numel(tHRF)-numel(tmpr);
                        try
                            tmpr = [tmpr;tmpr(end-nplus:end) ];
                        catch
                            msgbox('Too short regressor')
                        end
                    elseif numel(tmpr)>numel(tHRF)
                        tmpr = tmpr(1:numel(tHRF));
                    end
                    %we cut to fit the spar selection
                    tmpGLM.AUX.fs{iRegressor} = fsNIRS;
                    tmpGLM.AUX.data{iRegressor} = tmpr(tmpGLM.indt(1):tmpGLM.indt(end));
                    clear tmpr
                    tmpGLM.AUX.view{iRegressor} = 1;
                    iRegressor = iRegressor +1;
                end
            end
        end
        tmpGLM.AUX.label{1} = 'Constant';
        tmpGLM.AUX.fs{1} = fsNIRS;
        tmpGLM.AUX.data{1} = ones(numel(tmpGLM.indt(1):tmpGLM.indt(end) ),1);
        tmpGLM.AUX.view{1} = 1;
        Regressorlist = AUXid(ievent,:)';
        idreg = [];
        ListRegressorZone = [];
        for iReg = 1:numel(Regressorlist)
            for ilist = 1:numel(tmpGLM.AUX.label)
                if numel(strfind(tmpGLM.AUX.label{ilist},Regressorlist{iReg})) %strcmp(tmpGLM.AUX.label{ilist},Regressorlist{iReg})
                    idreg = [idreg,ilist];
                    disp([Regressorlist{iReg},' AUX regressor find'])
                end
            end
            
            if numel(Regressorlist{iReg})>4
                checkzone = Regressorlist{iReg};
                if strcmp(upper(checkzone(end-3:end)),'ZONE')
                    try
                        load(Regressorlist{iReg} ,'-mat'); %try the fullfile correct
                        disp([Regressorlist{iReg},' zone regressor find'])
                    catch
                        try
                            load(fullfile(dirxls, Regressorlist{iReg} ),'-mat'); %try in the excel folder
                        catch
                            disp(['Regressor zone :',  Regressorlist{iReg} ,' or ',fullfile(dirxls,Regressorlist{iReg}) , ' could not be load'])
                            return
                        end
                    end
                    idreg = [idreg,0]; %if idreg = 0 regressor is a zone of the current channel
                end
            end
        end
        
        %PERFORM REGRESSION
        tmpGLM.idreg = idreg;
        X = [];%ones(size(PMI{1}.tmpGLM.spar,1),1);
        Regressorzone = 0;
        for ireg = 1:numel(idreg)
            if idreg(ireg)==0
                Regressorzone = 1;
                ListRegressorZone= [];
                ListChannelZone= [];%(first column regresor channel to average, second row column ch to apply)
                for izone =1:numel(zone.label)
                    tmp = upper(zone.label{izone});
                    if numel(tmp)>9
                        if strcmp(tmp(1:9),'REGRESSOR')
                            ListRegressorZone = [ListRegressorZone,izone];
                        end
                    end
                end
                
                for iRegressor = 1:numel(ListRegressorZone)
                    tmp = upper(zone.label{ListRegressorZone(iRegressor)});
                    zoneidentification = tmp(10:end);
                    for izone = 1:numel(zone.label)
                        tmpzone = upper(zone.label{izone});
                        if strcmp( strtrim(zoneidentification), strtrim(tmpzone))
                            ListChannelZone = [ListChannelZone,izone];
                        end
                    end
                end
                if numel(ListChannelZone) ~= numel(ListRegressorZone)
                    msgbox('Regressor Zone and Zone must be in equal number in the zone list')
                    return
                end
                
            else
                X = [X,    tmpGLM.AUX.data{idreg(ireg)}];
            end
        end
        
        
        Xtmp = X;
        XmRegressor = zeros(size(tmpGLM.spar));
        if numel(ListRegressorZone)
            beta = zeros(size(X,2)+1,numel(tmpGLM.listgood));
            bstd = zeros(size(X,2)+1,numel(tmpGLM.listgood));
            for izoneRegressor = 1:numel(ListRegressorZone)
                for iconc = 1:2 %HbO and HbR
                    chlistRegressor = zone.plotLst{ListRegressorZone(izoneRegressor)};
                    chlistApply = zone.plotLst{ListChannelZone(izoneRegressor)};
                    if iconc==2
                        chlistRegressor =  chlistRegressor+ NC/2;
                        chlistApply = chlistApply + NC/2;
                    end
                    Xmean = nanmean(tmpGLM.spar(:,chlistRegressor),2);% -  nanmean(nanmean(tmpGLM.spar(:,chlistRegressor),2)); %center to zero
                    % figure;plot(Xmean)
                    X = [Xtmp,    Xmean];
                    for ich = 1:numel(chlistApply)
                        idch = chlistApply(ich);
                        y = tmpGLM.spar(:,idch);
                        %figure;plot(PMI{1}.tmpGLM.spar)
                        if sum(isnan(y))
                            b = zeros(size(X,2),1);
                            b = nan(size(X,2),1);
                            beta(:,idch) = nan; % 0;
                            bstd(:,idch) = nan; % 0;
                            R2(:,idch) = nan; %0
                        else
                            [b,bint,r,rint,stats]=  regress(y,X);
                            beta(:,idch) = b;
                            bstd(:,idch) = stats(4);
                            R2(:,idch) =  stats(1);
                        end
                        XmRegressor(:,idch) = Xmean*b(end);
                    end
                end
            end
        elseif ~numel(ListRegressorZone) %pas de zone de regression juste des auxiliaires
            %PERFORM REGRESSION
            tmpGLM.idreg = idreg;
            X = [];%ones(size(PMI{1}.tmpGLM.spar,1),1);
            for ireg = 1:numel(idreg)
                tmpGLM.AUX.label{idreg(ireg)};
                X = [X,    tmpGLM.AUX.data{idreg(ireg)}];
            end
            beta = zeros(size(X,2),numel(tmpGLM.listgood));
            bstd = zeros(size(X,2),numel(tmpGLM.listgood));
            for ich = 1:numel(tmpGLM.listgood)
                idch = tmpGLM.listgood(ich);
                y = tmpGLM.spar(:,idch);
                %figure;plot(PMI{1}.tmpGLM.spar)
                if sum(isnan(y))
                    b = zeros(size(X,2),1);
                    beta(:,idch) = 0;
                    bstd(:,idch) = 0;
                    R2(:,idch) = 0; 
                else
                    [b,bint,r,rint,stats]=  regress(y,X);
                    beta(:,idch) = b;
                    bstd(:,idch) = stats(4);
                    R2(:,idch) =  stats(1);
                end
                
            end
            
        end
        tmpGLM.beta = beta;
        tmpGLM.std = bstd;
        tmpGLM.R2 = R2;
        
        
        if  Regressorzone==1
            XmRegressorgood = XmRegressor(:,tmpGLM.listgood);
        end
        
        %Save in the comp list
        for iselected = 1:size(beta,1)
            beta =tmpGLM.beta ;
            idreg = tmpGLM.idreg;
            tmpGLM.selected  = iselected;
            if idreg(iselected)==0
                % Xmean(:) =1
                Xm = XmRegressorgood;
                label = Regressorlist{iselected};
                [filepath,name,ext] = fileparts(Regressorlist{iselected});
                label =  [name,ext];                
            else
               % Xm = tmpGLM.AUX.data{idreg(iselected)}*beta(iselected,:);
                Xm = tmpGLM.AUX.data{idreg(iselected)}; %compress multiply by beta later
                label = tmpGLM.AUX.label{idreg(iselected)};
            end
            try
                if isfield(job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport,'b_extractglmlist_autoexport_yes')                  
                   if job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.m_replaceglmlist == 1;
                         delete(fullfile(NIRSDtp{ievent},'SelectedFactors.mat'));
                         job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.m_replaceglmlist = 0;
                   end
                end
                load(fullfile(NIRSDtp{ievent},'SelectedFactors.mat'));
                newfile = 0;
            catch
                clear PARCOMP
                %donot exist create the stucture
                PARCOMP.file= fileDtp{ievent};
                PARCOMP.filestr =  sprintf('Bloc%03.0f',fileDtp{ievent});
                PARCOMP.module  = numel(NIRS.Dt.fir.pp);
                PARCOMP.modulestr = NIRS.Dt.fir.pp(end).pre;
                PARCOMP.listgood =  tmpGLM.listgood;
                PARCOMP.beta = tmpGLM.beta;
                PARCOMP.std =  tmpGLM.std;
                PARCOMP.R2 = tmpGLM.R2;
                PARCOMP.AUX =  tmpGLM.AUX;
                PARCOMP.indt = tmpGLM.indt(1):tmpGLM.indt(end); %indice de temps.
                PARCOMP.data = [];%d1(tmpGLM.indt(1):tmpGLM.indt(end),tmpGLM.listgood) ;%data(:,listgood,:);
                PARCOMP.Xm = Xm;
                %to large data file with data and xm remove them
                PARCOMP.ComponentToKeep = tmpGLM.selected;
                PARCOMP.idreg = tmpGLM.idreg;
                labelid  = labelDtp{ievent} ;
                PARCOMP.label= [ labelisbad, labelid,'_', label,'GLM', sprintf('%03.0f',size(PARCOMP,2))];
                disp([ labelisbad, labelid,'_',label, 'GLM', sprintf('%03.0f',size(PARCOMP,2))]);
                PARCOMP.type = 'GLM';
                PARCOMP.topo =  beta(iselected,:);
                newfile = 1;
            end
            if newfile == 0
                id = numel(PARCOMP);
                %donot exist create the stucture
                %donot exist create the stucture
                PARCOMP(id+1).file= fileDtp{ievent};
                PARCOMP(id+1).filestr =  sprintf('Bloc%03.0f',fileDtp{ievent});
                PARCOMP(id+1).module  = numel(NIRS.Dt.fir.pp);
                PARCOMP(id+1).modulestr = NIRS.Dt.fir.pp(end).pre;
                PARCOMP(id+1).listgood =  tmpGLM.listgood;
                PARCOMP(id+1).beta = tmpGLM.beta;
                PARCOMP(id+1).std =  tmpGLM.std;
                PARCOMP(id+1).R2 = tmpGLM.R2;
                PARCOMP(id+1).AUX =  tmpGLM.AUX;
                PARCOMP(id+1).indt = tmpGLM.indt(1):tmpGLM.indt(end); %indice de temps.
                PARCOMP(id+1).data = [];%d1(tmpGLM.indt(1):tmpGLM.indt(end),tmpGLM.listgood) ;%data(:,listgood,:);
                PARCOMP(id+1).Xm = Xm;
                PARCOMP(id+1).ComponentToKeep = tmpGLM.selected;
                PARCOMP(id+1).idreg = tmpGLM.idreg;
                labelid  = labelDtp{ievent} ;
                PARCOMP(id+1).label= [ labelisbad, labelid,'_',label,'GLM', sprintf('%03.0f',size(PARCOMP,2))];
                disp([ labelisbad, labelid,'_',label ,'GLM', sprintf('%03.0f',size(PARCOMP,2))]);
                PARCOMP(id+1).type = 'GLM';
                PARCOMP(id+1).topo =  beta(iselected,:);
                newfile = 1;
            end
            save(fullfile(NIRSDtp{ievent},'SelectedFactors.mat'),'PARCOMP');
            try 
                if isfield(job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport,'b_extractglmlist_autoexport_yes')                 
                    srsfile = NIRSmat;
                  
                    listgood=PARCOMP(end).listgood;
                      listbad_ch = ones(numel(listgood),1);
                     listbad_trial = ones(numel(listgood),1);
                    if job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.i_glmlist_autoexport_nan_chrejected == 1
                        listbad_ch=  NIRS.Cf.H.C.ok(listgood);
                    end
                    if job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.i_glmlist_autoexport_nan_int_rejected == 1
                      if sum(~ch_remove_for_int)>0
                        disp(['Event ', num2str(ievent),' Reject ', num2str(sum(~ch_remove_for_int)),'/',num2str(numel(ch_remove_for_int)), 'CH  noise  >' ,num2str(job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.i_glmlist_autoexport_labelbad_threshold),'%']);
                        listbad_trial =ch_remove_for_int(listgood)';
                      end
                    end 
                    
                    topo = PARCOMP(end).topo;
                    pathoutlist = job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.f_extractglmlist_autoexport_yes{1};
                    if ~isdir(pathoutlist)
                        mkdir(pathoutlist);
                        disp(['Create ' pathoutlist])
                    end
                    if isempty(job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.f_extractglmlist_autoexport_yes_ChannelList{1})
                        job.NIRSmat = job.c_extractcomponent.b_extractcomponent_glm.NIRSmat;
                        nirs_run_createseedlist(job);
                        job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.f_extractglmlist_autoexport_yes_ChannelList{1} = fullfile(currentpath,'channellist.txt');
                    end
                    ChannelListfile = job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.f_extractglmlist_autoexport_yes_ChannelList{1};
                    [listHBOch, listHBRch, listnameHbO, listnameHbR , zonelist]= findchinChannelList(NIRS, ChannelListfile,PARCOMP(end).listgood);
                    Xiselected = job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.i_glmlist_autoexport_Xi + 1;
                     listbad = listbad_trial(:)&listbad_ch(:);
                    if job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.m_glmlist_autoexport_HbO==1
                          if find(iselected==Xiselected) %only beta 1
                                A = nan(size(listHBOch,1),1);
                                idok = find(listbad(listHBOch)); 
                                A(idok,1) = topo(listHBOch(idok));
                                save(fullfile(pathoutlist,['TopoHbO',PARCOMP(end).label,'event',sprintf('%03.0f',ievent),'.mat']),'A' ,'zonelist','srsfile', 'pourcentagenoise', 'pourcentagenoisebychHbO' );
                                disp(['Create export: load(''', fullfile(pathoutlist,['TopoHbO',PARCOMP(end).label,'event',sprintf('%03.0f',ievent),'.mat'')'])]);
                          end
                    elseif job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.m_glmlist_autoexport_HbO==2
                          if find(iselected==Xiselected) %only beta 1
                            A = nan(size(listHBRch,1),1);
                            idok = find(listbad(listHBRch));
                            A(idok,1) = topo(listHBRch(idok));
                            save(fullfile(pathoutlist,['TopoHbR',PARCOMP(end).label,'event', sprintf('%03.0f',ievent),'.mat']),'A' ,'zonelist','srsfile', 'pourcentagenoise','pourcentagenoisebychHbR' );
                            disp(['Create export: load(''', fullfile(pathoutlist,['TopoHbR',PARCOMP(end).label,'event', sprintf('%03.0f',ievent),'.mat'')'])])
                          end
                    elseif job.c_extractcomponent.b_extractcomponent_glm.c_extractglmlist_autoexport.b_extractglmlist_autoexport_yes.m_glmlist_autoexport_HbO==3
                        
                          if find(iselected==Xiselected) %only beta 1
                                A = nan(size(listHBOch,1),1);
                                idok = find(listbad(listHBOch)); 
                                A(idok,1) = topo(listHBOch(idok));
                                save(fullfile(pathoutlist,['TopoHbO',PARCOMP(end).label,'event',sprintf('%03.0f',ievent),'.mat']),'A' ,'zonelist','srsfile','pourcentagenoise', 'pourcentagenoisebychHbO' );
                                disp(['Create export: load(''', fullfile(pathoutlist,['TopoHbO',PARCOMP(end).label,'event',sprintf('%03.0f',ievent),'.mat'')'])]);
                          end
                        
                        
                         if find(iselected==Xiselected) %only beta 1
                            A = nan(size(listHBRch,1),1);
                            idok = find(listbad(listHBRch));
                            A(idok,1) = topo(listHBRch(idok));
                            save(fullfile(pathoutlist,['TopoHbR',PARCOMP(end).label,'event', sprintf('%03.0f',ievent),'.mat']),'A' ,'zonelist','srsfile','pourcentagenoise', 'pourcentagenoisebychHbR' );
                            disp(['Create export: load(''', fullfile(pathoutlist,['TopoHbR',PARCOMP(end).label,'event', sprintf('%03.0f',ievent),'.mat'')'])])
                          end
                    end
                end
            catch
                disp('WARNING topo could not be exported verify channel list')
            end
        end
        %disp(['Error unable to GLM on ' , NIRSmat])
%        catch
%            disp('XLS line do not contain valid NIRS.mat file')
%        end
    end
    try
        filenamexls = fullfile(NIRSDtp{ievent},['Export ', labelid,'.xlsx']);  
        A = {'NIRS.mat folder', 'Type', 'Label', 'Channel List','Name'};
        %CREATE XLS here
        namefilter=[labelid,'_',Regressorlist{1}]; 
        disp(['Suggested filter name ',  namefilter,' to export'])
        VAL = {NIRSDtp{ievent},'GLM', namefilter, [NIRSDtp{ievent}, 'channellist.txt'],namefilter};
        try
        xlswrite(filenamexls,[A;VAL]);
        disp(['Create xls example file: ',filenamexls,' to help you to configure the Export list function to export component'])
        catch
        writetxtfile(filenamexls,[A;VAL]);
        disp(['Create xls example file: ',filenamexls,' to help you to configure the Export list function to export component'])
        end
    catch
    end
end
elseif isfield(job.c_extractcomponent,'b_extractcomponent_PARAFAC')
    [~,~,ext] =fileparts(job.c_extractcomponent.b_extractcomponent_PARAFAC.f_component_PARAFAClist{1});
    if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
        try
        [data, text, rawData] = xlsread(job.c_extractcomponent.b_extractcomponent_PARAFAC.f_component_PARAFAClist{1});
        catch
        [data, text, rawData] = readtxtfile_asxlsread(job.c_extractcomponent.b_extractcomponent_PARAFAC.f_component_PARAFAClist{1});
        end 
    elseif strcmp(ext,'.txt')
        [data, text, rawData] = readtxtfile_asxlsread(job.c_extractcomponent.b_extractcomponent_PARAFAC.f_component_PARAFAClist{1});
    end
    
    for icol=1:size(rawData,2)
        if strcmp(upper(deblank(rawData{1,icol})),deblank(upper('NIRS.mat folder')))
            id.NIRSDtp = icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('File'))
            id.fileDtp =  icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('tStart'))
            id.startDtp =  icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('tStop'))
            id.stopDtp  =  icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('Label'))
            id.labelDtp =  icol;
        end
    end
    
    NIRSDtp = rawData(2:end,id.NIRSDtp );
    fileDtp = rawData(2:end,2);
    chDtp = 'HBO' ; rawData(2:end,3); %do both to modify except ....
    % trigDtp = rawData(2:end,4); non utilis�
    %     chDtp = rawData(2:end,3);
    %     trigDtp = rawData(2:end,4);
    startDtp = rawData(2:end,id.startDtp);
    stopDtp = rawData(2:end,id.stopDtp);
    labelDtp = rawData(2:end,id.labelDtp);
    
    
    %Load data and aux for regression
    for ievent=1:size(NIRSDtp,1)
        NIRSDtp{ievent};
        NIRSmat = fullfile(NIRSDtp{ievent},'NIRS.mat');
        load(NIRSmat);
        NC = NIRS.Cf.H.C.N;
        fDtp = NIRS.Dt.fir.pp(end).p;
        d = fopen_NIR(fDtp{fileDtp{ievent}},NC)';
        
        listgood = 1:(NC/2);
        
        tHRF = 1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:size(d,1)*1/NIRS.Cf.dev.fs;
        fsNIRS = NIRS.Cf.dev.fs;
        tstart = find(tHRF<=startDtp{ievent});
        tstop = find(tHRF<= stopDtp{ievent});
        indt = [tstart(end),tstop(end)];%Time indice
        intensnorm = d(indt(1):indt(end),:);
        
        %Detrent DATA segment for centrering
        X = 1:1:size(intensnorm,1);
        Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
        Mb2 =  intensnorm(1,:)'; %offset
        A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
        spar = intensnorm - A;
        spar = cat(3,spar(:,1:end/2),spar(:,end/2+1:end));
        listgood=  find(NIRS.Cf.H.C.ok(listgood,fileDtp{ievent}));
        spar = spar(:,listgood,:);
        % figure;plot(spar(:,:,1))
        Nc = 1; %force le nombre de component a 2
        opt(1) = 1e-6;
        opt(2) = 0;    %initialisation 10 all methods
        opt(3) = 0;     %plotting
        opt(5) = 0;     %how often to show fit.
        opt(6) = 0;     %Max num iterations
        const = [0 0 0]; % constraints 0- nothing; 1-orthog; 2- nonneg; 3- unimod
        Oldload{2} = []; % rand(21,3);
        fixMode = [0 0 0]; %force 2 WAVELENGHT EQUAL %no fixe mode
        weights = []; %Mean for example put one at good time point and zero at noisy one
        
        [Factors,it,err,corcondia] = parafac(spar,Nc,opt,const,Oldload,fixMode,weights);
        
        PARAFAC.Factors = Factors;
        factorWavelength = Factors{3};
        
        
        A = Factors{1};
        B = Factors{2};
        C = Factors{3};
        % Try to keep HbO positive and HbR negative to respect HRF expected
        % behavior Parafac decomposition could not identify the sign
        % correctly as -1 * -1 = 1 try to keep meaning full sign to
        % avoid weird averaging.
        %Case A Sign A=1 B=1 C=1 correclty assign
        %Case B Sign A=-1 B=-1 C=1
        %Case C Sign A=1 B=-1 C=-1
        %Case D Sign A=-1  B=1 C=-1
        %on only one component
        
        if C(1) > C(2)
            if mean(A)>0
                cas = 'CaseA';
            end
            if mean(A)<=0
                cas = 'CaseB';
            end
        elseif C(1) < C(2)
            if mean(A)> 0
                cas = 'CaseC';
            elseif  mean(A)< 0
                cas = 'CaseD';
            end
        end
        
        switch cas
            case 'CaseA'
            case 'CaseB'
                A = -A;
                B = -B;
            case 'CaseC'
                B = -B;
                C = -C;
            case 'CaseD'
                A = -A;
                C = -C;
        end
        
        %         figure;subplot(2,2,1);plot(A);subplot(2,2,2);plot(B);subplot(2,2,3);plot(C);
        %         title(cas)
        
        ComponentToKeep=1;
        
        Ac = A(:,ComponentToKeep); Bc = B(:,ComponentToKeep); Cc = C(:,ComponentToKeep);
        [Xm]=nmodel(({Ac,Bc,Cc}));
        
        
        %check the get PARAFAC
        try
            load(fullfile(NIRSDtp{ievent},'SelectedFactors.mat'));
            newfile = 0;
        catch
            clear PARCOMP
            %donot exist create the stucture
            PARCOMP.file= fileDtp{ievent};
            PARCOMP.filestr =  fprintf('Bloc%03.0f,',fileDtp{ievent});
            PARCOMP.module  =numel(NIRS.Dt.fir.pp);
            PARCOMP.modulestr = NIRS.Dt.fir.pp(end).pre;
            PARCOMP.listgood =  listgood;
            PARCOMP.indt = indt; %indice de temps.
            PARCOMP.data =spar;
            PARCOMP.Xm = Xm;
            PARCOMP.FacA = A;
            PARCOMP.FacB = B;
            PARCOMP.FacC = C;
            PARCOMP.ComponentToKeep = ComponentToKeep ;
            labelid  = labelDtp{ievent} ;
            PARCOMP.label= [labelid,'PARFAC', sprintf('%03.0f',size(PARCOMP,2))];
            PARCOMP.type = 'PARAFAC';
            newfile = 1;
            FacSpatial = Factors{2};
            PARCOMP.topo = B(:,ComponentToKeep);
        end
        if newfile == 0
            id = numel(PARCOMP);
            PARCOMP(id+1).file= fileDtp{ievent};
            PARCOMP(id+1).filestr =  fprintf('Bloc%03.0f,',fileDtp{ievent});
            PARCOMP(id+1).module   =numel(NIRS.Dt.fir.pp);
            PARCOMP(id+1).modulestr  = NIRS.Dt.fir.pp(end).pre;
            PARCOMP(id+1).listgood =  listgood;
            PARCOMP(id+1).indt = indt; %indice de temps.
            PARCOMP(id+1).data = spar;
            PARCOMP(id+1).Xm = Xm;
            PARCOMP(id+1).FacA = A;
            PARCOMP(id+1).FacB = B;
            PARCOMP(id+1).FacC = C;
            PARCOMP(id+1).ComponentToKeep = ComponentToKeep;
            labelid  = labelDtp{ievent} ;
            PARCOMP(id+1).label= [labelid,'PARAFAC' sprintf('%03.0f',size(PARCOMP,2))];
            PARCOMP(id+1).type = 'PARAFAC';
            FacSpatial = Factors{2};
            PARCOMP(id+1).topo = B(:,ComponentToKeep);
        end
        
        save(fullfile(NIRSDtp{ievent},'SelectedFactors.mat'),'PARCOMP');
        
    end
elseif isfield(job.c_extractcomponent,'b_extractcomponent_AVG')
    
    %%
    % 
    % * ITEM1
    % * ITEM2
    % 
    [~,~,ext] =fileparts(job.c_extractcomponent.b_extractcomponent_AVG.f_component_AVGlist{1});
    if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
        [pathstr, name, ext]= fileparts(job.c_extractcomponent.b_extractcomponent_AVG.f_component_AVGlist{1});
         try
        [data, text, rawData] = xlsread(job.c_extractcomponent.b_extractcomponent_AVG.f_component_AVGlist{1});
        catch
        [data, text, rawData] = readtxtfile_asxlsread(job.c_extractcomponent.b_extractcomponent_AVG.f_component_AVGlist{1});
        end 

    elseif strcmp(ext,'.txt')
        [pathstr, name, ext]= fileparts(job.c_extractcomponent.b_extractcomponent_AVG.f_component_AVGlist{1});
        [data, text, rawData] = readtxtfile_asxlsread(job.c_extractcomponent.b_extractcomponent_AVG.f_component_AVGlist{1});
    end
    
    for icol=1:size(rawData,2)
        if strcmp(upper(deblank(rawData{1,icol})),deblank(upper('NIRS.mat folder')))
            id.NIRSDtp = icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('File'))
            id.fileDtp =  icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('tStart'))
            id.startDtp =  icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('tStop'))
            id.stopDtp  =  icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('Label'))
            id.labelDtp =  icol;
        elseif  strcmp(upper(deblank(rawData{1,icol})),upper('tStartAVG'))
            id.startwDtp =  icol;
        elseif  strcmp(upper(deblank(rawData{1,icol})),upper('tStopAVG'))
            id.stopwDtp =  icol;
        elseif strcmp(upper(deblank(rawData{1,icol})),upper('ZoneDisplay'))
            id.ZoneDisplay =  icol;
        end
    end
    NIRSDtp = rawData(2:end,id.NIRSDtp );
    fileDtp = rawData(2:end,2);
    startDtp = rawData(2:end,id.startDtp);
    stopDtp = rawData(2:end,id.stopDtp);
    startwDtp = rawData(2:end,id.startwDtp);
    stopwDtp = rawData(2:end,id.stopwDtp);
    labelDtp = rawData(2:end,id.labelDtp);
    if isfield(id,'ZoneDisplay')
        zoneDtp = rawData(2:end,id.ZoneDisplay);
    else
        zoneDtp = labelDtp;
    end
    
    %Load data and aux for regression
    for ievent=1:size(NIRSDtp,1)
        NIRSDtp{ievent};
        NIRSmat = fullfile(NIRSDtp{ievent},'NIRS.mat');
        try
            load(NIRSmat);
            NC = NIRS.Cf.H.C.N;
            fDtp = NIRS.Dt.fir.pp(end).p;
            try
                 d = fopen_NIR(fDtp{fileDtp{ievent}},NC)';
            catch
                jobfolderadjustment.NIRSmat = {fullfile(NIRSDtp{ievent},'NIRS.mat')};
                jobfolderadjustment.c_MultimodalPath.b_MultimodalPath_no = struct([]);
                outfolderadjustment = nirs_run_NIRSmatdiradjust(jobfolderadjustment)
                load(fullfile(NIRSDtp{ievent},'NIRS.mat'));  
                fDtp = NIRS.Dt.fir.pp(end).p;
                d = fopen_NIR(fDtp{fileDtp{ievent}},NC)';  
                disp('Folder adjustement apply')
            end
            listgood = 1:NC;
            tHRF = 1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:size(d,1)*1/NIRS.Cf.dev.fs;
            fsNIRS = NIRS.Cf.dev.fs;
            tstart = find(tHRF<=startDtp{ievent});
            if isempty(tstart)
                tstart = 1;
            end
            tstop = find(tHRF<= stopDtp{ievent});
            indt = [tstart(end),tstop(end)];%Time indice
            intensnorm = d(indt(1):indt(end),:);
            tstartw = find(tHRF<=startwDtp{ievent});
            tstopw = find(tHRF<= stopwDtp{ievent});
            try 
                load(zoneDtp{ievent},'-mat');
                %
            catch
                try
                load(fullfile(pathstr,[zoneDtp{ievent},'.zone']),'-mat');
                catch
                zone.plotLst{1} = 1;
                end
            end
            
            Xm = zeros(size(intensnorm,1), numel(zone.plotLst));
            AVG = zeros(NC,1 );
            for izone = 1:numel(zone.plotLst)
                plotLst = zone.plotLst{izone};
                Xm(:,izone) = nanmean(intensnorm(:,plotLst),2);
                AVG(plotLst,1) = nanmean(nanmean(intensnorm(:,plotLst),2));
            end
            if isfield(job.c_extractcomponent.b_extractcomponent_AVG.c_extractAVGlist_autoexport,'b_extractAVGlist_autoexport_yes')
                pathoutlist = job.c_extractcomponent.b_extractcomponent_AVG.c_extractAVGlist_autoexport.b_extractAVGlist_autoexport_yes.f_extractAVGlist_autoexport_yes{1};
                %review to improve channel list
                A = nanmean(d(tstartw(end) :tstopw(end),1:NC/2))';
                zonelist = []
                try
                save(fullfile(pathoutlist,['TopoHbO',labelDtp{ievent},'event',sprintf('%03.0f',ievent),'.mat']),'A','zonelist' );
                catch
                    disp(['Error ',fullfile(pathoutlist,['TopoHbO',labelDtp{ievent},'event',sprintf('%03.0f',ievent),'.mat']),' could not be saved. Please verify that the output folder exists' ]);
                end   
                disp(['Save :', fullfile(pathoutlist,['TopoHbO',labelDtp{ievent},'event',sprintf('%03.0f',ievent),'.mat'])])
                A = nanmean(d(tstartw(end) :tstopw(end),NC/2+1:end))';
                save(fullfile(pathoutlist,['TopoHbR',labelDtp{ievent},'event',sprintf('%03.0f',ievent),'.mat']),'A','zonelist' );
                disp(['Save :', fullfile(pathoutlist,['TopoHbR',labelDtp{ievent},'event',sprintf('%03.0f',ievent),'.mat'])])
            end
            
            %check the get PARAFAC
            try
                load(fullfile(NIRSDtp{ievent},'SelectedFactors.mat'));                
                newfile = 0;
            catch
                clear PARCOMP
                %donot exist create the stucture
                PARCOMP.file= fileDtp{ievent};
                PARCOMP.filestr =  fprintf('Bloc%03.0f,',fileDtp{ievent});
                PARCOMP.module  =numel(NIRS.Dt.fir.pp);
                PARCOMP.modulestr = NIRS.Dt.fir.pp(end).pre;
                PARCOMP.listgood =  listgood;
                PARCOMP.indt = indt; %indice de temps.
                PARCOMP.data =intensnorm;
                PARCOMP.Xm = Xm;
                labelid  = labelDtp{ievent} ;
                PARCOMP.label= [labelid,'AVG', sprintf('%03.0f',size(PARCOMP,2))];
                PARCOMP.type = 'AVG';
                newfile = 1;
                PARCOMP.topo = AVG;
            end
            if newfile == 0
                id = numel(PARCOMP);
                PARCOMP(id+1).file= fileDtp{ievent};
                PARCOMP(id+1).filestr =  fprintf('Bloc%03.0f,',fileDtp{ievent});
                PARCOMP(id+1).module  =numel(NIRS.Dt.fir.pp);
                PARCOMP(id+1).modulestr = NIRS.Dt.fir.pp(end).pre;
                PARCOMP(id+1).listgood =  listgood;
                PARCOMP(id+1).indt = indt; %indice de temps.
                PARCOMP(id+1).data =intensnorm;
                PARCOMP(id+1).Xm = Xm;
                labelid  = labelDtp{ievent} ;
                PARCOMP(id+1).label= [labelid,'AVG', sprintf('%03.0f',size(PARCOMP,2))];
                PARCOMP(id+1).type = 'AVG';
                newfile = 1;
                PARCOMP(id+1).topo = AVG;
            end
            
            save(fullfile(NIRSDtp{ievent},'SelectedFactors.mat'),'PARCOMP');
        catch
            disp(['Error unable to AVG on ' , NIRSmat])
        end  
        disp('Warning function to be revised')
    end
    
end
out.NIRSmat = {NIRSmat};


function  [y, h] = resample( x, p, q, N, bta )
%RESAMPLE  Change the sampling rate of a signal.
%   Y = RESAMPLE(X,P,Q) resamples the sequence in vector X at P/Q times
%   the original sample rate using a polyphase implementation.  Y is P/Q 
%   times the length of X (or the ceiling of this if P/Q is not an integer).  
%   P and Q must be positive integers.
%
%   RESAMPLE applies an anti-aliasing (lowpass) FIR filter to X during the 
%   resampling process, and compensates for the filter's delay.  The filter 
%   is designed using FIRLS.  RESAMPLE provides an easy-to-use alternative
%   to UPFIRDN, relieving the user of the need to supply a filter or
%   compensate for the signal delay introduced by filtering.
%
%   In its filtering process, RESAMPLE assumes the samples at times before
%   and after the given samples in X are equal to zero. Thus large
%   deviations from zero at the end points of the sequence X can cause
%   inaccuracies in Y at its end points.
%
%   Y = RESAMPLE(X,P,Q,N) uses a weighted sum of 2*N*max(1,Q/P) samples of X 
%   to compute each sample of Y.  The length of the FIR filter RESAMPLE applies
%   is proportional to N; by increasing N you will get better accuracy at the 
%   expense of a longer computation time.  If you don't specify N, RESAMPLE uses
%   N = 10 by default.  If you let N = 0, RESAMPLE performs a nearest
%   neighbor interpolation; that is, the output Y(n) is X(round((n-1)*Q/P)+1)
%   ( Y(n) = 0 if round((n-1)*Q/P)+1 > length(X) ).
%
%   Y = RESAMPLE(X,P,Q,N,BTA) uses BTA as the BETA design parameter for the 
%   Kaiser window used to design the filter.  RESAMPLE uses BTA = 5 if
%   you don't specify a value.
%
%   Y = RESAMPLE(X,P,Q,B) uses B to filter X (after upsampling) if B is a 
%   vector of filter coefficients.  RESAMPLE assumes B has odd length and
%   linear phase when compensating for the filter's delay; for even length 
%   filters, the delay is overcompensated by 1/2 sample.  For non-linear 
%   phase filters consider using UPFIRDN.
%
%   [Y,B] = RESAMPLE(X,P,Q,...) returns in B the coefficients of the filter
%   applied to X during the resampling process (after upsampling).
%
%   If X is a matrix, RESAMPLE resamples the columns of X.
%
%   See also UPFIRDN, INTERP, DECIMATE, FIRLS, KAISER, INTFILT,
%   MFILT/FIRSRC in the Filter Design Toolbox.

%   NOTE-1: digital anti-alias filter is desiged via windowing

%   Author(s): James McClellan, 6-11-93
%              Modified to use upfirdn, T. Krauss, 2-27-96
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.9.4.7 $  $Date: 2007/12/14 15:05:59 $

if nargin < 5,  bta = 5;  end   %--- design parameter for Kaiser window LPF
if nargin < 4,   N = 10;   end
if abs(round(p))~=p || p==0, error(generatemsgid('MustBePosInteger'),'P must be a positive integer.'), end
if abs(round(q))~=q || q==0, error(generatemsgid('MustBePosInteger'),'Q must be a positive integer.'), end

[p,q] = rat( p/q, 1e-12 );  %--- reduce to lowest terms 
   % (usually exact, sometimes not; loses at most 1 second every 10^12 seconds)
if (p==1) && (q==1)
    y = x; 
    h = 1;
    return
end
pqmax = max(p,q);
if length(N)>1      % use input filter
   L = length(N);
   h = N;
else                % design filter
   if( N>0 )
      fc = 1/2/pqmax;
      L = 2*N*pqmax + 1;
      h = p*firls( L-1, [0 2*fc 2*fc 1], [1 1 0 0]).*kaiser(L,bta)' ;
      % h = p*fir1( L-1, 2*fc, kaiser(L,bta)) ;
   else
      L = p;
      h = ones(1,p);
   end
end

Lhalf = (L-1)/2;
isvect = any(size(x)==1);
if isvect
    Lx = length(x);
else
    Lx = size(x, 1);
end

% Need to delay output so that downsampling by q hits center tap of filter.
nz = floor(q-mod(Lhalf,q));
z = zeros(1,nz);
h = [z h(:).'];  % ensure that h is a row vector.
Lhalf = Lhalf + nz;

% Number of samples removed from beginning of output sequence 
% to compensate for delay of linear phase filter:
delay = floor(ceil(Lhalf)/q);

% Need to zero-pad so output length is exactly ceil(Lx*p/q).
nz1 = 0;
while ceil( ((Lx-1)*p+length(h)+nz1 )/q ) - delay < ceil(Lx*p/q)
    nz1 = nz1+1;
end
h = [h zeros(1,nz1)];

% ----  HERE'S THE CALL TO UPFIRDN  ----------------------------
y = upfirdn(x,h,p,q);

% Get rid of trailing and leading data so input and output signals line up
% temporally:
Ly = ceil(Lx*p/q);  % output length
% Ly = floor((Lx-1)*p/q+1);  <-- alternately, to prevent "running-off" the
%                                data (extrapolation)
if isvect
    y(1:delay) = [];
    y(Ly+1:end) = [];
else
    y(1:delay,:) = [];
    y(Ly+1:end,:) = [];
end

h([1:nz (end-nz1+1):end]) = [];  % get rid of leading and trailing zeros 
                                 % in case filter is output
function [h,a]=firls(N,F,M,W,ftype)
% FIRLS Linear-phase FIR filter design using least-squares error minimization.
%   B=FIRLS(N,F,A) returns a length N+1 linear phase (real, symmetric
%   coefficients) FIR filter which has the best approximation to the
%   desired frequency response described by F and A in the least squares
%   sense. F is a vector of frequency band edges in pairs, in ascending
%   order between 0 and 1. 1 corresponds to the Nyquist frequency or half
%   the sampling frequency. A is a real vector the same size as F
%   which specifies the desired amplitude of the frequency response of the
%   resultant filter B. The desired response is the line connecting the
%   points (F(k),A(k)) and (F(k+1),A(k+1)) for odd k; FIRLS treats the
%   bands between F(k+1) and F(k+2) for odd k as "transition bands" or
%   "don't care" regions. Thus the desired amplitude is piecewise linear
%   with transition bands.  The integrated squared error is minimized.
%
%   For filters with a gain other than zero at Fs/2, e.g., highpass
%   and bandstop filters, N must be even.  Otherwise, N will be
%   incremented by one. Alternatively, you can use a trailing 'h' flag to
%   design a type 4 linear phase filter and avoid incrementing N.
%
%   B=FIRLS(N,F,A,W) uses the weights in W to weight the error. W has one
%   entry per band (so it is half the length of F and A) which tells
%   FIRLS how much emphasis to put on minimizing the integral squared error
%   in each band relative to the other bands.
%
%   B=FIRLS(N,F,A,'Hilbert') and B=FIRLS(N,F,A,W,'Hilbert') design filters
%   that have odd symmetry, that is, B(k) = -B(N+2-k) for k = 1, ..., N+1.
%   A special case is a Hilbert transformer which has an approx. amplitude
%   of 1 across the entire band, e.g. B=FIRLS(30,[.1 .9],[1 1],'Hilbert').
%
%   B=FIRLS(N,F,A,'differentiator') and B=FIRLS(N,F,A,W,'differentiator')
%   also design filters with odd symmetry, but with a special weighting
%   scheme for non-zero amplitude bands. The weight is assumed to be equal
%   to the inverse of frequency, squared, times the weight W. Thus the
%   filter has a much better fit at low frequency than at high frequency.
%   This designs FIR differentiators.
%
%   % Example of a length 31 lowpass filter.
%   h=firls(30,[0 .1 .2 .5]*2,[1 1 0 0]);
%   fvtool(h);
%
%   % Example of a length 45 lowpass differentiator.
%   h=firls(44,[0 .3 .4 1],[0 .2 0 0],'differentiator');
%   fvtool(h);
%
%   % Example of a length 26 type 4 highpass filter.
%   h=firls(25,[0 .4 .5 1],[0 0 1 1],'h');
%   fvtool(h);
%
%   See also FIRPM, FIR1, FIR2, FREQZ and FILTER.

%       Author(s): T. Krauss
%   History: 10-18-91, original version
%            3-30-93, updated
%            9-1-95, optimize adjacent band case
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.11.4.7 $  $Date: 2009/05/23 08:13:29 $

% check number of arguments, set up defaults.
error(nargchk(3,5,nargin,'struct'));

if (max(F)>1) || (min(F)<0)
    error(generatemsgid('InvalidRange'),'Frequencies in F must be in range [0,1].')
end
if (rem(length(F),2)~=0)
    error(generatemsgid('InvalidDimensions'),'F must have even length.');
end
if (length(F) ~= length(M))
    error(generatemsgid('InvalidDimensions'),'F and A must be equal lengths.');
end
if (nargin==3),
    W = ones(length(F)/2,1);
    ftype = '';
end
if (nargin==4),
    if ischar(W),
        ftype = W; W = ones(length(F)/2,1);
    else
        ftype = '';
    end
end
if (nargin==5),
    if isempty(W),
        W = ones(length(F)/2,1);
    end
end
if isempty(ftype)
    ftype = 0;  differ = 0;
else
    ftype = lower(ftype);
    if strcmpi(ftype,'h') || strcmpi(ftype,'hilbert')
        ftype = 1;  differ = 0;
    elseif strcmpi(ftype,'d') || strcmpi(ftype,'differentiator')
        ftype = 1;  differ = 1;
    else
        error(generatemsgid('InvalidEnum'),'Requires symmetry to be ''Hilbert'' or ''differentiator''.')
    end
end

% Check for valid filter length
[N,msg1,msg2] = firchk(N,F(end),M,ftype);
if ~isempty(msg1)
    error(generatemsgid('InvalidFilterOrder'),msg1);
end

if ~isempty(msg2),
    msg2 = sprintf([msg2,'\r',...
        '\nAlternatively, you can pass a trailing ''h'' argument,\r',...
        'as in firls(N,F,A,W,''h''), to design a type 4 linear phase filter.']);
    warning(generatemsgid('OrderIncreasedByOne'),msg2); 
end



N = N+1;                   % filter length
F=F(:)/2;  M=M(:);  W=sqrt(W(:));  % make these guys columns
dF = diff(F);

if (length(F) ~= length(W)*2)
    error(generatemsgid('InvalidDimensions'),'There should be one weight per band.');
end;
if any(dF<0),
    error(generatemsgid('InvalidFreqVec'),'Frequencies in F must be nondecreasing.')
end

% Fix for 67187
if all(dF(2:2:length(dF)-1)==0) && length(dF) > 1,
    fullband = 1;
else
    fullband = 0;
end
if all((W-W(1))==0)
    constant_weights = 1;
else
    constant_weights = 0;
end

L=(N-1)/2;

Nodd = rem(N,2);

if (ftype == 0),  % Type I and Type II linear phase FIR
    % basis vectors are cos(2*pi*m*f) (see m below)
    if ~Nodd
        m=(0:L)+.5;   % type II
    else
        m=(0:L);      % type I
    end
    k=m';
    need_matrix = (~fullband) || (~constant_weights);
    if need_matrix
        I1=k(:,ones(size(m)))+m(ones(size(k)),:);    % entries are m + k
        I2=k(:,ones(size(m)))-m(ones(size(k)),:);    % entries are m - k
        G=zeros(size(I1));
    end

    if Nodd
        k=k(2:length(k));
        b0=0;       %  first entry must be handled separately (where k(1)=0)
    end;
    b=zeros(size(k));
    for s=1:2:length(F),
        m=(M(s+1)-M(s))/(F(s+1)-F(s));    %  slope
        b1=M(s)-m*F(s);                   %  y-intercept
        if Nodd
            b0 = b0 + (b1*(F(s+1)-F(s)) + m/2*(F(s+1)*F(s+1)-F(s)*F(s)))...
                * abs(W((s+1)/2)^2) ;
        end
        b = b+(m/(4*pi*pi)*(cos(2*pi*k*F(s+1))-cos(2*pi*k*F(s)))./(k.*k))...
            * abs(W((s+1)/2)^2);
        b = b + (F(s+1)*(m*F(s+1)+b1)*sinc(2*k*F(s+1)) ...
            - F(s)*(m*F(s)+b1)*sinc(2*k*F(s))) ...
            * abs(W((s+1)/2)^2);
        if need_matrix
            G = G + (.5*F(s+1)*(sinc(2*I1*F(s+1))+sinc(2*I2*F(s+1))) ...
                - .5*F(s)*(sinc(2*I1*F(s))+sinc(2*I2*F(s))) ) ...
                * abs(W((s+1)/2)^2);
        end
    end;
    if Nodd
        b=[b0; b];
    end;

    if need_matrix
        a=G\b;
    else
        a=(W(1)^2)*4*b;
        if Nodd
            a(1) = a(1)/2;
        end
    end
    if Nodd
        h=[a(L+1:-1:2)/2; a(1); a(2:L+1)/2].';
    else
        h=.5*[flipud(a); a].';
    end;
elseif (ftype == 1),  % Type III and Type IV linear phase FIR
    %  basis vectors are sin(2*pi*m*f) (see m below)
    if (differ),      % weight non-zero bands with 1/f^2
        do_weight = ( abs(M(1:2:length(M))) +  abs(M(2:2:length(M))) ) > 0;
    else
        do_weight = zeros(size(F));
    end

    if Nodd
        m=(1:L);      % type III
    else
        m=(0:L)+.5;   % type IV
    end;
    k=m';
    b=zeros(size(k));

    need_matrix = (~fullband) || (any(do_weight)) || (~constant_weights);
    if need_matrix
        I1=k(:,ones(size(m)))+m(ones(size(k)),:);    % entries are m + k
        I2=k(:,ones(size(m)))-m(ones(size(k)),:);    % entries are m - k
        G=zeros(size(I1));
    end

    i = sqrt(-1);
    for s=1:2:length(F),
        if (do_weight((s+1)/2)),      % weight bands with 1/f^2
            if F(s) == 0, F(s) = 1e-5; end     % avoid singularities
            m=(M(s+1)-M(s))/(F(s+1)-F(s));
            b1=M(s)-m*F(s);
            snint1 = sineint(2*pi*k*F(s+1)) - sineint(2*pi*k*F(s));
            %snint1 = (-1/2/i)*(expint(i*2*pi*k*F(s+1)) ...
            %    -expint(-i*2*pi*k*F(s+1)) -expint(i*2*pi*k*F(s)) ...
            %    +expint(-i*2*pi*k*F(s)) );
            % csint1 = cosint(2*pi*k*F(s+1)) - cosint(2*pi*k*F(s)) ;
            csint1 = (-1/2)*(expint(i*2*pi*k*F(s+1))+expint(-i*2*pi*k*F(s+1))...
                -expint(i*2*pi*k*F(s))  -expint(-i*2*pi*k*F(s)) );
            b=b + ( m*snint1 ...
                + b1*2*pi*k.*( -sinc(2*k*F(s+1)) + sinc(2*k*F(s)) + csint1 ))...
                * abs(W((s+1)/2)^2);
            snint1 = sineint(2*pi*F(s+1)*(-I2));
            snint2 = sineint(2*pi*F(s+1)*I1);
            snint3 = sineint(2*pi*F(s)*(-I2));
            snint4 = sineint(2*pi*F(s)*I1);
            G = G - ( ( -1/2*( cos(2*pi*F(s+1)*(-I2))/F(s+1)  ...
                - 2*snint1*pi.*I2 ...
                - cos(2*pi*F(s+1)*I1)/F(s+1) ...
                - 2*snint2*pi.*I1 )) ...
                - ( -1/2*( cos(2*pi*F(s)*(-I2))/F(s)  ...
                - 2*snint3*pi.*I2 ...
                - cos(2*pi*F(s)*I1)/F(s) ...
                - 2*snint4*pi.*I1) ) ) ...
                * abs(W((s+1)/2)^2);
        else      % use usual weights
            m=(M(s+1)-M(s))/(F(s+1)-F(s));
            b1=M(s)-m*F(s);
            b=b+(m/(4*pi*pi)*(sin(2*pi*k*F(s+1))-sin(2*pi*k*F(s)))./(k.*k))...
                * abs(W((s+1)/2)^2) ;
            b = b + (((m*F(s)+b1)*cos(2*pi*k*F(s)) - ...
                (m*F(s+1)+b1)*cos(2*pi*k*F(s+1)))./(2*pi*k)) ...
                * abs(W((s+1)/2)^2) ;
            if need_matrix
                G = G + (.5*F(s+1)*(sinc(2*I1*F(s+1))-sinc(2*I2*F(s+1))) ...
                    - .5*F(s)*(sinc(2*I1*F(s))-sinc(2*I2*F(s)))) * ...
                    abs(W((s+1)/2)^2);
            end
        end;
    end

    if need_matrix
        a=G\b;
    else
        a=-4*b*(W(1)^2);
    end
    if Nodd
        h=.5*[flipud(a); 0; -a].';
    else
        h=.5*[flipud(a); -a].';
    end
    if differ, h=-h; end
end

if nargout > 1
    a = 1;
end

%----------------------------------------------------------------------------
function y = sineint(x)
% SINEINT (a.k.a. SININT)   Numerical Sine Integral
%   Used by FIRLS in the Signal Processing Toolbox.
%   Untested for complex or imaginary inputs.
%
%   See also SININT in the Symbolic Toolbox.

%   Was Revision: 1.5, Date: 1996/03/15 20:55:51

i1 = find(real(x)<0);   % this equation is not valid if x is in the
% left-hand plane of the complex plane.
% use relation Si(-z) = -Si(z) in this case (Eq 5.2.19, Abramowitz
%  & Stegun).
x(i1) = -x(i1);
y = zeros(size(x));
ind = find(x);
% equation 5.2.21 Abramowitz & Stegun
%  y(ind) = (1/(2*i))*(expint(i*x(ind)) - expint(-i*x(ind))) + pi/2;
y(ind) = imag(expint(i*x(ind))) + pi/2;
y(i1) = -y(i1);
function [n,msg1,msg2] = firchk(n,Fend,a,exception)
%FIRCHK   Check if specified filter order is valid.
%   FIRCHK(N,Fend,A) checks if the specified order N is valid given the
%   final frequency point Fend and the desired magnitude response vector A.
%   Type 2 linear phase FIR filters (symmetric, odd order) must have a
%   desired magnitude response vector that ends in zero if Fend = 1.  This
%   is because type 2 filters necessarily have a zero at w = pi.
%
%   If the order is not valid, a warning is given and the order
%   of the filter is incremented by one.
%
%   If A is a scalar (as when called from fircls1), A = 0 is
%   interpreted as lowpass and A = 1 is interpreted as highpass.
%
%   FIRCHK(N,Fend,A,EXCEPTION) will not warn or increase the order
%   if EXCEPTION = 1.  Examples of EXCEPTIONS are type 4 filters
%   (such as differentiators or hilbert transformers) or non-linear
%   phase filters (such as minimum and maximum phase filters).

%   Author : R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.7.4.5 $  $Date: 2007/12/14 15:15:06 $

error(nargchk(3,4,nargin,'struct'));

if nargin == 3,
    exception = false;
end

msg1 = '';
msg2 = '';
oddord = false; % Flag, initially we assume even order

if isempty(n) || length(n) > 1 || ~isnumeric(n) || ~isreal(n) || n~=round(n) || n<=0,
    msg1 = 'Filter order must be a real, positive integer.';
    return
end

if rem(n,2) == 1,
    oddord = true; % Overwrite flag
end
 
if (a(end) ~= 0) && Fend == 1 && oddord && ~exception,
    str = ['Odd order symmetric FIR filters must have a gain of zero \n'...
     'at the Nyquist frequency. The order is being increased by one.'];
    msg2 = sprintf(str);
    n = n+1;
end

function y=sinc(x)
%SINC Sin(pi*x)/(pi*x) function.
%   SINC(X) returns a matrix whose elements are the sinc of the elements 
%   of X, i.e.
%        y = sin(pi*x)/(pi*x)    if x ~= 0
%          = 1                   if x == 0
%   where x is an element of the input matrix and y is the resultant
%   output element.
%
%   % Example of a sinc function for a linearly spaced vector:
%   t = linspace(-5,5);
%   y = sinc(t);
%   plot(t,y);
%   xlabel('Time (sec)');ylabel('Amplitude'); title('Sinc Function')
%
%   See also SQUARE, SIN, COS, CHIRP, DIRIC, GAUSPULS, PULSTRAN, RECTPULS,
%   and TRIPULS.

%   Author(s): T. Krauss, 1-14-93
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.7.4.1 $  $Date: 2004/08/10 02:11:27 $

i=find(x==0);                                                              
x(i)= 1;      % From LS: don't need this is /0 warning is off                           
y = sin(pi*x)./(pi*x);                                                     
y(i) = 1;   
function w = kaiser(n_est,bta)
%KAISER Kaiser window.
%   W = KAISER(N) returns an N-point Kaiser window in the column vector W.
% 
%   W = KAISER(N,BTA) returns the BETA-valued N-point Kaiser window.
%       If omitted, BTA is set to 0.500.
%
%   See also CHEBWIN, GAUSSWIN, TUKEYWIN, WINDOW.

%   Author(s): L. Shure, 3-4-87
%   Copyright 1988-2005 The MathWorks, Inc.
%   $Revision: 1.17.4.4 $  $Date: 2007/12/14 15:05:16 $

error(nargchk(1,2,nargin,'struct'));

% Default value for the BETA parameter.
if nargin < 2 || isempty(bta), 
    bta = 0.500;
end

[nn,w,trivialwin] = check_order(n_est);
if trivialwin, return, end;

nw = round(nn);
bes = abs(besseli(0,bta));
odd = rem(nw,2);
xind = (nw-1)^2;
n = fix((nw+1)/2);
xi = (0:n-1) + .5*(1-odd);
xi = 4*xi.^2;
w = besseli(0,bta*sqrt(1-xi/xind))/bes;
w = abs([w(n:-1:odd+1) w])';

function [n_out, w, trivalwin] = check_order(n_in)
%CHECK_ORDER Checks the order passed to the window functions.
% [N,W,TRIVALWIN] = CHECK_ORDER(N_ESTIMATE) will round N_ESTIMATE to the
% nearest integer if it is not already an integer. In special cases (N is
% [], 0, or 1), TRIVALWIN will be set to flag that W has been modified.

%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.6.4.2 $  $Date: 2009/05/23 08:16:17 $

w = [];
trivalwin = 0;

if ~(isnumeric(n_in) & isfinite(n_in)),
    error(generatemsgid('InvalidOrder'),'The order N must be finite.');
end

% Special case of negative orders:
if n_in < 0,
   error(generatemsgid('InvalidOrder'),'Order cannot be less than zero.');
end

% Check if order is already an integer or empty
% If not, round to nearest integer.
if isempty(n_in) | n_in == floor(n_in),
   n_out = n_in;
else
   n_out = round(n_in);
   warning(generatemsgid('InvalidOrder'),'Rounding order to nearest integer.');
end

% Special cases:
if isempty(n_out) | n_out == 0,
   w = zeros(0,1);               % Empty matrix: 0-by-1
   trivalwin = 1; 
elseif n_out == 1,
   w = 1;
   trivalwin = 1;   
end
% [EOF] kaiser.m
function Y = upfirdn(x,h,varargin)
%UPFIRDN  Upsample, apply a specified FIR filter, and downsample a signal.
%   UPFIRDN(X,H,P,Q) is a cascade of three systems applied to input signal X:
%         (1) Upsampling by P (zero insertion).  P defaults to 1 if not 
%             specified.
%         (2) FIR filtering with the filter specified by the impulse response 
%             given in H.
%         (3) Downsampling by Q (throwing away samples).  Q defaults to 1 if not 
%             specified.
%   UPFIRDN uses an efficient polyphase implementation.
%
%   Usually X and H are vectors, and the output is a (signal) vector. 
%   UPFIRDN permits matrix arguments under the following rules:
%   If X is a matrix and H is a vector, each column of X is filtered through H.
%   If X is a vector and H is a matrix, each column of H is used to filter X.
%   If X and H are both matrices with the same number of columns, then the i-th
%      column of H is used to filter the i-th column of X.
%
%   Specifically, these rules are carried out as follows.  Note that the length
%   of the output is Ly = ceil( ((Lx-1)*P + Lh)/Q ) where Lx = length(X) and 
%   Lh = length(H). 
%
%      Input Signal X    Input Filter H    Output Signal Y   Notes
%      -----------------------------------------------------------------
%   1) length Lx vector  length Lh vector  length Ly vector  Usual case.
%   2) Lx-by-Nx matrix   length Lh vector  Ly-by-Nx matrix   Each column of X
%                                                            is filtered by H.
%   3) length Lx vector  Lh-by-Nh matrix   Ly-by-Nh matrix   Each column of H is
%                                                            used to filter X.
%   4) Lx-by-N matrix    Lh-by-N matrix    Ly-by-N matrix    i-th column of H is
%                                                            used to filter i-th
%                                                            column of X.
%
%   For an easy-to-use alternative to UPFIRDN, which does not require you to 
%   supply a filter or compensate for the signal delay introduced by filtering,
%   use RESAMPLE.
%
%   EXAMPLE: Sample-rate conversion by a factor of 147/160. It is used to
%            % downconvert from 48kHz to 44.1kHz.
%            L = 147; M = 160;                   % Interpolation/decimation factors.
%            Lp = 24;                            % Filter length of each phase
%            N = Lp*L-1;                         % Filter Order
%            h = fir1(N,1/M,kaiser(N+1,7.8562));
%            h = L*h; % Passband gain = L
%            Fs = 48e3;                          % Original sampling frequency: 48kHz
%            n = 0:10239;                        % 10240 samples, 0.213 seconds long
%            x  = sin(2*pi*1e3/Fs*n);            % Original signal, sinusoid at 1kHz
%            y = upfirdn(x,h,L,M);               % 9408 samples, still 0.213 seconds
%
%            % Overlay original (48kHz) with resampled signal (44.1kHz) in red.
%            stem(n(1:49)/Fs,x(1:49)); hold on 
%            stem(n(1:45)/(Fs*L/M),y(12:56),'r','filled'); 
%            xlabel('Time (sec)');ylabel('Signal value');
%    
%   See also RESAMPLE, INTERP, DECIMATE, FIR1, INTFILT, MFILT/FIRSRC in the
%   Filter Design Toolbox.
  
%   Author(s): Paul Pacheco
%   Copyright 1988-2008 The MathWorks, Inc.
%   $Revision: 1.6.4.5 $  $Date: 2008/09/13 07:14:26 $

%   This M-file validates the inputs, sets defaults, and then calls the C MEX-file.

% Validate number of I/O args.
error(nargchk(2,4,nargin,'struct'));
error(nargoutchk(0,1,nargout,'struct'));

% Force to be a column if input is a vector
[mx,nx] = size(x);
if find([mx nx]==1),
  x = x(:);  % columnize it.
end
[Lx,nChans] = size(x);

% Force to be a column if filter is a vector
if find(size(h)==1),
    h = h(:);  % columnize it.
end
[Lh,hCols] = size(h);

% Validate input args and define defaults.
[p,q,errid,errmsg] = validateinput(x,h,varargin);
if ~isempty(errmsg), error(errid,errmsg); end

% Call the MEX-file
Y = upfirdnmex(x,h,p,q,Lx,Lh,hCols,nChans);

% Convert output to be a row vector (if x was a row and H is NOT a matrix)
if (mx==1) && (hCols == 1)
    Y = Y(:).';
end


%----------------------------------------------------------------------
function [p,q,errid,errmsg] = validateinput(x,h,opts)

% Default values
p = 1;
q = 1;
errid = '';
errmsg = '';

% Validate 1st two input args: signal and filter.
if isempty(x) || issparse(x) || ~isa(x,'double'),
    errid = generatemsgid('invalidInput');
    errmsg = 'The input signal X must be a double-precision vector.';
    return;
end
if isempty(h) || issparse(h) || ~isa(h,'double'),
    errid = generatemsgid('invalidFilter');
    errmsg = 'The filter H must be a double-precision vector.';
    return;
end

% The following check is for case 4 (as seen on the reference page), i.e., 
% x and h are matrices, check that they both have the same number of
% columns. 
nChans = size(x, 2);
hCols  = size(h, 2);
if (nChans > 1) && (hCols > 1) && (hCols ~= nChans),
    errid = generatemsgid('xNhSizemismatch');
    errmsg = 'Signal X and filter H must have the same number of columns.';
    return;
end

% Validate optional input args: upsample and downsample factors.
nopts = length(opts);
if (nopts >= 1),
    p = opts{1};
    if isempty(p) || ~isa(p,'double') || p<1 || ~isequal(round(p),p),
        errid = generatemsgid('invalidP');
        errmsg = 'The upsample factor P must be a positive, double-precision, integer.';
        return;

    elseif (nopts == 2),
        q = opts{2};
        if isempty(q) || ~isa(q,'double') || q<1 || ~isequal(round(q),q),
            errid = generatemsgid('invalidQ');
            errmsg = 'The downsample factor Q must be a positive, double-precision, integer.';
            return;
        end
    end
    if p*q > intmax('int32'),
        errid = generatemsgid('ProdPNQTooLarge');
        errmsg = ['The product of the downsample factor Q and the upsample factor P must be',...
            ' less than 2^31.'];
        return;
    end
end

% [EOF]

