function out = nirs_run_E_extractcomponent(job)
% Batch detection of the components
if isfield(job.c_extractcomponent,'b_extractcomponent_PCA')
    NIRSmatlst = job.c_extractcomponent.b_extractcomponent_PCA.NIRSmat;    
   
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
        indconsecutifthreshold = 1;     %en durée en sample
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
                    mrks = find(ind_dur_ch(:,3)==Idx);
                    ind = ind_dur_ch(mrks,1);
                    indf = ind + ind_dur_ch(mrks,2) - 1;
                    if ~isempty(ind)
                        try
                            for i = 1:numel(ind)
                                noise(ind(i):indf(i),Idx) = 1;
                            end
                        catch
                            msgbox('Noise reading problem')
                        end
                    end
                end
            end
            
            
            % ici group channel with the same noise latency
            figureon = m_extractcomponentfigure;
            ind = find((sum(noise,2)./size(noise,2))>nbchminimum );
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
        indconsecutifthreshold = 1;     %en durée en sample
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
            0.
            mrks = [];
            ind = [];
            noise =  logical(zeros(size(d)));
            time  = 1/fs:1/fs:1/fs*size(d,1)
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
                            msgbox('Noise reading problem')
                        end
                    end
                end
            end
            
            
            % ici group channel with the same noise latency
            figureon = m_extractcomponentfigure;
            ind = find((sum(noise,2)./size(noise,2))>nbchminimum );
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
                        spar = cat(3,spar(:,1:end/2),spar(:,end/2+1:end));
                        spartmp = spar(:,listgood,:);
                        listgood = find(MeasListActplotLst);
        
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
                                            [Factors,it,err,corcondia] = parafac(spartmp,Nc,opt,const,Oldload,fixMode,weights);
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
                        
                        rejeterhightimecourse =(mean(sumA(:)))< sumA';                        
                        rejectwavelengthdistance = mean(abs(distC(:))) > abs(distC)';                        
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
                            try;
                            Ac = A(:,ComponentToKeep); Bc = B(:,ComponentToKeep); Cc = C(:,ComponentToKeep);
                            [Xm]=nmodel(({Ac,Bc,Cc}));
                             catch;
                            1
                            end
                            data = cat(3,d(indt,1:end/2),d(indt,end/2+1:end));
                            data(:,listgood,:) = data(:,listgood,:)-Xm;
                            if 0 %figureon==1
                                subplot(3,2,3);plot(Factors{1});subplot(3,2,4);plot(Factors{2});subplot(3,2,5);plot(Factors{3});
                                subplot(3,2,6);plot(Ac);hold on
                            end
                            %             %Optimise the component to keep as lowest distance between wavelength.
                            %figure;plot(std(data(:,listgood,1)./mean(data(:,listgood,1))))
                            %decomposed until zscore 
                            
                            
                            if sum(std(data(:,listgood,1)./mean(data(:,listgood,1)))>0.1) & mean(data(:,listgood,1))> max(data(:))*0.05 %test when not much variation except for low intensity channel
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
     indconsecutifthreshold = 1;     %en durée en sample
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
                            msgbox('Noise reading problem')
                        end
                    end
                end
            end
            
            
            % ici group channel with the same noise latency
 
            ind = find((sum(noise,2)./size(noise,2))>nbchminimum );
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
            
            % add event start 
          eventbadstartstop = [idstart,idstop] ;
          temp =   permute(eventbadstartstop, [2,1])
          temp= [1;temp(:); size(noise,1)];
          eventint=  reshape(temp,2,numel(temp)/2)
          eventgoodstartstop=permute(eventint, [2,1])               
              
                 
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
                       zoneidentification = tmp(10:end)
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
            for izoneRegressor = 1:numel(ListRegressorZone) 
                %do the wavelength 1 
                 chlistRegressor = zone.plotLst{ListRegressorZone(izoneRegressor)};
                 idbad = find(badch( chlistRegressor)==0); %remove exclude channel from regressor
                 if ~isempty(idbad)
                      chlistRegressor(idbad) = [];
                      if isempty(chlistRegressor)
                          disp(['No good channel in the regressor zone please verify your zone ', zone.label{ListRegressorZone(izoneRegressor)}])
                          break
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
    [~,~,ext] =fileparts(job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{1});
    if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
        [data, text, rawData] = xlsread(job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{1});
        [dirxls,filexls,extxls] = fileparts(job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{1});
        id.Regressor = [];

    elseif strcmp(ext,'.txt')   
        [data, text, rawData] = readtxtfile_asxlsread(job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{1});
        [dirxls,filexls,extxls] = fileparts(job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{1});
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
   % trigDtp = rawData(2:end,4); non utilisé
    startDtp = rawData(2:end,id.startDtp);
    stopDtp = rawData(2:end,id.stopDtp);
    labelDtp = rawData(2:end,id.labelDtp);
    AUXid = rawData(2:end,id.Regressor );
    catch
        msgbox(['Please verify the GLM extract xls file : ', job.c_extractcomponent.b_extractcomponent_glm.f_extractcomponent_glmlist{1},...
            ' have the following column information : NIRS.mat folder, File, tStart, tStop, label and Xn regressors'])
    end
    %Load data and aux for regression
    for ievent=1:size(NIRSDtp,1)
        %try
            warning off            
            NIRSmat = fullfile(NIRSDtp{ievent},'NIRS.mat');
            load(NIRSmat);        
            Regressorlist=AUXid{ievent,:};
            disp(NIRSmat); 
            NC = NIRS.Cf.H.C.N;
            fDtp = NIRS.Dt.fir.pp(end).p;
            d1 = fopen_NIR(fDtp{fileDtp{ievent}},NC)';
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
            tmpGLM.spar = d1(tmpGLM.indt(1):tmpGLM.indt(end),:);
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
                  
                    tmpr=resample( tmp , p, q);
                  
                        
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
                end
            end
       
                if numel(Regressorlist{iReg})>4
                    checkzone = Regressorlist{iReg};
                    if strcmp(upper(checkzone(end-3:end)),'ZONE')
                        try 
                            load(Regressorlist{iReg} ,'-mat'); %try the fullfile correct 
                        catch
                            try
                            load(fullfile(dirxls, Regressorlist{iReg} ),'-mat'); %try in the excel folder
                            catch
                            disp(['Regressor zone :', fullfile(dirxls, Regressorlist{iReg} ),' could not be load'])
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
                    Xmean = nanmean(tmpGLM.spar(:,chlistRegressor),2) -  nanmean(nanmean(tmpGLM.spar(:,chlistRegressor),2)); %center to zero
                   % figure;plot(Xmean)
                     X = [Xtmp,    Xmean];
                    for ich = 1:numel(chlistApply)
                        idch = chlistApply(ich);
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
                    label = Regressorlist{iReg};
                   [filepath,name,ext] = fileparts(Regressorlist{iReg});
                   label =  [name,ext];
              
                else
                Xm = tmpGLM.AUX.data{idreg(iselected)}*beta(iselected,:);
                label = tmpGLM.AUX.label{idreg(iselected)};
                end
                try
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
                    PARCOMP.data = d1(tmpGLM.indt(1):tmpGLM.indt(end),tmpGLM.listgood) ;%data(:,listgood,:);
                    PARCOMP.Xm = Xm;
                    PARCOMP.ComponentToKeep = tmpGLM.selected;
                    PARCOMP.idreg = tmpGLM.idreg;
                    labelid  = labelDtp{ievent} ;
                    PARCOMP.label= [labelid,'GLM',label, sprintf('%03.0f',size(PARCOMP,2))];
                    disp([labelid,'GLM',label, sprintf('%03.0f',size(PARCOMP,2))]);
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
                    PARCOMP(id+1).data = d1(tmpGLM.indt(1):tmpGLM.indt(end),tmpGLM.listgood) ;%data(:,listgood,:);
                    PARCOMP(id+1).Xm = Xm;
                    PARCOMP(id+1).ComponentToKeep = tmpGLM.selected;
                    PARCOMP(id+1).idreg = tmpGLM.idreg;
                    labelid  = labelDtp{ievent} ;
                     PARCOMP(id+1).label= [labelid,'GLM',label, sprintf('%03.0f',size(PARCOMP,2))];
                    disp([labelid,'GLM',label , sprintf('%03.0f',size(PARCOMP,2))]);
                    PARCOMP(id+1).type = 'GLM';
                    PARCOMP(id+1).topo =  beta(iselected,:);
                    newfile = 1;
                end 
                save(fullfile(NIRSDtp{ievent},'SelectedFactors.mat'),'PARCOMP');
            end
                %disp(['Error unable to GLM on ' , NIRSmat])
    end    
   
elseif isfield(job.c_extractcomponent,'b_extractcomponent_PARAFAC')
    [~,~,ext] =fileparts(job.c_extractcomponent.b_extractcomponent_PARAFAC.f_component_PARAFAClist{1});
    if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
        [data, text, rawData] = xlsread(job.c_extractcomponent.b_extractcomponent_PARAFAC.f_component_PARAFAClist{1});
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
   % trigDtp = rawData(2:end,4); non utilisé
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
    [~,~,ext] =fileparts(job.c_extractcomponent.b_extractcomponent_AVG.f_component_AVGlist{1});
    if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
        [pathstr, name, ext]= fileparts(job.c_extractcomponent.b_extractcomponent_AVG.f_component_AVGlist{1});
        [data, text, rawData] = xlsread(job.c_extractcomponent.b_extractcomponent_AVG.f_component_AVGlist{1});
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
        d = fopen_NIR(fDtp{fileDtp{ievent}},NC)';
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
             load(fullfile(pathstr,[zoneDtp{ievent},'.zone']),'-mat')
        catch
            zone.plotLst{1} = 1;
        end
       
        Xm = zeros(size(intensnorm,1), numel(zone.plotLst));
        for izone = 1:numel(zone.plotLst)
            plotLst = zone.plotLst{izone}
            Xm(:,izone) = nanmean(intensnorm(:,plotLst),2);
            %  figure;plot(Xm)
        end
        AVG = nanmean(d(tstartw(end) :tstopw(end),listgood));
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
    end
end
out.NIRSmat = {NIRSmat};