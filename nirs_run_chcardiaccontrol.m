function out = nirs_run_chcardiaccontrol(job)
%Look for the correlation matrix ch per ch to ensure that the cardiac
%coherence is higher than a treshold, else it set channel as rejected
%channel NIRS.Cf.H.C.ok



NIRS = [];
idxls = 1;
SNRfft = 0;
idxlsSNR = 1;
for filenb=1:size(job.NIRSmat,1) %do it one by one for the associate name 
    load(job.NIRSmat{filenb,1});
    ML_new= [NIRS.Cf.H.C.id(2:3,:)',...
        ones(size(NIRS.Cf.H.C.id,2),1),...
        [ones(size(NIRS.Cf.H.C.id,2)/2,1);ones(size(NIRS.Cf.H.C.id,2)./2,1).*2]];
    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N;
outliertreshold = 5;
COHtr = job.i_Freq_crossspectrum;
pourcentagetr = job.i_minch_cardiac/100;
    for f=1:size(rDtp,1) %Loop over all files of a NIRS.mat
        d1 = fopen_NIR(rDtp{f,1},NC);       
        [pathstrtmp, nametmp, ext] = fileparts(NIRS.Dt.fir.pp(lst).p{f});
       fs = NIRS.Cf.dev.fs;                         % Sample frequency (Hz)
       tseg = 20;
       t = 0:1/fs:tseg;
       Bsize = numel(t);
        n = size(d1',1);
        p =floor(n/Bsize);
       nb_random_sample = 100;
       indbloc = 1:Bsize:n;
       maxval = n-Bsize;
       idstart=randi(maxval,nb_random_sample,1);
            %definition des blocs on arrete
            for ibloc = 1:numel(idstart)
                %pstart and pstop for each bloc
                Bloc(ibloc,:) = [idstart(ibloc), idstart(ibloc)+Bsize];
            end
            dat = d1(:,Bloc(1,1):Bloc(1,2));
            [y, f_fft]= fft_EEGseries(dat,fs); 
            yall = zeros(size(y,1),size(y,2),size(Bloc,1));
            removetrial = [];
            for ibloc = 1:size(Bloc,1)
                dat = d1(:,Bloc(ibloc,1):Bloc(ibloc,2));
                dat=dat - mean(dat,2)*ones(1,Bsize+1) ;
                if sum(isnan(dat(:)))
                    removetrial = [removetrial,ibloc];
                else
                    [y, f_fft]= fft_EEGseries(dat,fs);
                    yall(:,:,ibloc)=y;
                end
            end
            if ~isempty(removetrial)
                yall(:,:,removetrial)=[];
                totaltrialgood = size(yall,3);
            else
                totaltrialgood = size(yall,3);
            end
            %Load zone for display spectrum
            
          
            power = yall.*conj(yall)/n;   
            idHBO = 1:size(power,2)/2;
            idHBR = (size(power,2)/2+1): size(power,2);
             freq = job.i_Freq_cardiac;
            startF =sum(f_fft<=freq(1));
            stopF  =sum(f_fft<=freq(end));     
           minval = min(reshape(nanmean(log10(power(3:end,:,:))),numel(nanmean(log10(power(3:end,:,:)))),1))-3;
           maxval =  max(reshape(nanmean(log10(power(3:end,:,:))),numel(nanmean(log10(power(3:end,:,:)))),1))+3;
                      
            hfig = figure ;subplot(5,4,[1,2,3,4,5,6,7,8]);hold on
            for ich=1:size(power,2)
                plot(f_fft(3:end),nanmean(log10(power(3:end,ich,:)),3),'displayname',num2str(ich))
            end
           
            xlabel('Frequency (Hz)' ,'fontsize',12);
            ylabel('Power spectrum (dB)' ,'fontsize',12);
            set(gca,'fontsize',12);
           [val,id]=max(nanmean(log10(power(8:end,:,:)),3));
         
 
    
           minval = min(reshape(nanmean(log10(power(3:end,:,:))),numel(nanmean(log10(power(3:end,:,:)))),1))-3;
           maxval =  max(reshape(nanmean(log10(power(3:end,:,:))),numel(nanmean(log10(power(3:end,:,:)))),1))+3;
                    ylim([minval,maxval]);
                    plot([f_fft(startF),f_fft(startF)],[minval,maxval],'r')
                    plot([f_fft(stopF),f_fft(stopF)],[minval,maxval],'r')
            matcorr= nan(size(d1,1)/2,size(d1,1)/2,1);          
            matcorrHbR= nan(size(d1,1)/2,size(d1,1)/2,1);
          [val,id] =max( nanmean(nanmean(log10(power(startF:stopF,:,:)),3),2));
          plot([f_fft(startF+id-1),f_fft(startF+id-1)],[minval,maxval],'k')
         ylim([minval, maxval])
         
         
         %try max by eletrode iss seem to have two close cardiac peak estimation... 
        [val,ideach] =max( nanmean(log10(power(startF:stopF,:,:)),3));
        
            % Outlier zscore by ch and fr
            for ich = 1:size(power,2)%channel
                removetrial = zeros(size(power,3),1);
                for ifr=startF:stopF%freq
%                     nmean = mean(power(ifr,ich,:)))
%                     nstr = std(power(ifr,ich,:)))
                    list =find(abs(zscore(power(ifr,ich,:)))>outliertreshold  );
                    if ~isempty(list)
                        removetrial(list)=1;
                    end
                end
                idbad = find(removetrial);
                yall(startF:stopF,ich, idbad)=nan;
            end
            
         tmp = nanmean(log10(power(2:end,:,:)),3);
         tmp(find(isinf(tmp)))=0;            
            maxval = max(max( tmp))+5;
            minval =  min(min( tmp))-5;            
            power = yall.*conj(yall)/n;       
            plot([f_fft(startF),f_fft(startF)],[minval,maxval],'r')
            plot([f_fft(stopF),f_fft(stopF)],[minval,maxval],'r')
            title(['Power spectrum :', nametmp])
          
            avgfft = nanmean(nanmean(log10(power(startF:stopF,:,:)),3),2);
            interval  = startF:stopF;
            [val,id] = max(avgfft);
            tablepeak(f,1) = f_fft(interval(id));
           nbpoint = numel(find(f_fft<job.i_cardiacwidth));
           if job.i_cardiacwidth==0
               nbpoint = 0;
           end
          
             startF = (interval(id))-nbpoint;
            stopF = (interval(id))+nbpoint;
            plot([f_fft(startF),f_fft(startF)],[minval,maxval],'k')
            plot([f_fft(stopF),f_fft(stopF)],[minval,maxval],'k')
            listHBO = 1:size(d1)/2;
            listHBR = (size(d1)/2+1) : size(d1);
            if SNRfft
                figure;plot(squeeze(mean(abs(yall),2)));hold on
                plot(startF:stopF,squeeze(mean(abs(yall(startF:stopF,:,:)),2)),'x')
            end
             startFr =sum(f_fft<=freq(1));
            stopFr  =sum(f_fft<=freq(end));
            if SNRfft
           [peakfft,peakindividual] = max(nanmean(power(startFr:stopFr,:,:),3))
           
           %idpeak = (startFr + peakindividual)
           % peakfft = squeeze(nanmean(nanmean(abs(yall(startF:stopF,:,:)),3),1));
           % peakstd =squeeze(std(nanmean(abs(yall(startF+4:end,:,:)),3),1));
           % snrpeak = peakfft./peakstd;  
            
            % peakfft = squeeze(nanmean(nanmean(power(startF:stopF,:,:),3),1));
            peakstd =squeeze(std(nanmean(power(startF+3:end,:,:),3),1));
            snrpeak = peakfft./peakstd;  
            
            
           % figure;plot(nanmean(power(10:end, 2 ,:),3))
         idbad = find(snrpeak<2.5)
         idgood = find(snrpeak>2.5)
        figure;hold on
        cla
         for ipeak=1:numel(snrpeak)
             if snrpeak(ipeak)>2.5
                plot(f_fft(10:end),nanmean(log10(power(10:end,  ipeak ,:)),3),'displayname',[num2str(ipeak),'SNR',num2str(snrpeak(ipeak))] )

             end
         end
         title('SNR > 2.5')
         figure;hold on
         for ipeak=1:numel(snrpeak)
             if snrpeak(ipeak)<2.5
                 plot(f_fft(10:end),nanmean(log10(power(10:end,  ipeak ,:)),3),'displayname',[num2str(ipeak),'SNR',num2str(snrpeak(ipeak))] )

             end
         end
         title('SNR < 2.5')
            end
            if 0
                startF = repmat((interval(id))-nbpoint,size(yall,2),1);
                stopF = repmat((interval(id))+nbpoint,size(yall,2),1);;
           else 
               startF = (interval(ideach))-nbpoint; %fine tuning take peak on channel individual
                stopF = (interval(ideach))+nbpoint;
           end 
         figure(hfig)
            for i=1:numel(listHBO)
                if listHBO(i)
                    j = 1; 
                    in1 =  nanmean(yall(startF(i):stopF(i),listHBO(i),:),1);
                    while j<i  
                        if listHBO(j)
                            in2 = nanmean(yall(startF(j):stopF(j) ,listHBO(j),:),1);                      
                            COVC1C2  = nansum(in1.*conj(in2),3);
                            COVC1 =nansum(in1.*conj(in1),3);
                            COVC2 =nansum(in2.*conj(in2),3);
                            matcorr(i,j,f) = abs(COVC1C2).^2 ./ (COVC1.*COVC2);
                            matcorr(j,i,f) = abs(COVC1C2).^2 ./ (COVC1.*COVC2);
                            id = id+1;
                            if 0% listHBO(i) ==3
                            figure;
                            subplot(2,2,1)
                            plot(squeeze(nanmean(yall(startF:stopF ,listHBO(i),:),1)),'x');
                            vallim = max([real(squeeze(nanmean(yall(startF:stopF ,listHBO(i),:),1)));imag(squeeze(nanmean(yall(startF:stopF ,listHBO(i),:),1)))]);
                            xlim([- vallim,  vallim])
                            ylim([- vallim,  vallim])
                            title(['CH',num2str(listHBO(i))])
                            subplot(2,2,2)
                            plot(squeeze(nanmean(yall(startF:stopF ,listHBO(j),:),1)),'x');
                              xlim([- vallim,  vallim])
                            ylim([- vallim,  vallim])
                            
                              %title(['CH',num2str(listHBO(j))])
                              title(['Coherence matix',num2str(NIRS.Cf.dev.fs),'nm'])
                            subplot(2,2,3);hold on
                            plot(squeeze(in1.*conj(in2)),'x')
                            vallim = max([real(squeeze(in1.*conj(in2)));imag(squeeze(in1.*conj(in2)))])
                            xlim([- vallim,  vallim])
                            ylim([- vallim,  vallim])
                            subplot(2,2,4)
                            plot([0,real(COVC1C2)./ (COVC1.*COVC2)], [0, imag(COVC1C2)./ (COVC1.*COVC2)])
                            val = abs(COVC1C2).^2 ./ (COVC1.*COVC2)
                            title([ num2str(val)])
                            vallim = max(abs([0,real(COVC1C2)./ (COVC1.*COVC2),  imag(COVC1C2)./ (COVC1.*COVC2)]));
                            xlim([- vallim,  vallim])
                            ylim([- vallim,  vallim])
                            end
                        end
                        j = j + 1;
                    end 
                end
            end
%             
            for i=1:numel(listHBR)
                if listHBR(i)
                    j = 1;
                    in1 =  nanmean(yall(startF:stopF,listHBR(i),:),1);
                    while j<i %1:numel(listelectrode)
                        if listHBR(j)
                            in2 = nanmean(yall(startF:stopF ,listHBR(j),:),1);
                            COVC1C2  = nansum(in1.*conj(in2),3);
                            COVC1 =nansum(in1.*conj(in1),3);
                            COVC2 =nansum(in2.*conj(in2),3);
                            matcorrHbR(i,j,f) = abs(COVC1C2).^2 ./ (COVC1.*COVC2);
                            matcorrHbR(j,i,f) = abs(COVC1C2).^2 ./ (COVC1.*COVC2);
                            id = id+1;
                        end
                        j = j + 1;
                    end
                end
            end  
            
            nbch = size(matcorr,1);
              measlistok = sum([matcorr(:,:,f)>COHtr])>(nbch*pourcentagetr )& sum([matcorrHbR(:,:,f)>COHtr])>(nbch*pourcentagetr );
              measlistok = reshape(measlistok,1,numel(measlistok));
             subplot(5,4, [13,14,17,18]);hold on
             xlabel('CH id','fontsize',12)
             ylabel('CH id','fontsize',12)
             set(gca,'fontsize',12)
                         title(['Coherence matrix ', num2str(NIRS.Cf.dev.wl(1)),'nm'])

            imagesc( matcorr(:,:,f) ), caxis([0 0.2]); %>COHtr
            if ~isempty(find( measlistok==0))
                    %title(['Rejected : ',num2str(find( measlistok==0))]);
            end
            subplot(5,4, [15,16,19,20]);hold on
            imagesc( matcorrHbR(:,:,f)), caxis([0 0.2]); % >COHtr
            title(['Coherence matrix ', num2str(NIRS.Cf.dev.wl(2)),'nm'])

            xlabel('CH id','fontsize',12)
            set(gca,'fontsize',12)
               [filepath,name,ext] = fileparts(rDtp{f,1});
            xlsall{idxls,1} =  name;
           idxls = idxls + 1;
            xlsall{idxls,1} = ['PEAK cardiac='];
            xlsall{idxls,2} = num2str(tablepeak(f,1));
            chremove = find( measlistok==0);
              
            for i=1:numel(chremove)
                idxls = idxls + 1;
                ich =  chremove(i);
                switch NIRS.Cf.dev.n  
                case 'ISS Imagent' 
                    strDet = SDDet2strboxy_ISS(ML_new( ich,2));
                    strSrs = SDPairs2strboxy_ISS(ML_new( ich,1));
                case 'NIRx'                     
                    strDet = SDDet2strboxy(ML_new( ich,2));
                    strSrs = SDPairs2strboxy(ML_new( ich,1));                  
                otherwise
                    strDet = SDDet2strboxy(ML_new( ich,2));
                    strSrs = SDPairs2strboxy(ML_new( ich,1));  
                end
                xlsall{idxls,1} =  sprintf('%s\t%s',strDet,strSrs);  
                xlsall{idxls,2} = ich;
            end
            if  SNRfft
                idbadsnr = find(snrpeak<2.5);
                xlsallSNR{idxlsSNR,1} = name;
                idxlsSNR = idxlsSNR + 1;
                xlsallSNR{idxlsSNR,1} = ['PEAK cardiac='];
                xlsallSNR{idxlsSNR,2} = num2str(tablepeak(f,1));
                chremove=find(snrpeak<2.5)
              
            for i=1:numel(chremove)
                idxlsSNR = idxlsSNR + 1;
                ich =  chremove(i);
                switch NIRS.Cf.dev.n  
                case 'ISS Imagent' 
                    strDet = SDDet2strboxy_ISS(ML_new( ich,2));
                    strSrs = SDPairs2strboxy_ISS(ML_new( ich,1));
                case 'NIRx'                     
                    strDet = SDDet2strboxy(ML_new( ich,2));
                    strSrs = SDPairs2strboxy(ML_new( ich,1));                  
                otherwise
                
                    strDet = SDDet2strboxy(ML_new( ich,2));
                    strSrs = SDPairs2strboxy(ML_new( ich,1));  
                end
                xlsallSNR{idxlsSNR,1} =  sprintf('%s\t%s',strDet,strSrs);  
                xlsallSNR{idxlsSNR,2} = ich;
            end
            end
            
            NIRS.Cf.H.C.ok(:,f)=[measlistok,measlistok];
            [filepath,name,ext] = fileparts(job.NIRSmat{1});
            saveas(hfig,fullfile(filepath,[nametmp 'CardiacCHCOH',num2str(f),'.fig']),'fig'); 
            saveas(hfig,fullfile(filepath,[nametmp 'CardiacCHCOH',num2str(f),'.jpg']),'jpg');
            close(hfig)
    end
    if ismac
        writetxtfile_asxlsfile(fullfile(filepath,['CardiacCHCOH.xls']),xlsall);
    else
        xlswrite(fullfile(filepath,['CardiacCHCOH.xls']),xlsall);
        if  SNRfft
            xlswrite(fullfile(filepath,['CardiacCHSNR.xls']),xlsallSNR);
        end
    end
    disp(['Cardiac report : ', fullfile(filepath,['CardiacCHCOH',num2str(f),'.jpg']), ' and ', fullfile(filepath,['CardiacCHCOH.xls']), ' are created'])
    save(job.NIRSmat{1},'NIRS');
    out.NIRSmat = job.NIRSmat;
end
end
