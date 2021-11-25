function out = nirs_run_average(job)
% Average over many epochs, following pre-time and post-time set by the
% user.
% avtype values : 2 -> average over separate files %to check
%                 1 -> average over multiple files %valid now !
%                 0 -> average over mutliple subjects %to check

prefix = 'AVE'; %for "averaged"
DelPreviousData = job.DelPreviousData;
% 
% if isfield(job.NewDirCopyNIRSTRUE,'CreateNIRSCopy')
%     NewNIRSdir = job.NewDirCopyNIRSTRUE.CreateNIRSCopy.NewNIRSdir;
%     disp(['Create directory for condition ',NewNIRSdir])
%     NewDirCopyNIRS = 1;
% else
%     NewDirCopyNIRS = 0;
% end
%  
trig = job.choiceave.trigger; %List of triggers for averaging
pretime = str2num(job.choiceave.pretime);
posttime = str2num(job.choiceave.posttime);
avtype = 1; %job.choiceave.avtype; remove the option 
%nirsformat = job.savenirs;
nirsformat = 0;
badintervalratio = job.choiceave.badintervalratio;
helpmemoryprob = 0 ;%job.choiceave.helpmemoryprob;


%For multiple subjects average, look if montage is the same for each NIRS.mat
% if avtype == 0 %NON VÉRIFIER JT
%     for filenb=1:size(job.NIRSmat,1) %For every specified NIRS.mat file
%         %Load NIRS.mat information
%         NIRS = [];
%         load(job.NIRSmat{filenb,1});
%         if filenb > 1
%             if size(previous_montage,2) ~= size(NIRS.Cf.H.C.id,2)
%                 errordlg('All NIRS.mat do not contain the same montage. Multiple subjects average requires that all subjects have the same montage. Average over multiple files will be used instead.','Averaging Error');
%                 avtype = 1;
%                 break
%             elseif previous_montage ~= NIRS.Cf.H.C.id
%                 errordlg('All NIRS.mat do not contain the same montage. Multiple subjects average requires that all subjects have the same montage. Average over multiple files will be used instead.','Averaging Error');
%                 avtype = 1;
%                 break
%             end
%         end
%         previous_montage = NIRS.Cf.H.C.id;
%     end
% end
% A = [];



for filenb = 1:size(job.NIRSmat,1) %For every specified NIRS.mat file
    %Load NIRS.mat information
    %     try
    NIRS = [];
    unused_blocks = [];
    
    load(job.NIRSmat{filenb,1});
    [dir2,tmp,tmp] = fileparts(job.NIRSmat{filenb,1});
    report_temp = ones(NIRS.Cf.H.C.N);
    if avtype == 0 %For multiple subjects average, look if montage is the same for each NIRS.mat
        if filenb > 1
            if previous_montage ~= NIRS.Cf.H.S.r.o.mm.p
                errordlg('All NIRS.mat do not contain the same montage. Multiple subjects average requires that all subjects have the same montage. Average over multiple files will be used instead.','Averaging Error');
                avtype = 1;
            end
        end
        previous_montage = NIRS.Cf.H.S.r.o.mm.p;
    end
    %use last step of preprocessing
    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    [pathstr, ~, ~] = fileparts(rDtp{1});
    NC = NIRS.Cf.H.C.N;
    fs = NIRS.Cf.dev.fs;
    %erase matrices from previous NIRS.mat
    av_tot = [];
    trialnb = 1;
    for f = 1:size(rDtp,1) %For every path of NIRS.mat file
        %obtain indices of triggers
        ind_trig = [];
        aux_trig = NIRS.Dt.fir.aux5{f};
        for b = 1:size(trig,2)
            ind_trig = [ind_trig; aux_trig(aux_trig == trig(b),2)];
        end
        if isempty(ind_trig)
%             errordlg(['The selected trigger(s) are not avalaible for :',rDtp{f,1}],'Averaging Error');
%             break
%             answer = inputdlg('Enter valid trigger','Enter valid trigger');
%             answer = str2num(answer{1});
%             for b = 1:size(trig,2)
%                 ind_trig = [ind_trig aux_trig(aux_trig==answer,2)];
%             end
     %       disp(['No trig found in file : ',rDtp{f,1}])
            report_temp (:,f) = 0; 
            
        else
        %erase matrices from previous filepaths
        av = [];
        stdav = [];
        av_ind = [];
        ind_trig = ind_trig(:);        
        istart = ind_trig - abs((round(pretime*fs)*ones(size(ind_trig))));
        if istart<=0
            istart = 1;
        end            
        istop = ind_trig + abs((round(posttime*fs)*ones(size(ind_trig))));
        
        av_ind = zeros(2,numel(istart));
        av_ind(1,:) = istart;
        av_ind(2,:) = istop;
        %             try
        if ~helpmemoryprob %Help memory problem == no
            try
            if job.choiceave.avg_datatype == 2 % Open DC                
                namefile = rDtp{f,1};
                d = fopen_NIR(namefile,NC); %Load whole block
            elseif job.choiceave.avg_datatype == 1 % Open AC
                [pathstr, name, ext] = fileparts(rDtp{f,1});
                namefile = fullfile(pathstr,[name,'AC',ext]);
                d = fopen_NIR(namefile,NC); %Load whole block
            elseif job.choiceave.avg_datatype == 0 % Open PH
                [pathstr, name, ext] = fileparts(rDtp{f,1});
                namefile = fullfile(pathstr,[name,'PH',ext]);
                d = fopen_NIR(namefile,NC); %Load whole block
            end
            catch
                msgbox(['The file : ', namefile,' does''nt exist'])
                return
            end 
            
            %Enlever les indices des trigs qui sont en dehors des donnees
            maxsample = size(d,2);            
            listtoremove = find(maxsample<av_ind(2,:)|av_ind(1,:)<1);
            if ~isempty(listtoremove)
                disp('Trig timing out of the data')
                av_ind(:,listtoremove)=[];
                ind_trig(listtoremove)=[];
                if isempty(ind_trig)
                    disp(['No trig found in file : ',rDtp{f,1}])
                end
            end
            %av_ind is a 2 x N matrix : start and end of each average
            %epoch
            for atrig = 1:length(ind_trig)
                [pathstr, name, ext] = fileparts(NIRS.Dt.fir.pp(lst).p{f,1});
                InfoTrig(trialnb).file = name;
                InfoTrig(trialnb).filenb = f;
                InfoTrig(trialnb).nbpoint = size(d,2);
                InfoTrig(trialnb).startpoint = av_ind(1,atrig); %pre
                InfoTrig(trialnb).stoppoint = av_ind(2,atrig);  %post
                InfoTrig(trialnb).indtrig   = ind_trig(atrig);  %ind trig
                for Ich = 1:NC
               %     try
                        if NIRS.Cf.H.C.ok(Ich,f) == 1
                            %Default baseline correction
                            if isfield(job.choiceave.c_baseline_corr, 'm_defaultbaseline_corr') %Remove mean pretime to baseline
                                A(Ich,:,trialnb) = d(Ich,av_ind(1,atrig):av_ind(2, atrig))...
                                    - nanmean(d(Ich,av_ind(1, atrig):av_ind(1,atrig)+round(pretime*fs)));
                            elseif isfield(job.choiceave.c_baseline_corr,'b_manualbaseline_corr') %Manual timing baseline correction
                                nbs_pretrig = round(job.choiceave.c_baseline_corr.b_manualbaseline_corr.e_base_pretime*fs);
                                nbs_posttrig = round(job.choiceave.c_baseline_corr.b_manualbaseline_corr.e_base_posttime*fs); 
                                if job.choiceave.c_baseline_corr.b_manualbaseline_corr.m_mean == 0            %Mean Nan
                                    A(Ich,:,trialnb) = d(Ich,av_ind(1,atrig):av_ind(2, atrig))...
                                        - nanmean(d(Ich,(ind_trig(atrig)+nbs_pretrig):(ind_trig(atrig)+nbs_posttrig)));       
                                elseif job.choiceave.c_baseline_corr.b_manualbaseline_corr.m_mean == 1      %Median NaN
                                    A(Ich,:,trialnb) = d(Ich,av_ind(1,atrig):av_ind(2, atrig))...
                                        - nanmean(d(Ich,(ind_trig(atrig)+nbs_pretrig):(ind_trig(atrig)+nbs_posttrig)));   
                                end    
                                if isnan(nanmean(d(Ich,(ind_trig(atrig)+nbs_pretrig):(ind_trig(atrig)+nbs_posttrig))))
                                    1;
                                end
                                %attention si zero NaN                               
                            elseif isfield(job.choiceave.c_baseline_corr,'m_nobaseline_corr') %No baseline correction
                                 A(Ich,:,trialnb) = d(Ich,av_ind(1,atrig):av_ind(2, atrig));
                                
                            end
                            
                            if numel(find(isnan(A(Ich,:,trialnb))))/numel(A(Ich,:,trialnb)) >= badintervalratio %For epochs with too many artefacts
                                A(Ich,:,trialnb) = NaN;
                                disp(['File ', num2str(f), ', Channel ', num2str(Ich), ', Epoch ', num2str(atrig), ' was not used for averaging due to poor quality.']);
                                 
                            end
                        else %If the channel is bad (flagged in NIRS.Cf.H.C.ok)
                            A(Ich,:,trialnb) = ones(1,numel(av_ind(1,atrig):av_ind(2,atrig)))*NaN;
                            report_temp(Ich,f) = 0;
                           
                        end
                end
                trialnb = trialnb+1;
            end             
        end   
       %1;
       
        if avtype == 2 %save data file for separate files average %NON VÉRIfIER JT
            av = nanmean(A,3);
            stdav = nanstd(A,0,3);
            nbevents = sum(~isnan(A),3);
            [dir1,fil1,ext1] = fileparts(rDtp{f,1});
            %.nirs format
            if nirsformat
                location = [dir1, filesep, 'av_ep_subject',num2str(f),'file', num2str(f), '_trig', num2str(trig), '.nirs'];
                fwrite_NIR_homer(location,av,NIRS);
            end
            %.nir format
%             if NewDirCopyNIRS
%                 dir2 = [dir1 filesep NewNIRSdir];
                if ~exist(dir2,'dir'), mkdir(dir2); end;
                outfile = fullfile(dir2,[prefix fil1 ext1]);
                outfile2 = fullfile(dir2,[prefix fil1 '_std' ext1]);
                outfile3 = fullfile(dir2,[prefix fil1 '_events' ext1]);
                outfilevmrk = fullfile(dir2,[prefix fil1 '.vmrk']);
                fileOut_nir = fullfile(dir2,[prefix fil1 '.nir']);                
%             else
%                 outfile = fullfile(dir1,[prefix fil1 ext1]);
%                 outfile2 = fullfile(dir1,[prefix fil1 '_std' ext1]);
%                 outfile3 = fullfile(dir1,[prefix fil1 '_events' ext1]);               
%                 outfilevmrk = fullfile(dir1,[prefix fil1 '.vmrk']);
%                 fileOut_nir = fullfile(dir1,[prefix fil1 '.nir']);
%             end
            if DelPreviousData
                delete(rDtp{f,1});
            end
            fwrite_NIR(outfile,av);
            fwrite_NIR(outfile2,stdav);
            fwrite_NIR(outfile3,nbevents);
            
            
            temp_markers{1}.Type = 'NewSegment';
            temp_markers{1}.Description = '';
            temp_markers{1}.Position = 1;
            temp_markers{1}.Size = 1;
            temp_markers{1}.ChNumber = 0;
            nirs_boxy_write_markers(outfilevmrk,... %Output file
                fileOut_nir,... %DataFile
                temp_markers);
            ind_dur_ch = [abs((round(pretime*fs)*ones(size(ind_trig)))),0,0];
            typemarker = trig;
            write_vmrk(outfilevmrk,'trigger',num2str(typemarker),ind_dur_ch);
            
            
            %add outfile name to NIRS
            if f == 1
                NIRS.Dt.fir.pp(lst+1).pre = 'Epoch averaging (single file)';
                NIRS.Dt.fir.pp(lst+1).job = job;
            end
            NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
            %show graph
            if 0
                t = linspace(-pretime,posttime,size(av,2));
                figure
                plot(t,av)
                title('Averaged Intensity')
            end
            A = [];
            trialnb = 1;
        end
        
        end     %fin filenb .nir
    end         %Fin boucle des fichier
end             %fin filenb nirs.mat

[dir1,fil1,ext1] = fileparts(rDtp{f,1});
%  if NewDirCopyNIRS
%         dir2 = [dir1 filesep NewNIRSdir];
%         if ~exist(dir2,'dir'), mkdir(dir2); end;
%  end


%%%%%%%%processing for multiple files average%%%%%%%%%%
if avtype == 1 %save data file for average over many files
    if isempty(A)
        disp('No event was found for this condition')
        out.NIRSmat = job.NIRSmat;
        return 
    end
        if isfield(job.choiceave.c_rejecttrial,'b_reject_trial')
            maxz = zeros(size(A,1),1);
            minz = zeros(size(A,1),1);
            nbtrialrejected = zeros(size(A,1),1);
            zthresh = abs(job.choiceave.c_rejecttrial.b_reject_trial.e_reject_outlier_threshold);
                for ich = 1:size(A,1)
                    x = A(ich,:,:);
                    z = (x(:)-nanmean(x(:)))./nanstd(x(:));               
                    minz(ich) = nanmin(z);
                    maxz(ich) = nanmax(z);
                    z = reshape(z,size(x,2),size(x,3));                 
                    tmin = nanmin(z);
                    tmax = nanmax(z);                
                    indNaN = find(tmin < -zthresh | tmax > zthresh);
                    if ~isempty(indNaN)
                        A(ich,:,indNaN) = NaN;
                        nbtrialrejected(ich) = numel(indNaN);
                    end              
                end
                if job.choiceave.c_rejecttrial.b_reject_trial.m_reject_outlier_printreport
                    hreport=figure;
                    subplot(2,1,1);
                    bar(maxz);
                    hold on
                    bar(minz);
                    xlabel('Channel number')
                    ylabel('Extreme value zscore')
                    subplot(2,1,2);
                    bar(nbtrialrejected);
                    xlabel('Channel number')
                    ylabel('Rejected trial')
                    title(['Total of ', num2str(sum(nbtrialrejected)), ' trial rejected'])
%                     if NewDirCopyNIRS
%                         dir2 = [dir1 filesep NewNIRSdir];
%                         filereport = fullfile(dir2,['Report_Reject_trial_zscore_',num2str(zthresh)]);
%                     else
%                         dir2 = [dir1 filesep];
                        filereport = fullfile(dir2,['Report_Reject_trial_zscore_',num2str(zthresh)]);
%                    end
                    saveas(hreport,[filereport,'.jpg'],'jpg');  
                    saveas(hreport,[filereport,'.fig'],'fig');   
                    disp(['Save event report: ',filereport])
                end
                close(hreport)
        end
       
        if size(A,3)>1
            av = nanmean(A,3);   
            stdav = nanstd(A,0,3);
            nbevents = sum(~isnan(A),3);
            mediannbevents = median(nbevents');
        else
            av = A(:,:,1);   
            stdav = zeros(size(squeeze(A(:,:,1))));
            nbevents = sum(~isnan(A),3);
            mediannbevents = median(nbevents');
        end
        disp(['We find ',num2str(size(A,3)), ' events in the data'])
 
        modeTvalue = job.choiceave.m_Tvalueoption; %0 against zero, 1 unpaired vs baseline Worstcase , 2 Avg baseline 

        NIRS.Cf.H.C.okavg = (mediannbevents >= ((trialnb-1)*job.choiceave.badchannelratio))';
    %   NIRS.Dt.fir.pp(lst+1).chok replace by NIRS.Cf.H.C.okavg
        zeronbevent = find(nbevents == 0);
        tval = av./(stdav./sqrt(nbevents));
        tval(zeronbevent) = 0;  

        
               
        %Mettre les NAN à zero
        if modeTvalue >= 1
            tval(:)=0;%initialisation tval
            for ich = 1:size(tval,1)      
                for istop = round(abs((round(pretime*fs)*ones(size(ind_trig))))/2);
                    ind_baseline_start = 1;
                    ind_baseline_stop = istop;
                    Baseline = A(ich,ind_baseline_start:ind_baseline_stop,:);        

                    if modeTvalue == 1 %worstcasetvalue %On garde le temps avec la plus fort variance dans le baseline
                        mean2 = nanmean(squeeze(Baseline)'); %Mean
                        [v2,id] = max(nanvar(squeeze(Baseline)'));
                        m2 = mean2(id);
                        nt2 = sum(~isnan(Baseline(:,id)));
                    elseif modeTvalue == 2
                        m2 = nanmean(nanmean(Baseline,2)); %Mean
                        dismeantime = nanmean(Baseline,2);
                        v2 = max(nanvar(dismeantime));            
                        nt2 = sum(~isnan(dismeantime));
                    end        
                    nt1 = nbevents(ich,:);
                    m1 = nanmean(A(ich,:,:),3);
                    v1 = nanstd(A(ich,:,:),0,3);     
                    %Equal variance
                    gl = (nt1+nt2)-2; %degrees of freedom
                    s = ((nt1-1).*v1 + (nt2-1).*v2)./((nt1+nt2)-2); %combined variance
                    denom = sqrt((s./nt1+s./nt2));
                    dm = m1-m2;
                    tvalue = dm./denom; %t value unpaired
                    tval(ich,:) = tvalue;
                end
            end     
        end
    %test   
%     h=plot(tval(chlist,:))
%     set(h,'displayname', num2str(modeTvalue));

    %Repport nb of trials
    hnbtrial = figure;
    set(gca,'fontsize',14)
    hold on
    for ichannel = 1:numel(mediannbevents)/2
        if NIRS.Cf.H.C.okavg(ichannel)
            plot(ichannel,(mediannbevents(ichannel)/(trialnb-1))*100,'xk','markersize',10,'linewidth',3)
        else
            plot(ichannel,(mediannbevents(ichannel)/(trialnb-1))*100,'xg','markersize',10,'linewidth',3)
        end
    end
    plot([1,numel(mediannbevents)/2],[job.choiceave.badchannelratio*100,job.choiceave.badchannelratio*100])
    xlabel('Channel')
    ylabel(['% Events = ',num2str(size(A,3))])
    ylim([-10,110])
%     if NewDirCopyNIRS
%         dir2 = [dir1 filesep NewNIRSdir];
%         filereport = fullfile(dir2,['Report_Nb_trial_median_per_channel_reject_threshold_',num2str(job.choiceave.badchannelratio)]);
%     else
%         dir2 = [dir1 filesep];
        filereport = fullfile(dir2,['Report_Nb_trial_median_per_channel_reject_threshold_',num2str(job.choiceave.badchannelratio)]);
%    end 
    title(['Channel with nb Trial >= ', num2str(job.choiceave.badchannelratio*100), '% are keeped (green one are rejected)' ])
    saveas(hnbtrial,[filereport,'.jpg'],'jpg');  
    saveas(hnbtrial,[filereport,'.fig'],'fig'); 
    disp(['Save event report: ',filereport])
    close(hnbtrial)
    
%     figure;plot(tval')
%     figure;plot(av')
%     figure;plot(stdav')
%     figure;plot(nbevents)
    %.nirs format
    if nirsformat
        location = [dir1, filesep, 'av_ep_subject',num2str(f),'file1-', num2str(f), '_trig', num2str(trig), '.nirs'];
        fwrite_NIR_homer(location,av,NIRS);
    end
    %.nir format
%     if NewDirCopyNIRS
%         dir2 = [dir1 filesep NewNIRSdir];
    if ~exist(dir2,'dir'), mkdir(dir2); end;
        outfile = fullfile(dir2,[prefix fil1 ext1]);
        outfile2 = fullfile(dir2,[prefix fil1 '_std' ext1]);
        outfile3 = fullfile(dir2,[prefix fil1 '_events' ext1]);
        outfile4 = fullfile(dir2,[prefix fil1 '_tval' ext1]);
        outfilevmrk = fullfile(dir2,[prefix fil1 '.vmrk']);
        outfilevhdr = fullfile(dir2,[prefix fil1 '.vhdr']);

        fileOut_nir = fullfile(dir2,[prefix fil1 '.nir']);
%     else
%         outfile = fullfile(dir1,[prefix fil1 ext1]);
%         outfile2 = fullfile(dir1,[prefix fil1 '_std' ext1]);
%         outfile3 = fullfile(dir1,[prefix fil1 '_events' ext1]);
%         outfile4 = fullfile(dir1,[prefix fil1 '_tval' ext1]);
%         outfilevmrk = fullfile(dir1,[prefix fil1 '.vmrk']);
%         outfilevhdr = fullfile(dir2,[prefix fil1 '.vhdr']);
%         fileOut_nir = fullfile(dir1,[prefix fil1 '.nir']);
%     end
    if DelPreviousData
        delete(rDtp{f,1});
    end
    save([outfile(1:end-3),'mat'],'A','InfoTrig','-mat');
    fwrite_NIR(outfile,av);
    fwrite_NIR(outfile2,stdav);
    fwrite_NIR(outfile3,nbevents);
    fwrite_NIR(outfile4,tval);
    disp(['Save average data: ',outfile])


    temp_markers{1}.Type = 'NewSegment';
    temp_markers{1}.Description = '';
    temp_markers{1}.Position = 1;
    temp_markers{1}.Size = 1;
    temp_markers{1}.ChNumber = 0;
    nirs_boxy_write_markers(outfilevmrk,... %Output file
        fileOut_nir,... %DataFile
        temp_markers);
 
        %ind_dur_ch = [abs((round(pretime*fs)*ones(size(ind_trig)))),0,0];
        ind_dur_ch = [abs((round(pretime*fs))),0,0];
    
    typemarker = trig;
    write_vmrk(outfilevmrk,'trigger',num2str(typemarker),ind_dur_ch);
    
    
    %add outfile name to NIRS
    NIRS.Dt.fir.pp(lst+1).pre = 'Epoch averaging (multiple files)';
    NIRS.Dt.fir.pp(lst+1).job = job;
    for i = 1:1 %size(rDtp,1)
        NIRS.Dt.fir.pp(lst+1).p{i,1} = outfile;
    end
    %show graph
    if 0
        t = linspace(-pretime,posttime,size(av_tot,2));
        figure
        plot(t,av_tot)
        title('Multiple files Averaged Intensity')
    end
    
    
    %Report of blocks used for averaging
    report = [];
  
    for jtrig = 1:f
        report{1,(jtrig+1)} = sprintf('%s%d', 'Block ',jtrig);
    end
    
    for jchannel = 1:(size(A,1)/2)
        report{(jchannel+1),1} = sprintf('%s%d', 'Channel ',jchannel);
        for jtrig = 1:f
            if report_temp(jchannel,jtrig) == 0
                report{(jchannel+1),(jtrig+1)} = sprintf('%d',0); 
            else
                report{(jchannel+1),(jtrig+1)} = sprintf('%d',1);
            
            end
        end
    end
    
    if filenb > 1
        for jfile = 1:filenb
            [dir2,~,~] = fileparts(job.NIRSmat{jfile,1});
            path = fullfile(dir2,['Report_Blocks_Used_for_Averaging_File',num2str(jfile)]);
            if ismac
                writetxtfile_asxlswrite(path, report)
            else
                xlswrite(path, report);
            end
        end
    else
        path = fullfile(dir2,'Report_Blocks_Used_for_Averaging');
        if ismac
            writetxtfile_asxlswrite(path, report)
        else
            xlswrite(path, report);
        end
        xlswrite(path, report);
    end
    A = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save NIRS.mat
% if NewDirCopyNIRS
if isfield(NIRS.Dt,'Video')
   NIRS.Dt = rmfield(NIRS.Dt,'Video');
end
if isfield(NIRS.Dt,'AUX')
   NIRS.Dt = rmfield(NIRS.Dt,'AUX');
end
if isfield(NIRS.Dt,'EEG')
   NIRS.Dt = rmfield(NIRS.Dt,'EEG');
end
 
% else
%     save(job.NIRSmat{filenb,1},'NIRS');
% end
if 1 %newtrig
    aux5 = [];
    for i = 1:numel(trig)
        aux5 = [aux5; trig(i),round(pretime*NIRS.Cf.dev.fs)];
    end
    NIRS.Dt.fir.aux5{1} = aux5;
end
   save(fullfile(dir2,'NIRS.mat'),'NIRS');
    job.NIRSmat{1} =fullfile(dir2,'NIRS.mat');


%%%%%%%%processing for multiple subjects average%%%%%%%%%%
if avtype==0 %save data file for multiple subjects average %NON vérifier
    try
        av = nanmean(A,3);
        stdav = nanstd(A,0,3);
        nbevents = sum(~isnan(A),3);
        [dir1,fil1,ext1] = fileparts(rDtp{f,1});
        %.nirs format
        if nirsformat
            location = [dir1, filesep, 'av_ep_subjects1-',num2str(filenb),'_trig', num2str(trig), '.nirs'];
            fwrite_NIR_homer(location,av_subj,NIRS);
        end
        %.nir format
%         if NewDirCopyNIRS
%             dir2 = [dir1 filesep NewNIRSdir];
            if ~exist(dir2,'dir'), mkdir(dir2); end;
            outfile = fullfile(dir2,[prefix fil1 ext1]);
            outfile2 = fullfile(dir2,[prefix fil1 '_std' ext1]);
            outfile3 = fullfile(dir2,[prefix fil1 '_events' ext1]);
            outfilevhdr = fullfile(dir2,[prefix fil1 '.vhdr']);

%         else
%             outfile = fullfile(dir1,[prefix fil1 ext1]);
%             outfile2 = fullfile(dir1,[prefix fil1 '_std' ext1]);
%             outfile3 = fullfile(dir1,[prefix fil1 '_events' ext1]);
%             outfilevhdr = fullfile(dir1,[prefix fil1 '.vhdr']);
%         end
        if DelPreviousData
            delete(rDtp{f,1});
        end
        fwrite_NIR(outfile,av);
        fwrite_NIR(outfile2,stdav);
        fwrite_NIR(outfile3,nbevents);
        try
            ChannelLabels = ConvertmlIDsrs2label(NIRS);
            SamplingInterval =floor(1000000/NIRS.Cf.dev.fs);
            nirs_boxy_write_vhdr(outfilevhdr,... %Output file
                        outfile,... %DataFile
                        outfilevmrk,... %MarkerFile,...
                        'nirs_run_average',... %Function that created the header
                        '',... %Channel Resolution
                        '',... %Channel Units
                        ChannelLabels,... %names given as a column of cells
                        SamplingInterval,...
                        size(d,2)); %SamplingInterval in microseconds
        catch; end
        %write new .vmrk file
        outfilevmrk = fullfile(dir1,[prefix fil1 '.vmrk']);
        fileOut_nir = fullfile(dir1,[prefix fil1 '.nir']);
        temp_markers{1}.Type = 'NewSegment';
        temp_markers{1}.Description = '';
        temp_markers{1}.Position = 1;
        temp_markers{1}.Size = 1;
        temp_markers{1}.ChNumber = 0;
        nirs_boxy_write_markers(outfilevmrk,... %Output file
            fileOut_nir,... %DataFile
            temp_markers);
        ind_dur_ch = [abs((round(pretime*fs)*ones(size(ind_trig)))),0,0];
        typemarker = trig;
        write_vmrk(outfilevmrk,'trigger',num2str(typemarker),ind_dur_ch);
        
        %show graph
        if 0
            t = linspace(-pretime,posttime,size(av,2));
            figure
            plot(t,av_subj)
            title('Multiple subjects Averaged Intensity')
        end
    catch
        disp('Failed to average over multiple subjects. Check if all subjects have the same montage and/or triggers.');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


out.NIRSmat = job.NIRSmat;
