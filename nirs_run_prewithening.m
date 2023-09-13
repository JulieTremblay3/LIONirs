%___________________________________________________________________
% Copyright (C) 2019 LION Lab, Centre de recherche CHU Sainte-Justine 
% www.lionlab.umontreal.ca
%___________________________________________________________________
function out = nirs_run_prewithening(job)
%filename prefix 
prefix = 'w'; %for "whithening"
DelPreviousData  = job.DelPreviousData;

% [lowcut,applylowcut] = str2num(job.lowcutfreq);
% [ highcut ,applyhighcut] = str2num(job.highcutfreq);
% 
% filt_ord = job.filterorder;
% paddingsym = job.paddingsymfilter; %symetrie padding on the signal to avoid edge on the filtering
% interpolate = job.interpolatebadfilter;  %interpolate bad interval
% for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
%     %Load NIRS.mat information
%         NIRS = [];
%         load(job.NIRSmat{filenb,1});
%    
%         %use last step of preprocessing
%         lst = length(NIRS.Dt.fir.pp);
%         rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
%         NC = NIRS.Cf.H.C.N;
%         fs = NIRS.Cf.dev.fs;
%         
%          %Verify job.NIRS.mat location 
%          %dir2 will be the new location 
%         [dir2,tmp,tmp] = fileparts(job.NIRSmat{filenb,1});
% 
%         
%         fprintf('%s\n',['File processed using filter low pass=', num2str(lowcut) ,' high pass=' , num2str(highcut),' order=', num2str(filt_ord),'  butterword zero phase filter (function butter.m and filtfilt.m)']);
%         if paddingsym
%             fprintf('%s\n','Using mirror padding')
%         end
%         for f=1:size(rDtp,1) %Loop over all files of a NIRS.mat 
% %             try
%         figon=0;
% 
%                 d = fopen_NIR(rDtp{f,1},NC);
%                 if figon
%                 figure;subplot(3,1,1);plot(d')
%                 end
%                  [dir1,fil1,~] = fileparts(rDtp{f});
%                 vmrk_path = fullfile(dir1,[fil1 '.vmrk']);
%                 [ind_dur_ch] = read_vmrk_find(vmrk_path,'bad_step');
%                 %d = NaN*ones(size(d));
%                 
%                 if paddingsym
%                     dtmp = [fliplr(d),d,fliplr(d)];
%                     tstart = size(d,2)+1;
%                     tstop = size(d,2)*2;
%                     %figure;plot(dtmp')
%                     d = dtmp;
%                 end 
%                 id = find(isnan(d));
%                 if ~isempty(id)
%                   d(id) = 0;
%                end
%               
%                 if ~isempty(ind_dur_ch)&&  interpolate == 1
%                     if paddingsym
%                         ind_dur_chpre  = ind_dur_ch;
%                         ind_dur_chpre(:,1)  = tstart-ind_dur_ch(:,1)-ind_dur_ch(:,2);
%                         ind_dur_chmid = ind_dur_ch;
%                         ind_dur_chmid(:,1) =ind_dur_ch(:,1)+tstart-1;
%                         ind_dur_chpost = ind_dur_ch;
%                         ind_dur_chpost(:,1) = (tstart-ind_dur_ch(:,1)-ind_dur_ch(:,2))+tstop;
%                         ind_dur_chtot = [ind_dur_chpre;ind_dur_chmid;ind_dur_chpost];
%                     end
%                     [dinterp] = interpolate_bad(d,ind_dur_chtot);
%                 else
%                     dinterp = d;
%                 end 
%                 if figon
%                     subplot(3,1,2);plot(d')
%                     subplot(3,1,3);plot(dinterp')
%                 end
% 
%                 for Idx = 1:NC %Loop over all channels
%                     if 1 %NIRS.Cf.H.C.ok(Idx,f) ~= 0;  %If the channel was not previously flagged. %DO ALL
%                         %LOWPASS
%                         if applylowcut
%                             Wn = lowcut*2/fs;
%                             [fb,fa]=butter(filt_ord,Wn);
%                             dfilt(Idx,:) = filtfilt(fb,fa,dinterp(Idx,:));                 
%                         else
%                             dfilt(Idx,:) = dinterp(Idx,:);
%                         end
%                         %HIGHPASS
%                         if applyhighcut
%                             Wn = highcut*2/fs;
%                             [fb,fa]=butter(filt_ord,Wn,'high');
%                             dfilt(Idx,:) = filtfilt(fb,fa,dfilt(Idx,:));
%                             dfilt(Idx,:) = dfilt(Idx,:);
%                         end
%                     end
%                 end
%                 %plot response
%                 if paddingsym
%                     dfilt = dfilt(:,tstart:tstop);
%                 end
%                 [dir1,fil1,ext1] = fileparts(rDtp{f});
%                 infilevmrk = fullfile(dir1,[fil1 '.vmrk']);
%                 try
%                 infilevhdr = fullfile(dir1,[fil1 '.vhdr']);    
%                 catch
%                 end
%                 outfile = fullfile(dir2,[prefix fil1 ext1]);
%                 outfilevmrk = fullfile(dir2,[prefix fil1 '.vmrk']);
%                 outfilevhdr = fullfile(dir2,[prefix fil1 '.vhdr']);
% 
%                
%                 fwrite_NIR(outfile,dfilt); 
%                 clear dfilt
%                 fprintf('%s\n',outfile);
%  
%                 %write new .vmrk file 
%                 try
%                 copyfile(infilevmrk,outfilevmrk);
%                 catch
%                 end
%                 try
%                 copyfile(infilevhdr,outfilevhdr);
%                     try
%                         info = read_vhdr_brainvision((fullfile(dir1,[fil1,'.vhdr'])));
%                         ChannelLabels = info.label;
%                     catch
%                         ChannelLabels = ConvertmlIDsrs2label(NIRS);
%                     end
%                 SamplingInterval =floor(1000000/NIRS.Cf.dev.fs);
%                 nirs_boxy_write_vhdr(outfilevhdr,... %Output file
%                             fileOut,... %DataFile
%                             outfilevmrk,... %MarkerFile,...
%                             'nirs_run_filter',... %Function that created the header
%                             '',... %Channel Resolution
%                             '',... %Channel Units
%                             ChannelLabels,... %names given as a column of cells
%                             SamplingInterval,...
%                             size(dfilt,2)); %SamplingInterval in microseconds
%     
%                 
%                 catch
%                 end
%                  if DelPreviousData
%                     try
%                     delete(rDtp{f,1});                    
%                     delete(infilevmrk);
%                     delete(infilevhdr);
%                     disp(['Delete previous .nir data file: ',rDtp{f,1}]);
%                     catch                        
%                     end
%                      
%                 end
%                 %add outfile name to NIRS
%                 if f == 1
%                     NIRS.Dt.fir.pp(lst+1).pre = 'PreWhitening';
%                     NIRS.Dt.fir.pp(lst+1).job = job;
%                 end 
%                 NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
%         end
%             save(fullfile(dir2,'NIRS.mat'),'NIRS');
%             job.NIRSmat{1} =fullfile(dir2,'NIRS.mat');
% 
% end
 out.NIRSmat = job.NIRSmat;

