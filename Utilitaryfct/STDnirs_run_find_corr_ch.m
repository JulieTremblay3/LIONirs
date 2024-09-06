function out = nirs_run_find_corr_ch(job)
% Displays the correlation between the coupled values of the .nir file
%

%filename prefix 
prefix = 'corr'; %for "correlated"
DelPreviousData  = job.DelPreviousData;
try 
    NewNIRSdir = job.NewNIRSdir;
    NewDirCopyNIRS = 1;
catch
    NewDirCopyNIRS = 0;
end

threshold = job.corr_thr;

for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
    %Load NIRS.mat information
%     try
        NIRS = [];
        load(job.NIRSmat{filenb,1});
        
        ind_dur_ch = [];
        
        %use last step of preprocessing
        lst = length(NIRS.Dt.fir.pp);
        rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
        NC = NIRS.Cf.H.C.N;
        fs = NIRS.Cf.dev.fs;
        win = floor(2.5*fs); %floor(2.5*fs); %2 times win = 5s moving window
        
        for f=1:size(rDtp,1) %Loop over all files of a NIRS.mat
%             try
                d = fopen_NIR(rDtp{f,1},NC);
                samp_length = size(d,2);
                %Erase previous matrices
                rho = []; 
                ind_dur_ch = [];
                a = 1;
                
                %new .vmrk file path
                [dir1,fil1,ext1] = fileparts(rDtp{f});
                infilewmrk = fullfile(dir1,[fil1 '.vmrk']);
                outfilevmrk = fullfile(dir1,[prefix fil1 '.vmrk']);
                try
                    copyfile(infilewmrk,outfilevmrk);
                catch %In the case of a module of nirs10 has been used just before (no .vmrk file have been updated).
                    rDtp2 = NIRS.Dt.fir.pp(lst-1).p;
                    [dir2,fil2,ext2] = fileparts(rDtp2{f});
                    infilewmrk = fullfile(dir2,[fil2 '.vmrk']);
                    copyfile(infilewmrk,outfilevmrk);
                end
                
                hwaitbar = waitbar(0);  
                for Idx = 1:NC/2 %Loop over all channel couples
                    waitbar(Idx/NC*2,hwaitbar,['Looking for correlated channels',' (',num2str(f),'/',num2str(size(rDtp,1)),')']);
%                     if NIRS.Cf.H.C.ok(Idx,f) ~= 0;  %If the channel was not previously flagged.
                        d1 = d(Idx,:)'; 
                        d2 = d(Idx+NC/2,:)';
                            %alpha = std(d1)/std(d2);
                        for t = 1:win
                            rho(Idx,t) = std(d1(1:t*2),d2(1:t*2));
                        end
                        for t = 1+win:samp_length-win
                            rho(Idx,t) = std(d1(t-win:t+win),d2(t-win:t+win));
                        end
                        for t = size(d,2)-win+1:samp_length
                            rho(Idx,t) = std(d1(samp_length-(samp_length-t)*2:samp_length),d2(samp_length-(samp_length-t)*2:samp_length));
                        end
                        %     ave_rho(1,Idx) = mean(rho(Idx,:)); %Mean and
                        %     Global correlations 
                        %     ave_rho(2,Idx) = corr(d1,d2);
                        
                        %%%%%%%%%%%TEST Xu Cui's method%%%
                        % x0(Idx,:) = 0.5.*(d1-alpha.*d2);
                        % y0(Idx,:) = -x0(Idx,:)./alpha;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        %Detection of very high correlation between concentrations
                        threshold = 0.8*max(rho(Idx,:));
                        highrho = find(rho(Idx,:) >= threshold);
                        diff_highrho = diff(highrho);
                        dur = 0;
                        for i = 1:numel(diff_highrho)
                            if diff_highrho(i) == 1
                                dur = dur + 1;
                            else
                                if dur > 2 % Neglect too small intervals
                                    ind_dur_ch(a,1) = highrho(i)-dur;
                                    ind_dur_ch(a,2) = dur+1;
                                    ind_dur_ch(a,3) = Idx;
                                    a = a+1;
                                end
                                dur = 0;
                            end
                        end
%                     end
                end
                close(hwaitbar);
                
                %write high correlation as artefacts in .vmrk file
                write_vmrk(outfilevmrk,'bad_step','high_corr',ind_dur_ch);
                               
                %%%%%%%%%%%%%%%%%%%%%
                if NewDirCopyNIRS
                    dir2 = [dir1 filesep NewNIRSdir];
                    if ~exist(dir2,'dir'), mkdir(dir2); end; 
                    outfile = fullfile(dir2,[prefix fil1 ext1]);
                else
                    outfile = fullfile(dir1,[prefix fil1 ext1]);
                end
                if DelPreviousData
                    delete(rDtp{f,1});
                end
                fwrite_NIR(outfile,d);
                               
                %add outfile name to NIRS
                if f == 1
                    NIRS.Dt.fir.pp(lst+1).pre = 'Find correlation';
                    NIRS.Dt.fir.pp(lst+1).job = job;
                end
                NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
        end
        
        if NewDirCopyNIRS
            save(fullfile(dir2,'NIRS.mat'),'NIRS');           
        else
            save(job.NIRSmat{filenb,1},'NIRS');
        end
%      catch
%          disp(['Could not find correlation for subject' int2str(filenb)])
%     end
end
out.NIRSmat = job.NIRSmat;

                
                
%CODE DANS ARTEFACT_REJECTION (HOMER)
% for i = 1:numel(NC)   %Parcourt tous les canaux et calcule leur corrélation avec les autres
%     d1ok = d(i,:)';    %On sélectionne un canal valide
%     for j = i:numel(NC)
%         d2ok = d(:,j);
%          matcorr(ch,ch2) = corr(d1ok,d2ok); %Corrélation entre les deux canaux
%          matcorr(ch2,ch) = matcorr(ch,ch2); %Matrice symétrique
%     end
% end