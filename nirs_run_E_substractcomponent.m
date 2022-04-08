function out = nirs_run_E_substractcomponent(job)
% Automatic step detection following a threshold set by the user. Bad
% points indice, duration and related channel are written in the .vmrk
% file.

%filename prefix
prefix = 's'; %for "substrac"
DelPreviousData  = 0; %job.DelPreviousData;
figureon = 0;
% if isfield(job.NewDirCopyNIRSTRUE,'CreateNIRSCopy')
%     NewNIRSdir = job.NewDirCopyNIRSTRUE.CreateNIRSCopy.NewNIRSdir;
%     disp(['Create directory for condition ',NewNIRSdir])
%     NewDirCopyNIRS = 1;
% else
%     NewDirCopyNIRS = 0;
% end

for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
    %Load NIRS.mat information
    tic
    NIRS = [];
    clear PARCOMP
    clear PARCORR
    clear NIRS
    load(job.NIRSmat{filenb,1});
    [dirmat,~,~] = fileparts(job.NIRSmat{filenb,1});
    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N;
    fs = NIRS.Cf.dev.fs;
    fprintf('%s\n','File processed');
    %create a copy of all file to apply the correction
    for f = 1:size(rDtp,1)
        d = fopen_NIR(rDtp{f,1},NC);
        [dir1,fil1,ext1] = fileparts(rDtp{f});
   
            
            if ~exist(dirmat,'dir'), mkdir(dirmat); end
                outfile = fullfile(dirmat,[prefix fil1 ext1]);
                infilevmrk = fullfile(dir1,[fil1 '.vmrk']);
                outfilevmrk = fullfile(dirmat,[prefix fil1 '.vmrk']);
                copyfile(infilevmrk,outfilevmrk);    
                infilevhdr = fullfile(dir1,[fil1 '.vhdr']);
                outfilevhdr = fullfile(dirmat,[prefix fil1 '.vhdr']);
                copyfile(infilevhdr,outfilevhdr);  
            
        
        if f == 1
            NIRS.Dt.fir.pp(lst+1).pre = 'Manual Gui Substract Component ';
            NIRS.Dt.fir.pp(lst+1).job = job;
        end
        fwrite_NIR(outfile,d);     
        NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
  
    end
    %Label to find and substract
    [outfolder,fil1,ext1] = fileparts(outfile);
    [pathstr, name, ext] = fileparts(job.NIRSmat{filenb,1});
    try
    infileSelectedfactors = fullfile(pathstr,'SelectedFactors.mat');
    outfileSelectedfactors = fullfile(outfolder,'SelectedFactors.mat');
    copyfile(infileSelectedfactors,outfileSelectedfactors);
    catch
    end
    try
    infileCorrectionApply = fullfile(pathstr,'CorrectionApply.mat');
    outfileCorrectionApply = fullfile(outfolder,'CorrectionApply.mat');
    copyfile(infileCorrectionApply,outfileCorrectionApply);
    catch
    end
    
    labelid = job.i_substractcomponent_label;
    [pathstr, name, ext] = fileparts(job.NIRSmat{filenb,1});
    try
        load(outfileSelectedfactors);
    catch
        disp('Substract component, no component exist')
        out.NIRSmat = job.NIRSmat;
        return
    end
    
    ifile = zeros(size(PARCOMP));
    for icomp = 1:numel(PARCOMP) 
        ifile(icomp) = PARCOMP(icomp ).file;
    end
    [ifile, iorder]=sort(ifile);
    %SortallFactorbyfile
    moduleout = lst+1;
    currentfile = 0;
    %Apply the correction and save the file
    toremove = [];
    for icomporder=1:numel(PARCOMP) %Loop over all files of a NIRS.mat
        icomp= iorder(icomporder);
        file = ifile(icomporder);
        if currentfile~=file
            rDtp = NIRS.Dt.fir.pp(moduleout).p;
            try
                d = fopen_NIR(rDtp{file,1},NC)';
                currentfile = file;
            catch
                disp(['File : ',rDtp{f,1},' cannot be open'])
            end
        end
        type = job.i_substractcomponent_label;

        if ~isempty(strfind(PARCOMP(icomp).label, type))%si bon label
           fprintf('%s',['Substract: ', PARCOMP(icomp).label]);
            if strfind(PARCOMP(icomp).type,'PARAFAC')
                toremove = [toremove,icomp];
                samp_length = size(d,2);
                time = 1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:samp_length*1/NIRS.Cf.dev.fs;
                indt =PARCOMP(icomp).indt;%Time indice
                
                
                intensnorm = d(indt,:);
                X = 1:1:size(intensnorm,1);
                Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
                Mb2 =  intensnorm(1,:)'; %offset
                A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
                spar = intensnorm - A;
                A = PARCOMP(icomp).FacA;
                B = PARCOMP(icomp).FacB;
                C = PARCOMP(icomp).FacC;
                listgood = PARCOMP(icomp).listgood;
                ComponentToKeep = PARCOMP(icomp).ComponentToKeep;
                Ac = A(:,ComponentToKeep); Bc = B(:,ComponentToKeep); Cc = C(:,ComponentToKeep);
                [Xm]=nmodel(({Ac,Bc,Cc}));
                if figureon
                    figure;plot(d(indt, listgood));hold on
                    plot(Xm(:,:,1),'k')
                end
                data = cat(3,d(indt,1:end/2),d(indt,end/2+1:end));
                data(:,listgood,:) = data(:,listgood,:)-Xm;
                d(indt,:) = reshape(data,[numel(indt),size(d,2)]);
                 %new value substracted with xm
                try
                    load(outfileCorrectionApply)
                    newfile = 0;
                catch
                    %donot exist create the stucture
                    PARCORR.file= PARCOMP(icomp).file;
                    PARCORR.filestr =  PARCOMP(icomp).filestr;
                    PARCORR.module  = PARCOMP(icomp).module+1;  
                    PARCORR.modulestr = 'Substract Component'; 
                    PARCORR.listgood =  PARCOMP(icomp).listgood;  
                    PARCORR.indt = PARCOMP(icomp).indt; %indice de temps.
                    PARCORR.data = PARCOMP(icomp).data;
                    PARCORR.Xm = PARCOMP(icomp).Xm;
                    PARCORR.FacA = PARCOMP(icomp).FacA;
                    PARCORR.FacB = PARCOMP(icomp).FacB;
                    PARCORR.FacC = PARCOMP(icomp).FacC;
                    PARCORR.ComponentToKeep = PARCOMP(icomp).ComponentToKeep;
                    PARCORR.label= PARCOMP(icomp).label;
                    PARCORR.type = PARCOMP(icomp).type;
                    PARCORR.topo =  PARCOMP(icomp).topo;
                    newfile = 1;
                end
                if newfile == 0
                    id = numel(PARCORR);
                    PARCORR(id+1).file= PARCOMP(icomp).file;
                    PARCORR(id+1).filestr =  PARCOMP(icomp).filestr;
                    PARCORR(id+1).module  = PARCOMP(icomp).module +1; 
                    PARCORR(id+1).modulestr = 'Substract Component';  
                    PARCORR(id+1).listgood =  PARCOMP(icomp).listgood;  
                    PARCORR(id+1).indt = PARCOMP(icomp).indt; %indice de temps.
                    PARCORR(id+1).data = PARCOMP(icomp).data;
                    PARCORR(id+1).Xm = PARCOMP(icomp).Xm;
                    PARCORR(id+1).FacA = PARCOMP(icomp).FacA;
                    PARCORR(id+1).FacB = PARCOMP(icomp).FacB;
                    PARCORR(id+1).FacC = PARCOMP(icomp).FacC;
                    PARCORR(id+1).ComponentToKeep = PARCOMP(icomp).ComponentToKeep;
                    PARCORR(id+1).label= PARCOMP(icomp).label;
                    PARCORR(id+1).type = PARCOMP(icomp).type;
                    PARCORR(id+1).topo =  PARCOMP(icomp).topo;
                end 
                save(outfileCorrectionApply,'PARCORR');

            elseif strcmp(PARCOMP(icomp).type,'PCA')
                toremove = [toremove,icomp];
                samp_length = size(d,2);
                indt =PARCOMP(icomp).indt;%Time indice
                intensnorm = d(indt,:);    
                X = 1:1:size(intensnorm,1);
                Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
                Mb2 =  intensnorm(1,:)'; %offset
                A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
                spar = intensnorm - A;
                
                data = d(indt,:);
                u =PARCOMP(icomp).u ;
                s =PARCOMP(icomp).s ;
                v =PARCOMP(icomp).v ;
                listgood = [PARCOMP(icomp).listgood] ;
                lstSV = 1;
                temp = u(:,lstSV)*s(lstSV,lstSV)*v(:,lstSV)';
                data(:,listgood) = data(:,listgood)- temp;
                 d(indt,listgood) = data(:,listgood) ;
                 %new value substracted with xm
                try
                    load(outfileCorrectionApply)
                    newfile = 0;
                catch
                    %donot exist create the stucture
                    PARCORR.file= PARCOMP(icomp).file;
                    PARCORR.filestr =  PARCOMP(icomp).filestr;
                    PARCORR.module  = PARCOMP(icomp).module+1;  
                    PARCORR.modulestr = 'Substract Component';
                    PARCORR.listgood =  PARCOMP(icomp).listgood;  
                    PARCORR.indt = PARCOMP(icomp).indt; %indice de temps.
                    PARCORR.data = PARCOMP(icomp).data;
                    PARCORR.Xm = PARCOMP(icomp).Xm;
                    PARCORR.u = PARCOMP(icomp).u;
                    PARCORR.s = PARCOMP(icomp).s;
                    PARCORR.v = PARCOMP(icomp).v;
                    PARCORR.ComponentToKeep = PARCOMP(icomp).ComponentToKeep;
                    PARCORR.label= PARCOMP(icomp).label;
                    PARCORR.type = PARCOMP(icomp).type;
                    PARCORR.topo =  PARCOMP(icomp).topo;
                    newfile = 1;
                end
                if newfile == 0
                    id = numel(PARCORR);
                    PARCORR(id+1).file= PARCOMP(icomp).file;
                    PARCORR(id+1).filestr =  PARCOMP(icomp).filestr;
                    PARCORR(id+1).module  = PARCOMP(icomp).module +1; 
                    PARCORR(id+1).modulestr = 'Substract Component';  
                    PARCORR(id+1).listgood =  PARCOMP(icomp).listgood;   
                    PARCORR(id+1).indt = PARCOMP(icomp).indt; %indice de temps.
                    PARCORR(id+1).data = PARCOMP(icomp).data;
                    PARCORR(id+1).Xm = PARCOMP(icomp).Xm;
                    PARCORR(id+1).u = PARCOMP(icomp).u;
                    PARCORR(id+1).s = PARCOMP(icomp).s;
                    PARCORR(id+1).v = PARCOMP(icomp).v;
                    PARCORR(id+1).ComponentToKeep = PARCOMP(icomp).ComponentToKeep;
                    PARCORR(id+1).label= PARCOMP(icomp).label;
                    PARCORR(id+1).type = PARCOMP(icomp).type;
                    PARCORR(id+1).topo =  PARCOMP(icomp).topo;
                end
                save(outfileCorrectionApply,'PARCORR');

            elseif strcmp(PARCOMP(icomp).type,'GLM')
                   toremove = [toremove,icomp];                  
                try
                    load(outfileCorrectionApply);
                    newfile = 0;
                catch
                    %donot exist create the stucture                 
                    PARCORR.file= PARCOMP(icomp).file;
                    PARCORR.filestr= PARCOMP(icomp).filestr;
                    PARCORR.module  = PARCOMP(icomp).module; 
                    PARCORR.modulestr = 'Substract Component'; 
                    PARCORR.listgood  = PARCOMP(icomp).listgood;
                    PARCORR.beta =  PARCOMP(icomp).beta;
                    PARCORR.std =  PARCOMP(icomp).std;
                    PARCORR.AUX = PARCOMP(icomp).AUX;
                    PARCORR.indt = PARCOMP(icomp).indt; %indice de temps.
                    PARCORR.data = PARCOMP(icomp).data  ;
                    PARCORR.Xm = PARCOMP(icomp).Xm;
                    PARCORR.ComponentToKeep = PARCOMP(icomp).ComponentToKeep ;
                    PARCORR.idreg = PARCOMP(icomp).idreg;
                    PARCORR.label= PARCOMP(icomp).label;
                    PARCORR.type = 'GLM';
                    PARCORR.topo = PARCOMP(icomp).topo;
                    newfile = 1;
                end
                if newfile == 0
                    id = numel(PARCORR);
                    PARCORR(id+1).file= PARCOMP(icomp).file;
                    PARCORR(id+1).filestr= PARCOMP(icomp).filestr;
                    PARCORR(id+1).module  = PARCOMP(icomp).module;
                    PARCORR(id+1).modulestr = 'Substract Component'; 
                    PARCORR(id+1).listgood  = PARCOMP(icomp).listgood;
                    PARCORR(id+1).beta =  PARCOMP(icomp).beta;
                    PARCORR(id+1).std =  PARCOMP(icomp).std;
                    PARCORR(id+1).AUX = PARCOMP(icomp).AUX;
                    PARCORR(id+1).indt = PARCOMP(icomp).indt; %indice de temps.
                    PARCORR(id+1).data = PARCOMP(icomp).data  ;
                    PARCORR(id+1).Xm = PARCOMP(icomp).Xm;
                    PARCORR(id+1).ComponentToKeep = PARCOMP(icomp).ComponentToKeep ;
                    PARCORR(id+1).idreg = PARCOMP(icomp).idreg;
                    PARCORR(id+1).label= PARCOMP(icomp).label;
                    PARCORR(id+1).type = 'GLM';
                    PARCORR(id+1).topo = PARCOMP(icomp).topo;
                end
                 save(outfileCorrectionApply,'PARCORR');
                d(PARCOMP(icomp).indt(1):PARCOMP(icomp).indt(end) ,PARCOMP(icomp).listgood) = d(PARCOMP(icomp).indt(1):PARCOMP(icomp).indt(end) ,PARCOMP(icomp).listgood) -PARCOMP(icomp).Xm;
                %new value substracted with xm
            end
            
            if job.m_substract_and_offsetcorrection==2
                    fprintf('%s\r',[' done']); 
            end
            
            if job.m_substract_and_offsetcorrection==1   
                %listgood must be for the 2 wavelenght for the offset
                %correction 
                if strfind(PARCOMP(icomp).type,'PARAFAC')
                    listgood2wv = [ PARCOMP(icomp).listgood; PARCOMP(icomp).listgood+NC/2];
                else
                    listgood2wv = PARCOMP(icomp).listgood;
                end
                indstart = PARCOMP(icomp).indt(1);
                indstop = PARCOMP(icomp).indt(end); 
                intensnorm = d(PARCOMP(icomp).indt(1):PARCOMP(icomp).indt(end) ,listgood2wv);
                X = 1:1:size(intensnorm,1);
                Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
                Mb2 =  intensnorm(1,:)'; %offset
                A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
                spar = intensnorm - A;
                Xm = A-1;              
                if indstart == 1 & (indstop+1)< size(d,1)
                    d(PARCOMP(icomp).indt(1):PARCOMP(icomp).indt(end) ,listgood2wv) = spar + ones(numel(PARCOMP(icomp).indt),1)*d(indstart,listgood2wv) ;
                    d(indstop+1:end,listgood2wv) =  d(indstop+1:end,listgood2wv)-...
                        ones(size(d(indstop+1:end,1)))* (d(indstop+1,listgood2wv) - d(indstart,listgood2wv) ); %offset adjustement
                    dureetime = (PARCOMP(icomp).indt(end)-PARCOMP(icomp).indt(1))*1/fs;
                    fprintf(' %3.1f%s\r',dureetime,['s done and offset adjustement']); 
                elseif indstart == 1 & indstop== size(d,1)
                     d(PARCOMP(icomp).indt(1):PARCOMP(icomp).indt(end) ,listgood2wv) = spar + ones(numel(PARCOMP(icomp).indt),1)*d(indstart,listgood2wv) ;
                    d(indstop+1:end,listgood2wv) =  d(indstop+1:end,listgood2wv)-...
                        ones(size(d(indstop:end,1)))* (d(indstop,listgood2wv) - d(indstart,listgood2wv) ); %offset adjustement
                    dureetime = (PARCOMP(icomp).indt(end)-PARCOMP(icomp).indt(1))*1/fs;
                    fprintf(' %3.1f%s\r',dureetime,['s done and offset adjustement']); 
                elseif (indstop+1)>= size(d,1)
                    d(PARCOMP(icomp).indt(1):PARCOMP(icomp).indt(end) ,listgood2wv) = spar + ones(numel(PARCOMP(icomp).indt),1)*d(indstart-1,listgood2wv) ;
                    dureetime = (PARCOMP(icomp).indt(end)-PARCOMP(icomp).indt(1))*1/fs;
                    fprintf(' %3.1f%s\r',dureetime,['s  done ']);
                else
                    try
                    d(PARCOMP(icomp).indt(1):PARCOMP(icomp).indt(end) ,listgood2wv) = spar + ones(numel(PARCOMP(icomp).indt),1)*d(indstart-1,listgood2wv) ;
                    d(indstop+1:end,listgood2wv) =  d(indstop+1:end,listgood2wv)-...
                    ones(size(d(indstop+1:end,1)))* (d(indstop+1,listgood2wv) - d(indstart-1,listgood2wv) );
                    dureetime = (PARCOMP(icomp).indt(end)-PARCOMP(icomp).indt(1))*1/fs;
                    fprintf(' %3.1f%s\r ',dureetime,['s  done and offset adjustement']); 
                    catch
                        fprintf(' %3.1f%s\r ',dureetime,['s  Error']); 
                    end
                end
               
            
                  try
                    load(outfileCorrectionApply)
                    newfile = 0;
                catch
                    PARCORR.file=  PARCOMP(icomp).file;
                    PARCORR.filestr =  PARCOMP(icomp).filestr;
                    PARCORR.module  = PARCOMP(icomp).module;
                    PARCORR.modulestr  = 'Substract Component'; 
                    PARCORR.listgood = listgood2wv;
                    PARCORR.indt = indstart:indstop;%indice de temps.
                    PARCORR.Xm = Xm; 
                    if indstart == 1 & (indstop+1)< size(d,1)
                        PARCORR.Offset =  (d(indstop+1,:) - d(indstart,:) );
                        fprintf('%s\r',[' Offset Corrected ']); 
                    elseif  indstart == 1 & indstop== size(d,1)
                        PARCORR.Offset  =  (d(indstop,:) - d(indstart,:) );
                    elseif (indstop+1)>size(d,1)
                        PARCORR.Offset =  (d(end,:) - d(indstart-1,:) );
                    else
                        PARCORR.Offset =  (d(indstop+1,:) - d(indstart-1,:) );
                    end
                    newfile = 1;
                    PARCORR.label= ['Offset Adjustment' sprintf('%03.0f',size(PARCORR,2)),' ',PARCOMP(icomp).filestr];
                    PARCORR.type = 'Offset Adjustment';
                end
                if newfile == 0
                    id = numel(PARCORR);
                    PARCORR(1+id).file=  PARCOMP(icomp).file;           
                    PARCORR(1+id).filestr =  PARCOMP(icomp).filestr;
                    PARCORR(1+id).module = PARCOMP(icomp).module;             
                    PARCORR(1+id).modulestr = 'Substract Component'; 
                    PARCORR(1+id).listgood =  listgood2wv;
                    PARCORR(1+id).indt = indstart:indstop;
                    PARCORR(1+id).Xm = Xm; %indice de temps.
                    if indstart == 1 & (indstop+1)< size(d,1)
                        PARCORR(1+id).Offset =  (d(indstop+1,:) - d(indstart,:) );
                    elseif indstart == 1 & indstop== size(d,1)
                        PARCORR(1+id).Offset =  (d(indstop,:) - d(indstart,:) );
                    elseif (indstop+1)>size(d,1)
                        PARCORR(1+id).Offset =  (d(end,:) - d(indstart-1,:) );
                    else
                        PARCORR(1+id).Offset =  (d(indstop+1,:) - d(indstart-1,:) );
                    end
                    newfile = 1;
                    PARCORR(1+id).label= ['Offset Adjustment' sprintf('%03.0f',size(PARCORR,2)),' ',PARCOMP(icomp).filestr];
                    PARCORR(1+id).type = 'Offset Adjustment';
                end
                save(outfileCorrectionApply,'PARCORR');
            end
        end
        
       
        if icomporder < numel(PARCOMP) %save if the next one is different
            if ifile(icomporder)~= ifile(icomporder+1)
                outfile = rDtp{file,1};
                fwrite_NIR(rDtp{file},d');   
                   disp(['Save: ', outfile]);
            end
        end
        
        if icomporder == numel(PARCOMP)%save if last one
            outfile = rDtp{file,1};
            fwrite_NIR(rDtp{file,1},d');
             disp(['Save: ', outfile]);
        end

    end
    PARCOMP(toremove)=  [];
    %Apply module to ensure that the oparation is accessible form the
    %manual step
    for id=1:numel(PARCORR)
        PARCORR(id).module = length(NIRS.Dt.fir.pp);
        PARCORR(id).modulestr = 'Manual Gui Substract Component';
    end
     for id=1:numel(PARCOMP)
        PARCOMP(id).module =length(NIRS.Dt.fir.pp);
        PARCOMP(id).modulestr = 'Manual Gui Substract Component';
    end
   
    save(outfileSelectedfactors,'PARCOMP');
    save(outfileCorrectionApply,'PARCORR');
    
    
    %save the NIRS.mat new structure
    [dir1,~,~] = fileparts(rDtp{f});
    save(fullfile(dir1,'NIRS.mat'),'NIRS');
    job.NIRSmat{1} = fullfile(dir1,'NIRS.mat');
    %figure;plot(d);
end
out.NIRSmat = job.NIRSmat;

