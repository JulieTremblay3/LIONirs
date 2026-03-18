function out = nirs_run_writenirshomer(job)
% Write toolbox data from the last node as data as .nirs format
% Field present in the structure
% d : Array of raw intensity
% t : Array of data time point
% s : Array of stimulus onset N time points x N conditions
% SD structure
% SD.Lambda Wavelength used for data acquisition
% SD.nSrcs
% SD.nDets
% SD.SrcPos Array of source coordinates nSrcs x 3
% SD.DetPos Array of detector coordinates nDets x 3
% SD.MeasList %list of measurement channel: source idx/ detector id / unused/ wavelength
%NIRS_exportoption,NIRSname
load(job.NIRSmat{1})     % Matrice NIRS
%Big loop over all subjects

for Idx=1:size(job.NIRSmat,1)
    nbch = numel(NIRS.Cf.H.C.n);
    dall = [];
    Sessionid = numel(NIRS.Dt.fir.pp); %prendre la dernière session
    try %create a defauld batchHistory file with lionirs pipeline used
        if isfield(job,'m_NIRSBATCHhistory') %print session history a .m not all step could be save...
            if job.m_SNIRFBATCHhistory
                if isfield(job,'c_NIRSname')
                    if isfield(job.c_NIRSname,'c_NIRSnamespecific')
                        name =  job.c_NIRSname.b_NIRSnamespecific.e_NIRSnamespecific ;
                        outfilebatch= [job.f_writeNIRSdir{1},filesep,name,'BatchHistory.m']   ;
                    else
                        outfilebatch= fullfile(job.f_writeNIRSdir{1},'BatchHistory.m');
                    end
                else
                    outfilebatch= fullfile(job.f_writeNIRSdir{1},'BatchHistory.m');
                end

                fid = fopen(outfilebatch,'w');
                for ibatch=1:numel(NIRS.Dt.fir.pp)
                    jobinstruction = convertpretobatchname(NIRS.Dt.fir.pp(ibatch).pre);
                    jobtext = ['matlabbatch{',num2str(ibatch),'}.',jobinstruction,'.'];
                    fieldlist = fieldnames(NIRS.Dt.fir.pp(ibatch).job);
                    for idjob = 1:numel( fieldlist)
                        id = 1;
                        if ~isstruct(eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]))
                            if iscell(eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]))
                                logjob=    [jobtext, [fieldlist{idjob}],'={''',   eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]),'''}'];
                                fprintf(fid,'%s',string(logjob));
                                fprintf(fid,'\n');
                                %sprintf('%s',string(logjob))
                            elseif isstr(eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]))
                                logjob=    [jobtext, [ fieldlist{idjob}],'=''',   eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]),''''];
                                fprintf(fid,'%s',string(logjob));
                                fprintf(fid,'\n');
                                %sprintf('%sn',string(logjob))
                            elseif isnumeric(eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]))
                                tmp = eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]);
                                if numel(tmp) == 1
                                    logjob=    [jobtext, [fieldlist{idjob}],'=',   num2str(tmp)];
                                    fprintf(fid,'%s',string(logjob));
                                    fprintf(fid,'\n');
                                elseif size(tmp,1)<size(tmp,2)
                                    logjob=    [jobtext, [fieldlist{idjob}],'=[',   num2str(tmp),']'];
                                    fprintf(fid,'%s',string(logjob));
                                    fprintf(fid,'\n');
                                elseif size(tmp,1)>size(tmp,2)
                                    logjob=    [jobtext, [fieldlist{idjob}],'= [',   num2str(tmp'),']'];
                                    fprintf(fid,'%s',string(logjob));
                                    fprintf(fid,'\n');
                                end
                            end
                        elseif isstruct(eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]))
                            firstpart = ['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' ];
                            temp =[fieldlist{idjob}];
                            fieldlist2 = fieldnames( eval([firstpart,temp]));
                            temp = [temp '.' fieldlist2{1}];
                            while isstruct(eval([firstpart,temp])) & isempty(isstruct(eval([firstpart,temp])));
                                fieldlist2 = fieldnames( eval([firstpart,temp]));
                                temp = [temp '.' fieldlist2{1}];
                            end


                            if iscell(eval([firstpart,temp]))
                                logjob=    [jobtext, [temp],'={''',   eval([firstpart,temp]),'''}'];
                                fprintf(fid,'%s',string(logjob));
                                fprintf(fid,'\n');
                                %sprintf('%s',string(logjob))
                            elseif isstr(eval([firstpart,temp]))
                                logjob=    [jobtext, [temp],'=''',   eval([firstpart,temp]),''''];
                                fprintf(fid,'%s',string(logjob));
                                fprintf(fid,'\n');
                                %sprintf('%sn',string(logjob))
                            elseif isnumeric(eval([firstpart,temp]))
                                tmp = eval(eval([firstpart,temp]));
                                if numel(tmp) == 1
                                    logjob=    [jobtext, [temp],'=',   num2str(tmp)];
                                    fprintf(fid,'%s',string(logjob));
                                    fprintf(fid,'\n');
                                elseif size(tmp,1)<size(tmp,2)
                                    logjob=    [jobtext, [temp],'=[',   num2str(tmp),']'];
                                    fprintf(fid,'%s',string(logjob));
                                    fprintf(fid,'\n');
                                elseif size(tmp,1)>size(tmp,2)
                                    logjob=    [jobtext, [temp],'= [',   num2str(tmp'),']'];
                                    fprintf(fid,'%s',string(logjob));
                                    fprintf(fid,'\n');
                                end
                            end
                        end
                    end
                end
                fclose(fid);
                disp(['Create log LIONirs as a Batch.m file : ', outfilebatch ])
            end
        end
    catch

    end


    for f = 1:numel(NIRS.Dt.fir.pp(1,Sessionid).p)
        [pathstr, name, ext] = fileparts(NIRS.Dt.fir.pp(1,Sessionid).p{f});
        %nbf = sprintf('%02.0f',i)%nom du fichier data .nir
        d = fopen_NIR([pathstr,filesep,name,ext],nbch);
        if ~isdir(job.f_writeNIRSdir{1})
            mkdir(job.f_writeNIRSdir{1})
        end
           if isfield(job,'c_NIRSname')
            if isfield(job.c_NIRSname,'b_NIRSnamespecific')
                  name=  job.c_NIRSname.b_NIRSnamespecific.e_NIRSnamespecific; 
                  outfile= [job.f_writeNIRSdir{1},filesep,name,sprintf('b%02.0f',f),'.nirs'];        
            else
                 if ~isdir([job.f_writeNIRSdir{1},filesep,name])
                     mkdir([job.f_writeNIRSdir{1},filesep,name]);
                 end 
                outfile= [job.f_writeNIRSdir{1},filesep,name,filesep,name,'.nirs'];
            end
        else
            outfile= [job.f_writeNIRSdir{1},filesep,name,filesep,name,'.nirs'];
             if ~isdir([job.f_writeNIRSdir{1},filesep,name])
                mkdir([job.f_writeNIRSdir{1},filesep,name]);
            end 
        end

       
        % check for channel list 
        if isfield(job.c_NIRSchannellist,'b_NIRSchannellist')
                 ML_actuel= [NIRS.Cf.H.C.id(2:3,:)',...
                    ones(size(NIRS.Cf.H.C.id,2),1),...
                    [ones(size(NIRS.Cf.H.C.id,2)/2,1);ones(size(NIRS.Cf.H.C.id,2)./2,1).*2]];
                fid = fopen(job.c_NIRSchannellist.b_NIRSchannellist.f_NIRSchannellist{1});
                chlist = textscan(fid, '%s%s');
                fclose(fid);
                DetL= chlist{1};
                SrsL= chlist{2};
                name =  DetL{1};
                if numel(name)>1
                    if strcmp(name(1:2),'D0')
                        Devicename = 'NIRx';
                    else
                        Devicename  = 'ISS Imagent';
                    end
                else
                    Devicename  = 'NIRx';
                end


         list = zeros(numel(chlist{1}),1);
         ML_new_hbo = [];
         ML_new_hbr = [];
        for i=1:numel(chlist{1})
            SDdetL = StrBoxy2SDDet_ISS(DetL{i});
            SDsrsL = StrBoxy2SDPairs(SrsL{i});
            switch  Devicename
                case 'ISS Imagent'
                    SDdetL = StrBoxy2SDDet_ISS(DetL{i});
                    SDsrsL = StrBoxy2SDPairs(SrsL{i});
                case 'NIRx'
                    SDdetL = StrBoxy2SDDet(DetL{i});
                    tmp = SrsL{i};
                    SDsrsL =str2num(tmp(2:end));
                case 'nirs'
                    SDdetL = StrBoxy2SDDet(DetL{i});
                    tmp = SrsL{i};
                    SDsrsL =str2num(tmp(2:end));
                otherwise
                    SDdetL = StrBoxy2SDDet_ISS(DetL{i});
                    SDsrsL = StrBoxy2SDPairs(SrsL{i});
            end
               L1 = find(ML_actuel(:,1)==SDsrsL & ML_actuel(:,2)==SDdetL & ML_actuel(:,4)==1);
               L2 = find(ML_actuel(:,1)==SDsrsL & ML_actuel(:,2)==SDdetL & ML_actuel(:,4)==2);
               ML_new_hbo = [ML_new_hbo; SDsrsL,SDdetL,1,1 ];
               ML_new_hbr = [ML_new_hbr; SDsrsL,SDdetL,1,2 ];
            if isempty(L1)
                sprintf(['check ', DetL{i},' ', SrsL{i}]);
                listname{i,1} = [DetL{i} ' ' SrsL{i}];
                listHBOch(i,1)= 0;
                listHBRch(i,1)= 0;
                listname{i,1} = [DetL{i} ' ' SrsL{i}];
            else
                listHBOch(i,1)= L1; 
                listHBRch(i,1)= L2;
                listname{i,1} = [DetL{i} ' ' SrsL{i}];
            end
        end
            idabsent = find(listHBOch==0);
            listHBOtmp =listHBOch;
            listHBOtmp( idabsent) = 1;
            idokHBO = NIRS.Cf.H.C.ok(listHBOtmp,f);
            idokHBO(idabsent)=0;
            idabsent = find(listHBRch==0);
            listHBRtmp =listHBRch;
            listHBRtmp( idabsent) = 1;
            idokHBR = NIRS.Cf.H.C.ok(listHBRtmp,f);
            idokHBR(idabsent)=0;
            listHBO = listHBOch;%(idok);
            listHBR = listHBRch;%(idok);
            ML_new = [ML_new_hbo;ML_new_hbr];
        end
          
        %%%



        d = d'; %data time point x channels
        SD = [];
        SD.Lambda = NIRS.Cf.dev.wl;
        SD.SrcPos = NIRS.Cf.H.S.r.o.mm.p';
        SD.DetPos = NIRS.Cf.H.D.r.o.mm.p';
        SD.nSrcs = size(SD.SrcPos,1);
        SD.nDets = size(SD.DetPos,1);
        if numel(NIRS.Cf.H.C.id(2,:)) == size(d,2)
            ml = [NIRS.Cf.H.C.id(2,:)',NIRS.Cf.H.C.id(3,:)', ones(numel(NIRS.Cf.H.C.id(2,:)),1),NIRS.Cf.H.C.wl']; %srs,det,1,wav
            t = 1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:size(d,1)*1/NIRS.Cf.dev.fs;
        end

        aux5 =  NIRS.Dt.fir.aux5{1}; %trigger info to set in s structure
        %Example; [2, 9;
        %2, 1733;
        %1,1]
        %Trig 2 sample of time 2 et 1733,
        %Trig 1 sample
        nbmaxtrig =  max(aux5(:,1));
        s = zeros(numel(t), nbmaxtrig);
        for i = 1:size(aux5,1)
            s(aux5(i,2),aux5(i,1))=1;
        end
        t = t';
        aux = zeros(size(s,1),1);
        SD.MeasList = ml;
        %force channel list order 
        %not change position or optode,  change channel order and ml
        if  isfield(job.c_NIRSchannellist,'b_NIRSchannellist')
               % ML_new         
                d1ok = nan(size(d,1),numel(listname)*2);        
               d1ok(:,find(listHBO)) = d(:,listHBO(find(listHBO))); 
               idHBR = find(listHBR)+ size(ML_new,1)/2;
               d1ok(:,idHBR)= d(:,listHBR(find(listHBR))); 
               SD.MeasList = ML_new ;            
                d = d1ok;
        end




        save(outfile,'d','SD','t','s','aux','-mat');
        fprintf('Saving %s ...\n', outfile);
    end
end



out.NIRSmat = job.NIRSmat;