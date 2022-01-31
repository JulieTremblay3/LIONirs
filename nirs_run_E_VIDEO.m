function out = nirs_run_E_VIDEO(job)
TOPOmat = []; %initialisation du fichier de sortie contenant le nom de fichier TOPO.mat créer lors de la video
out  = [];
fprintf('%s\n','File processed');
m_prjfile_mode = job.b_option.m_prjfile_mode; %prj paired == 0, prj unique == 1 
if isfield(job.typedata,'b_HMRdatad1manual')
    job.typedata.b_HMRdata = job.typedata.b_HMRdatad1manual;
end

label.projectiontype = job.b_option.m_projection_mode;

if isfield(job.typedata,'b_HMRdata')
    %try
    if m_prjfile_mode==0
        if numel(job.b_option.prjfile)~=numel(job.typedata.b_HMRdata.HMRdata)
            msgbox('Please enter same number of .prj than .HMR data file they should be corresponding if you used Option prj : paired mode');
            return
        end
    end      
        for ifile=1:numel(job.typedata.b_HMRdata.HMRdata)  
            clear PMI
            clear oMRI
            clear PrjStruct
%             if ifile > 1
%                 set(guiHOMER,'UserData',[]);
%             end
            load([job.typedata.b_HMRdata.HMRdata{ifile}],'-mat');
            if ~iscell(DOT)
                DOT1={};
                DOT1{1}=DOT;
                clear DOT;
                PMI{1}=DOT1{1};
                clear DOT1;
            else
                PMI = DOT;
                clear DOT;
            end
            if m_prjfile_mode==0    %Paired
                PMI{1}.prj_name = job.b_option.prjfile{ifile};
            elseif m_prjfile_mode==1 %Unique
                PMI{1}.prj_name = job.b_option.prjfile{1};
            end
            [p,file_prj_name,ext]=fileparts(PMI{1}.prj_name);

            cfifind = 0;
            filenmlist = [];
            for cfi=1:numel(PMI{1}.data)
                filenmlist = [filenmlist,' ',PMI{1}.data(cfi).filenm];
                if strcmp(PMI{1}.data(cfi).filenm,job.typedata.b_HMRdata.HMR_cffile)
                    cf=cfi;
                    PMI{1}.currentFile = cf;
                    cfifind =1;
                end
            end
            if cfifind==0
                msgbox(['File condition : ', job.typedata.b_HMRdata.HMR_cffile,' was not found in the data. ',...
                    'Process stop at subject : ',job.typedata.b_HMRdata.HMRdata{ifile}, '. The condition list for this subject is :', filenmlist ]);
                return
            end
            
            disp(['File condition :' PMI{1}.data(cf).filenm])
            %Overture du .prj IO_Helmet pour créer l'affichage
            global currentsub
            currentsub=1;
            hfig = figure;
            set(hfig,'visible','off')
            %PRJ
            setappdata(0,'gui_SPMvideo',hfig);
            guiHOMER = getappdata(0,'gui_SPMvideo');
            set(guiHOMER,'UserData',PMI); %PMI de base
            
            
            handles = [];
            option.conc = 1;
            option.scalemin = job.b_option.v_cmin;
            option.scalemax = job.b_option.v_cmax;
            option.scaletr  = job.b_option.v_ctr;
            option.skintype = job.b_option.v_skin;
            option.cover = job.b_option.v_showcover;
            step_time = job.b_option.v_step;
            option.video = 1; % la fonction qui appelle le casque est le mode video
            option.cover = job.b_option.v_showcover;
            hfighelmet = IO_HelmetMTG_Display(handles,option);
            set(hfighelmet,'visible','off');
            
            
            %possibilité d'ouvrir des fichier déjà créer en vcolor ?
            if isfield(job.typedata,'b_HMRdata')
                [pathavi, name, ext] = fileparts(job.typedata.b_HMRdata.HMRdata{ifile});
                file = [name];
            end
            
            start = job.b_option.v_start;
            postTime= job.b_option.v_stop-job.b_option.v_start;
            if  postTime < 0
                postTime = 1;
                msgbox('Please look at that the start time it should start before the stop time !')
            end
            
            [tok,rem]=strtok(job.b_option.v_outpath,':');
            if isempty(rem) %si le nom du repertoire est juste un nom à ajouter dans le fichier de donné
                label.pathout = [pathavi,'/',job.b_option.v_outpath,'/'];
            else %si le nom du repertoire est un répertoire complet
                try
                    if strcmp(job.b_option.v_outpath(end),filesep)||strcmp(job.b_option.v_outpath(end),filesep)
                        label.pathout = [job.b_option.v_outpath];
                    else
                        label.pathout = [job.b_option.v_outpath,filesep];
                    end
                catch
                    msgbox(['The name of the output directory : ',job.b_option.v_outpath,' is invalid']);
                end
            end
            if ~isdir(label.pathout) 
                   [status,message,messageid]= mkdir(label.pathout);
                   if status == 0
                       msgbox(['ERROR the path : ', label.pathout,'can''t be create']);
                   end
            end
            
            %ICI LES DONNÉES À UTILISER
            %On choisit de tout storer dans d1 pour mieux structure le code et les
            %nouveaux ajout.
            if job.typedata.b_HMRdata.HMR_datatype==1 || job.typedata.b_HMRdata.HMR_datatype==3 ||  job.typedata.b_HMRdata.HMR_datatype==5
                flag_avgtime = 0; %each time point
            elseif  job.typedata.b_HMRdata.HMR_datatype==2 || job.typedata.b_HMRdata.HMR_datatype==4 ||  job.typedata.b_HMRdata.HMR_datatype==6
                flag_avgtime = 1; %avg periode of time point
            end
            chok  = find(PMI{currentsub}.data(cf).MeasListAct(1:end/2));
            nbch = numel(PMI{currentsub}.data(cf).MeasListAct(1:end/2));
            %dConc
            if job.typedata.b_HMRdata.HMR_datatype==1 || job.typedata.b_HMRdata.HMR_datatype==2
                dHbO = PMI{currentsub}.data(cf).HRF.AvgC(:,chok ,1);
                dHbR = PMI{currentsub}.data(cf).HRF.AvgC(:,chok ,2);
                dHbT = PMI{currentsub}.data(cf).HRF.AvgC(:,chok ,3);
                ndHbO = PMI{currentsub}.data(cf).HRF.AvgC(:,chok ,1);
                ndHbR = PMI{currentsub}.data(cf).HRF.AvgC(:,chok ,2);
                ndHbT = PMI{currentsub}.data(cf).HRF.AvgC(:,chok ,3);
                typenorm = 'No';
                %zscore
            elseif job.typedata.b_HMRdata.HMR_datatype==3 || job.typedata.b_HMRdata.HMR_datatype==4
                dHbO = PMI{currentsub}.data(cf).HRF.AvgC(:,chok ,1);
                dHbR = PMI{currentsub}.data(cf).HRF.AvgC(:,chok ,2);
                dHbT = PMI{currentsub}.data(cf).HRF.AvgC(:,chok ,3);
                ndHbO = zscore(dHbO(:));
                ndHbO = reshape(ndHbO,size(dHbO));
                ndHbR = zscore(dHbR(:));
                ndHbR = reshape(ndHbR,size(dHbR));
                ndHbT = zscore(dHbT(:));
                ndHbT = reshape(ndHbT,size(dHbT));
                typenorm = 'Zscore';
                %min max
            elseif  job.typedata.b_HMRdata.HMR_datatype==5 || job.typedata.b_HMRdata.HMR_datatype==6
                typenorm = 'min max';
            end
            
            
            if flag_avgtime==0  %each data point
                itime = 1;
                clear d1;
                d1 = zeros(numel(start:step_time:(start+postTime)),nbch,3);
                for find_time = start:step_time:(start+postTime)
                    x = find_time;
                    echantillon_time = find(find_time >= PMI{currentsub}.data(cf).HRF.tHRF);
                    echantillon_time = echantillon_time(end);
                    labeltime{itime} = fixdecimal2string(PMI{currentsub}.data(cf).HRF.tHRF(echantillon_time),job.b_option.v_outtimeprecisionleft,job.b_option.v_outtimeprecision);
                    d1(itime,chok,1) = ndHbO(echantillon_time,:);
                    d1(itime,chok,2) = ndHbR(echantillon_time,:);
                    d1(itime,chok,3) = ndHbT(echantillon_time,:);
                    itime = itime + 1;
                end
            elseif flag_avgtime==1;  %average data time point
                itime = 1;
                clear d1;
                d1 = zeros(numel(start:step_time:(start+postTime)),nbch,3);
                for find_time = start:step_time:(start+postTime)
                    x = find_time;
                    echantillon_time = find(x > PMI{currentsub}.data(cf).HRF.tHRF);
                    echantillon_time = echantillon_time(end);
                    xstop = find_time+step_time;
                    echantillon_timestop = find(xstop > PMI{currentsub}.data(cf).HRF.tHRF);
                    echantillon_timestop = echantillon_timestop(end);
                    labeltime{itime} = [fixdecimal2string(PMI{currentsub}.data(cf).HRF.tHRF(echantillon_time),job.b_option.v_outtimeprecisionleft,job.b_option.v_outtimeprecision),...
                        'to',fixdecimal2string(PMI{currentsub}.data(cf).HRF.tHRF(echantillon_timestop),job.b_option.v_outtimeprecisionleft,job.b_option.v_outtimeprecision)] ;
                    if numel(echantillon_time:echantillon_timestop)>1
                        d1(itime,chok,1)= mean(ndHbO(echantillon_time:echantillon_timestop,:));
                        d1(itime,chok,2)= mean(ndHbR(echantillon_time:echantillon_timestop,:));
                        d1(itime,chok,3)= mean(ndHbT(echantillon_time:echantillon_timestop,:));
                    else
                        d1(itime,chok,1)= ndHbO(echantillon_time:echantillon_timestop,:);
                        d1(itime,chok,2)= ndHbR(echantillon_time:echantillon_timestop,:);
                        d1(itime,chok,3)= ndHbT(echantillon_time:echantillon_timestop,:);
                    end
                    itime = itime + 1;
                end
            end
            %sortie label time pour la sauvegarde et d1(time, ch, HbO, HbR
            %et HbT)

            if isfield(job.typedata.b_HMRdata,'HMR_datad1mat')
                d1 = zeros(1,nbch,3);
                load(job.typedata.b_HMRdata.HMR_datad1mat{ifile})
                d1(1,:,1) = A;
                d1(1,:,2) = A;  
                d1(1,:,3) = A;
            end
            
            label.file = file;
            label.labeltime = labeltime ;
            label.skintype = job.b_option.v_skin;  %0 = SKIN 1=cortex
            angle = job.b_option.v_view;
            toponew = plot_video_HSJHelmet(angle,d1,label,PMI);
            TOPOmat  = [TOPOmat,toponew];
            for i=1:numel(toponew)
                fprintf('%s\n',toponew{i});
            end
            fprintf('%s\n',file_prj_name);
            if job.typedata.b_HMRdata.HMR_savenormfig
                hfigurezscore = figure;
                subplot(2,1,1)
                plot(PMI{currentsub}.data(cf).HRF.tHRF,ndHbO)
                title(['Normalised',typenorm,'HbO',file])
                subplot(2,1,2)
                plot(PMI{currentsub}.data(cf).HRF.tHRF,dHbO)
                title(['HbO',file])
                saveas(hfigurezscore,[label.pathout,file, 'HbO',typenorm,'.jpg'],'jpg');
                saveas(hfigurezscore,[label.pathout,file, 'HbO',typenorm,'.fig'],'fig');
                close(hfigurezscore);
                hfigurezscore = figure;
                subplot(2,1,1)
                plot(PMI{currentsub}.data(cf).HRF.tHRF,ndHbR)
                title(['Normalised',typenorm,'HbR',file])
                subplot(2,1,2)
                plot(PMI{currentsub}.data(cf).HRF.tHRF,dHbR)
                title(['HbR',file])
                saveas(hfigurezscore,[label.pathout,file, 'HbR',typenorm,'.jpg'],'jpg');
                saveas(hfigurezscore,[label.pathout,file, 'HbR',typenorm,'.fig'],'fig');
                close(hfigurezscore);
                hfigurezscore = figure;
                subplot(2,1,1)
                plot(PMI{currentsub}.data(cf).HRF.tHRF,ndHbT)
                title(['Normalised',typenorm,'HbT',file])
                subplot(2,1,2)
                plot(PMI{currentsub}.data(cf).HRF.tHRF,dHbT)
                title(['HbT',file])
                saveas(hfigurezscore,[label.pathout,file, 'HbT',typenorm,'.jpg'],'jpg');
                saveas(hfigurezscore,[label.pathout,file, 'HbT',typenorm,'.fig'],'fig');
                close(hfigurezscore);
            end                     
            clear PMI
            set(guiHOMER,'UserData',[]);
            close(hfig) 
            close(hfighelmet)
        end %Fin fichier hmr
%     catch
%         clear PMI
%         set(guiHOMER,'UserData',[]);
%         close(hfig)  
%     end
end

if isfield(job.typedata,'b_NIRSdata') %NIRS.mat FILE !!!
    if isfield(job.typedata.b_NIRSdata,'NIRSmat')
        if m_prjfile_mode==0 %Paired
            if numel(job.b_option.prjfile)~=numel(job.typedata.b_NIRSdata.NIRSmat)
                msgbox('Please enter same number of .prj than .NIRS data file they should be corresponding if you used Option prj : paired mode');
                return
            end        
        end
        if isfield(job.typedata.b_NIRSdata,'NIRS_dataHemo')
            label.dataHemo = job.typedata.b_NIRSdata.NIRS_dataHemo;
        end
        
        for ifile=1:numel(job.typedata.b_NIRSdata.NIRSmat)
            global currentsub;
            currentsub = 1;
            PMI{currentsub}.currentFile = 1;
            cf = 1;
            if m_prjfile_mode==0 %Paired
                PMI{1}.prj_name = job.b_option.prjfile{ifile};
            elseif m_prjfile_mode==1
                PMI{1}.prj_name = job.b_option.prjfile{1};
            end
            [p,file_prj_name,ext]=fileparts(PMI{1}.prj_name);

            load(job.typedata.b_NIRSdata.NIRSmat{ifile});
            for imodule = numel(NIRS.Dt.fir.pp):-1:1
                if strfind(NIRS.Dt.fir.pp(imodule).pre(1:15),'Epoch averaging')
                    PMI{currentsub}.data(cf).MeasListAct = NIRS.Cf.H.C.ok; %NIRS.Dt.fir.pp(imodule).chok;
                    break
                end
            end
            PMI{currentsub}.data(cf).MeasList = [NIRS.Cf.H.C.id(2:3,:)',...
                ones(size(NIRS.Cf.H.C.id,2),1),...
                [ones(size(NIRS.Cf.H.C.id,2)/2,1);ones(size(NIRS.Cf.H.C.id,2)./2,1).*2]];
            PMI{currentsub}.color = lines(size(NIRS.Cf.H.C.id,2));
            PMI{currentsub}.plotLst = [1];
            PMI{currentsub}.plot = [1,1];
            
            %PRJ
            hfig = figure;
            set(hfig,'visible','off')
            setappdata(0,'gui_SPMvideo',hfig);
            guiHOMER = getappdata(0,'gui_SPMvideo');
            set(guiHOMER,'UserData',PMI); %PMI de base
            
            %Ouverture de IO_Helmet pour affichage
            handles = [];
            option.conc = 1;
            option.scalemin = job.b_option.v_cmin;
            option.scalemax = job.b_option.v_cmax;
            option.scaletr  = job.b_option.v_ctr;
            option.skintype = job.b_option.v_skin;
            option.cover = job.b_option.v_showcover;
            step_time = job.b_option.v_step;
            option.video = 1; % la fonction qui appelle le casque est le mode video
            hfighelmet=IO_HelmetMTG_Display(handles,option);
            set(hfighelmet,'visible','off')
            % LOAD DATA D1
            if job.typedata.b_NIRSdata.NIRS_datatype == 1 ||...%Moyennage
                    job.typedata.b_NIRSdata.NIRS_datatype == 2 %Moyennage avgtime
                [pathstr, name, ext] = fileparts(NIRS.Dt.fir.pp(imodule).p{1});
                xname =fullfile(pathstr,[name,ext]);
            elseif job.typedata.b_NIRSdata.NIRS_datatype == 3||... %Tvalue
                    job.typedata.b_NIRSdata.NIRS_datatype == 4 %Tvalue avgtime
                [pathstr, name, ext] = fileparts(NIRS.Dt.fir.pp(imodule).p{1});
                xname =fullfile(pathstr,[name,'_tval',ext]);
            end
            
            pretime = -abs(str2num(NIRS.Dt.fir.pp(imodule).job.choiceave.pretime));
            PMI{currentsub}.data(cf).HRF.AvgC = fopen_NIR(xname,NIRS.Cf.H.C.N)';
              PMI{currentsub}.data(cf).MeasListAct = NIRS.Cf.H.C.ok;
            %PMI{currentsub}.data(cf).MeasListAct = NIRS.Dt.fir.pp(imodule).chok;
            
            if job.typedata.b_NIRSdata.NIRS_datatype==1 || job.typedata.b_NIRSdata.NIRS_datatype==3
                flag_avgtime = 0; %each time point
            elseif  job.typedata.b_NIRSdata.NIRS_datatype==2 || job.typedata.b_NIRSdata.NIRS_datatype==4
                flag_avgtime = 1; %avg periode of time point
            end
            
            PMI{currentsub}.data(cf).HRF.tHRF = (1:size(PMI{currentsub}.data(cf).HRF.AvgC,1))*1/NIRS.Cf.dev.fs+pretime(1);
            step_time = job.b_option.v_step;
            start = job.b_option.v_start;
            postTime= job.b_option.v_stop-job.b_option.v_start;
            if  postTime < 0
                postTime = 1;
                msgbox('Please look at that the start time it should start before the stop time !')
            end
            
            
            
            idHbO = 1:NIRS.Cf.H.C.N/2;
            idHbR = NIRS.Cf.H.C.N/2+1:NIRS.Cf.H.C.N;
            chokHbO  = find(PMI{currentsub}.data(cf).MeasListAct(1:end/2));
            nbch = numel(PMI{currentsub}.data(cf).MeasListAct(1:end/2));
            chokHbR = chokHbO + nbch;
            if job.typedata.b_NIRSdata.NIRS_datatype == 1 ||...  %dConc
                    job.typedata.b_NIRSdata.NIRS_datatype == 2 ||... %dConc avg time
                    job.typedata.b_NIRSdata.NIRS_datatype == 3 ||... %Tval
                    job.typedata.b_NIRSdata.NIRS_datatype == 4      %Tval avg time
                %No normalisation
                dHbO = PMI{currentsub}.data(cf).HRF.AvgC(:,chokHbO);
                dHbR = PMI{currentsub}.data(cf).HRF.AvgC(:,chokHbR);
                dHbT = PMI{currentsub}.data(cf).HRF.AvgC(:,chokHbO)+PMI{currentsub}.data(cf).HRF.AvgC(:,chokHbR);
                ndHbO = PMI{currentsub}.data(cf).HRF.AvgC(:,chokHbO);
                ndHbR = PMI{currentsub}.data(cf).HRF.AvgC(:,chokHbR);
                ndHbT = PMI{currentsub}.data(cf).HRF.AvgC(:,chokHbO)+PMI{currentsub}.data(cf).HRF.AvgC(:,chokHbR);
                if job.typedata.b_NIRSdata.NIRS_datatype ==1||job.typedata.b_NIRSdata.NIRS_datatype ==2
                    typenorm = 'No mean';
                elseif job.typedata.b_NIRSdata.NIRS_datatype ==3||job.typedata.b_NIRSdata.NIRS_datatype ==4
                    typenorm = 'No Tval';
                end
            end
            
            d1 = zeros(numel(start:step_time:(start+postTime)),NIRS.Cf.H.C.N/2,3);
            itime = 1;          
            if flag_avgtime==0 %Each time point
                for find_time = start:step_time:(start+postTime)
                    x = find_time;
                    echantillon_time = find(find_time > PMI{currentsub}.data(cf).HRF.tHRF);
                    echantillon_time = echantillon_time(end);
                    labeltime{itime} = fixdecimal2string(PMI{currentsub}.data(cf).HRF.tHRF(echantillon_time),job.b_option.v_outtimeprecisionleft,job.b_option.v_outtimeprecision);
                    d1(itime,chokHbO,1) = ndHbO(echantillon_time,:,1);
                    d1(itime,chokHbO,2) = ndHbR(echantillon_time,:,1);
                    d1(itime,chokHbO,3) = ndHbT(echantillon_time,:,1);
                    itime = itime + 1;
                end
            elseif  flag_avgtime==1  %avg time point
                for find_time = start:step_time:(start+postTime)
                    x = find_time;
                    echantillon_time = find(x > PMI{currentsub}.data(cf).HRF.tHRF);
                    echantillon_time = echantillon_time(end);
                    xstop = find_time+step_time;
                    echantillon_timestop = find(xstop > PMI{currentsub}.data(cf).HRF.tHRF);
                    echantillon_timestop = echantillon_timestop(end);
                    labeltime{itime} = [fixdecimal2string(PMI{currentsub}.data(cf).HRF.tHRF(echantillon_time),job.b_option.v_outtimeprecisionleft,job.b_option.v_outtimeprecision),...
                        'to',fixdecimal2string(PMI{currentsub}.data(cf).HRF.tHRF(echantillon_timestop),job.b_option.v_outtimeprecisionleft,job.b_option.v_outtimeprecision)] ;
                    if numel(echantillon_time:echantillon_timestop)>1
                        d1(itime,chokHbO,1)= mean(ndHbO(echantillon_time:echantillon_timestop,:));
                        d1(itime,chokHbO,2)= mean(ndHbR(echantillon_time:echantillon_timestop,:));
                        d1(itime,chokHbO,3)= mean(ndHbT(echantillon_time:echantillon_timestop,:));
                    else
                        d1(itime,chokHbO,1)= ndHbO(echantillon_time:echantillon_timestop,:);
                        d1(itime,chokHbO,2)= ndHbR(echantillon_time:echantillon_timestop,:);
                        d1(itime,chokHbO,3)= ndHbT(echantillon_time:echantillon_timestop,:);
                    end
                    itime = itime + 1;
                end
            end
            
            
            
            [pathavi, name, ext] = fileparts(NIRS.Dt.fir.pp(imodule).p{1});
            file = [fliplr(strtok(fliplr(pathavi),'\')),name];
            [tok,rem]=strtok(job.b_option.v_outpath,':');
            if isempty(rem) %si le nom du repertoire est juste un nom à ajouter dans le fichier de donné
                label.pathout = [pathavi,filesep,job.b_option.v_outpath,filesep];
            else %si le nom du repertoire est un répertoire complet
                try
                    if ~isdir(job.b_option.v_outpath)
                        mkdir(job.b_option.v_outpath)
                    end
                    label.pathout = [job.b_option.v_outpath,filesep];
                catch
                    msgbox(['The name of the output directory : ',job.b_option.v_outpath,' is invalid'])
                end
            end
            
            label.file = file;
            label.labeltime = labeltime ;
            label.skintype = job.b_option.v_skin;  %0 = SKIN 1=cortex
            angle = job.b_option.v_view;
            toponew = plot_video_HSJHelmet(angle,d1,label,PMI);
            TOPOmat  = [TOPOmat,toponew];
            for i=1:numel(toponew)
                fprintf('%s\n',toponew{i});
            end
            fprintf('%s\n',file_prj_name);
            
            if job.typedata.b_NIRSdata.HMR_savenormfig
                hfigurezscore = figure;
                subplot(2,1,1)
                plot(PMI{currentsub}.data(cf).HRF.tHRF,ndHbO)
                title(['Normalised',typenorm,'HbO',file])
                subplot(2,1,2)
                plot(PMI{currentsub}.data(cf).HRF.tHRF,dHbO)
                title(['HbO',file])
                saveas(hfigurezscore,[label.pathout,file, 'HbO',typenorm,'.jpg'],'jpg');
                saveas(hfigurezscore,[label.pathout,file, 'HbO',typenorm,'.fig'],'fig');
                close(hfigurezscore)
                hfigurezscore = figure;
                subplot(2,1,1)
                plot(PMI{currentsub}.data(cf).HRF.tHRF,ndHbR)
                title(['Normalised',typenorm,'HbR',file])
                subplot(2,1,2)
                plot(PMI{currentsub}.data(cf).HRF.tHRF,dHbR)
                title(['HbR',file])
                saveas(hfigurezscore,[label.pathout,file, 'HbR',typenorm,'.jpg'],'jpg');
                saveas(hfigurezscore,[label.pathout,file, 'HbR',typenorm,'.fig'],'fig');
                close(hfigurezscore);
                hfigurezscore = figure;
                subplot(2,1,1)
                plot(PMI{currentsub}.data(cf).HRF.tHRF,ndHbT)
                title(['Normalised',typenorm, 'HbT', file])
                subplot(2,1,2)
                plot(PMI{currentsub}.data(cf).HRF.tHRF,dHbT)
                title(['HbT',file])
                saveas(hfigurezscore,[label.pathout,file, 'HbT',typenorm,'.jpg'],'jpg');
                saveas(hfigurezscore,[label.pathout,file, 'HbT',typenorm,'.fig'],'fig');
                close(hfigurezscore);
            end
            close(hfig)   
            close(hfighelmet)
        end %FIN  liste de FICHIER NIRS.mat
    end
end %Fin NIRS.mat (epoch average data)
if isfield(job.typedata,'b_SPMdata') %SPM analysis NIRS.mat FILE
    if isfield(job.typedata.b_SPMdata,'NIRSmat')
        if m_prjfile_mode==0 %Paired
            if numel(job.b_option.prjfile)~=numel(job.typedata.b_SPMdata.NIRSmat)
            msgbox('Please enter same number of .prj than SPM NIRS.mat data file they should be corresponding if you used Option prj : paired mode');
            end        
        end
        for ifile=1:numel(job.typedata.b_SPMdata.NIRSmat)
            global currentsub;
            currentsub = 1;
            PMI{currentsub}.currentFile = 1;
            cf = 1;
            if m_prjfile_mode==0
                PMI{1}.prj_name = job.b_option.prjfile{ifile};
            elseif m_prjfile_mode==1
                PMI{1}.prj_name = job.b_option.prjfile{1};
            end
            [p,file_prj_name,ext]=fileparts(PMI{1}.prj_name);
            
            load(job.typedata.b_SPMdata.NIRSmat{ifile});
            PMI{currentsub}.data(cf).MeasList = [NIRS.Cf.H.C.id(2:3,:)',...
                ones(size(NIRS.Cf.H.C.id,2),1),...
                [ones(size(NIRS.Cf.H.C.id,2)/2,1);ones(size(NIRS.Cf.H.C.id,2)./2,1).*2]];
            PMI{currentsub}.color = lines(size(NIRS.Cf.H.C.id,2));
            PMI{currentsub}.plotLst = [1];
            PMI{currentsub}.plot = [1,1];
            PMI{currentsub}.data(cf).MeasListAct = ones(size(NIRS.Cf.H.C.id,2),1);
            %PRJ
            hfig = figure;
            setappdata(0,'gui_SPMvideo',hfig);
            guiHOMER = getappdata(0,'gui_SPMvideo');
            set(guiHOMER,'UserData',PMI); %PMI de base
            %Ouverture de IO_Helmet pour affichage
            handles = [];
            option.conc = 1;
            option.scalemin = job.b_option.v_cmin;
            option.scalemax = job.b_option.v_cmax;
            option.scaletr  = job.b_option.v_ctr;
            step_time = job.b_option.v_step;
            option.skintype = job.b_option.v_skin;
            option.cover = job.b_option.v_showcover;
            option.video = 1; % la fonction qui appelle le casque est le mode video
            hfighelmet=IO_HelmetMTG_Display(handles,option);
            set(hfighelmet,'visible','off')
            %LOAD SPM RESULT
            if isfield(NIRS,'SPM')
                load([NIRS.SPM{1},filesep,'SPM.mat'],'-mat');
            else
                msgbox('No SPM-GLM estimation have been process for this NIRS.mat file, please review your selection maybe is it in the dataSPM folder')
                return
            end
            for indbeta = 1:1%size(SPM.xXn{1}.beta,1)
                Ball = [];
                Tall = [];
                for ibloc = 1:numel(SPM.xXn)
                    Tall = [Tall;SPM.xXn{ibloc}.t(indbeta ,:)];
                    Ball = [Ball;SPM.xXn{ibloc}.beta(indbeta,:)];
                end
            end
            %DERIVATIVE ? NO
            nbch = numel(NIRS.Cf.H.C.n);
            half = numel(NIRS.Cf.H.C.n)/2;
            %tavgHbO.mat
            d1(:,:,1) = mean(Ball(:,1:half))./(std(Ball(:,1:half))./sqrt(size(Ball,1)));
            d1(:,:,2) = mean(Ball(:,half+1:2*half))./(std(Ball(:,half+1:2*half))./sqrt(size(Ball,1)));
            d1(:,:,3) =  mean(Ball(:,1:half))./(std(Ball(:,1:half))./sqrt(size(Ball,1)))+...
                mean(Ball(:,half+1:2*half))./(std(Ball(:,half+1:2*half))./sqrt(size(Ball,1)));
            
            
            [pathname, name, ext] = fileparts(NIRS.Dt.fir.pp(1).p{1});
            file = [name];
            [tok,rem]=strtok(job.b_option.v_outpath,':');
            if isempty(rem) %si le nom du repertoire est juste un nom à ajouter dans le fichier de donné
                label.pathout = [NIRS.SPM{1},filesep,job.b_option.v_outpath,filesep];
            else %si le nom du repertoire est un répertoire complet
                try
                    if ~isdir(job.b_option.v_outpath)
                        mkdir(job.b_option.v_outpath);
                    end
                    label.pathout = job.b_option.v_outpath;
                catch
                    msgbox(['The name of the output directory : ',job.b_option.v_outpath,' is invalid'])
                end
            end
            label.file = file;
            label.labeltime = {'T'};%labeltime ;
            label.skintype = job.b_option.v_skin;  %0 = SKIN 1=cortex
            angle = job.b_option.v_view;
            toponew = plot_video_HSJHelmet(angle,d1,label,PMI);
            TOPOmat  = [TOPOmat,toponew];
            for i=1:numel(toponew)
                fprintf('%s\n',toponew{i});
            end
            fprintf('%s\n',file_prj_name);
            close(hfig);
            close(hfighelmet);
        end
    end
end
if isfield(job.typedata,'b_vColordata')
    TOPOmat = [];

    [pathavi, name, exttest] = fileparts(job.typedata.b_vColordata.vColordata{1});
    label.vcolortype  = 1;
    if strcmp(exttest,'.mat')
        loopoverTOPOmat = numel(job.typedata.b_vColordata.vColordata); %Loop over file of topo.mat
        mode_TOPOmat = 1; %mode TOPOmat file    
        if m_prjfile_mode==0 %Paired
            if numel(job.b_option.prjfile)~=numel(job.typedata.b_vColordata.vColordata)
                msgbox('Usualy you can use the same project for all Topo.mat file. ''Option prj'' : Same .prj for all file''',...
                    'But paired is selected now and you don''t have the same number of subject.')
                return
            end        
        end
    else
        loopoverTOPOmat = 1; %Loop over file .img
        mode_TOPOmat = 0; %mode Vcolor file
    end
    
    for iTOPOmat=1:loopoverTOPOmat
        global currentsub
        currentsub = 1;
        PMI{1}.currentFile = 1;
        if m_prjfile_mode==1||mode_TOPOmat==0
            PMI{1}.prj_name = job.b_option.prjfile{1}; 
        elseif m_prjfile_mode==0
            PMI{1}.prj_name = job.b_option.prjfile{iTOPOmat};
        end
        [p,file_prj_name,ext]=fileparts(PMI{1}.prj_name);
        hfig = figure;
        set(hfig,'visible','off')
        setappdata(0,'gui_SPMvideo',hfig);
        guiHOMER = getappdata(0,'gui_SPMvideo');
        set(guiHOMER,'UserData',PMI); %PMI de base           
        %Ouverture de IO_Helmet pour affichage
        handles = [];
        option.conc = 1;
        option.scalemin = job.b_option.v_cmin;
        option.scalemax = job.b_option.v_cmax;
        option.scaletr  = job.b_option.v_ctr;
        option.skintype = job.b_option.v_skin;
        option.cover = job.b_option.v_showcover;
        step_time = job.b_option.v_step;
        option.video = 1; % la fonction qui appelle le casque est le mode video
        hfighelmet = IO_HelmetMTG_Display(handles,option);

        if strcmp(exttest,'.mat')
            [pathavi, name, exttest] = fileparts(job.typedata.b_vColordata.vColordata{iTOPOmat});
            load(job.typedata.b_vColordata.vColordata{iTOPOmat});
            filelist_vColordata = TOPO.pp(end).p;
            label.file = name; 
            toponew = job.typedata.b_vColordata.vColordata{iTOPOmat};%Dans ce cas aucun topomat n'est créer on utilise simplement l'ancient ou aucun si vcolor file
            fprintf('%s\n',toponew);                
            fprintf('%s\n',file_prj_name);
        else 
            [pathavi, name, exttest] = fileparts(job.typedata.b_vColordata.vColordata{iTOPOmat});
            filelist_vColordata = job.typedata.b_vColordata.vColordata;   
            for i=1:numel(filelist_vColordata)
                fprintf('%s\n',filelist_vColordata{i});
            end
            fprintf('%s\n',file_prj_name);
            label.file = 'ManualSelection';
        end
        

        file = [fliplr(strtok(fliplr(pathavi),filesep)),name];
        [tok,rem]=strtok(job.b_option.v_outpath,':');
        if isempty(rem) %si le nom du repertoire est juste un nom à ajouter dans le fichier de donné
            label.pathout = [pathavi,'/',job.b_option.v_outpath,filesep];
        else %si le nom du repertoire est un répertoire complet
            try
                if ~isdir(job.b_option.v_outpath)
                    mkdir(job.b_option.v_outpath);
                end
                label.pathout = [job.b_option.v_outpath,filesep];
            catch
                msgbox(['The name of the output directory : ',job.b_option.v_outpath,' is invalid']);
            end
        end
        [pathavi, name, ext] = fileparts(filelist_vColordata{1});
        label.skintype = job.b_option.v_skin;
        angle = job.b_option.v_view;
        PMI = [];
        for ifile=1:numel(filelist_vColordata)
            filename = filelist_vColordata{ifile};
            time = strfind(filename,'(s)');
            if isempty(time)
                label.labeltime{ifile} = ['s',num2str(ifile),'-na'];                
            else
                label.labeltime{ifile} = ['s',num2str(ifile),'-',filename(time-8:time-1)];
            end
            vcolor = opentopo(filename)';
            indna = isnan(vcolor);
            if ~isempty(indna)
                vcolor(indna)=0;
            end
            vcolortemp = vcolor;
            try
            d1(ifile,:) = vcolor;
            catch
                msgbox('Check input file they should have the same format')
            end
        end
        plot_video_HSJHelmet(angle,d1,label,PMI);        
        a_cote = [];
        a_face = [];
        rem=angle;
        i=1;
        while ~isempty(rem)
            i=i+1;
            [token, rem] = strtok(rem,',');
            if isempty(rem)
                a_cote =[a_cote 90];
                a_face = [a_face,0];
            else
                a_cote = [a_cote,str2num(token)];
                [token, rem] = strtok(rem(2:end),';');
                a_face = [a_face,str2num(token)];
            end
        end
        for idaxes = 1:numel(a_cote)
            definefirsttime  = 0;
            Total = [];
            for itime=1:numel(label.labeltime)
                labelangle = [num2str(a_cote(idaxes)) ,'_', num2str(a_face(idaxes))];
                pathavi = label.pathout;
                n = label.labeltime{itime};
                filename= [pathavi,labelangle,filesep,labelangle,label.file, 'Time', label.labeltime{itime},'(s)','.tif'];
                [X] = imread(filename,'tif');
                if definefirsttime==0
                    X2d = sum(X,3);
                    idx = find(765~=sum(X2d,1)./size(X2d,1));
                    idy = find(765~=sum(X2d,2)./size(X2d,2));
                    idxtrim = (idx(1)-30):(idx(end)+30);
                    idytrim = (idy(1)-30):(idy(end)+30);
                end
                
                hfigurestd = figure;
                set(gca,'visible','off');
                text(0.5,0.5,[label.labeltime{itime},' s'],'FontSize',30,'HorizontalAlignment','left','VerticalAlignment','bottom','FontName','arial','unit','normalized');
                saveas(hfigurestd,'temp.tif','tif');
                imgall = imread(['temp.tif'],'tif');
                close(hfigurestd);
                X2d = sum(imgall,3);
                if definefirsttime==0
                    idxch = find(765~=sum(X2d,1)./size(X2d,1));
                    idych = find(765~=sum(X2d,2)./size(X2d,2));
                    idxchok = (idxch(1)-120):(idxch(end)+120);
                    idychok  = (idych(1)):(idych(end));
                    definefirsttime = 1;
                end
                titrenb = imgall(idychok, idxchok,:);
                img=imageresize(titrenb,[],numel(idxtrim));
                Trim = [img;X(idytrim,idxtrim,:)];
                Total = [Total,Trim];                
            end
            [pathavi,'all',labelangle,label.file,label.labeltime{1},'to',label.labeltime{end},'.tif']
            imwrite(Total,[pathavi,'all',labelangle,label.file,label.labeltime{1},'to',label.labeltime{end},'.tif'],'tif');
        end   
        close(hfig)
        close(hfighelmet)
    end
    if mode_TOPOmat==0
        fprintf('%s\n','No dependancy available for img file');
    end
    
 end
    out.TOPOmat = TOPOmat;    
end
