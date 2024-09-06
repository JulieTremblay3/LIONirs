function TOPOmat = plot_video_HSJHelmet(varargin)
%angle % string contenant les angles de vue à afficher 
% Pour appeler cette fonction le programme IO_HelmetMTG_Display doit au
% préalable être ouvert avec les options d'affichage voulue préselectionné,
% correspondant 
warning off all
if numel(varargin) >=1
    angle = varargin{1};%90;30,-90;30,50;20
end
if numel(varargin) >=2
    d1 = varargin{2};
end
if numel(varargin) >=3
    label = varargin{3};
end
if numel(varargin) >=4
    PMI = varargin{4}; %BASIC PMI structure 
end
%test essenstial field if absent define them
try;label.vcolortype;catch;label.vcolortype = 0;end
try;label.projectiontype; catch; label.projectiontype = 0; end
try;label.dataHemo;catch; label.dataHemo = 1; end
%non de sortir 
labeltime = label.labeltime;
file = label.file;
skintype = label.skintype;
% label.pathout 
% label.file 
% [pathout,file, 'HbTTime   ', labeltime{itime},'(s)','.vcolor']
%NAME FILE ABSOLUT RÉFÉRENCE AU NOM DE FICHIER DU SUJET 

if ~isdir([label.pathout])
     mkdir([label.pathout]);
 end

%try
    %Separation des angles de vue ! 
            a_cote = [];
            a_face = [];
            rem=angle;
            if isempty(rem)
                flag_viewangle = 0;
            end
            i=1;
            while ~isempty(rem)
                i=i+1;
                [token, rem] = strtok(rem,',');
                if isempty(rem)
                    a_cote =[a_cote 90];
                    a_face = [a_face,0];
                    flag_viewangle = 0;
                else
                    a_cote = [a_cote,str2num(token)];
                    [token, rem] = strtok(rem(2:end),';');
                    a_face = [a_face,str2num(token)];
                    flag_viewangle = 1;
                end
            end
      %Fin separation des angles de vue
      
      
      %Section ouvrir le video
      %FAIT POUR VIEILLE VERSION 2009 
%         if file~=0
%             try
%                 if flag_viewangle
%                     mov = avifile([label.pathout,filesep,file,'.avi'],'fps',1);
%                 end
%                 nameavi =[file,'.avi'];
%                 pathavi = [label.pathout,filesep];
%             catch
%                 [nameavi,pathavi]=uiputfile([label.pathout,file,'.avi']);%si nom invalide selection manuel du nom
%                 if flag_viewangle
%                     mov = avifile([pathavi,nameavi],'fps',1);
%                 end
%             end
%         end
        nameavi =[file,'.avi'];
        pathavi = [label.pathout,filesep];
    if label.vcolortype == 0 %Travail data original pas avec les fichiers déjà en topo
        if ~isdir([pathavi,'Topo',filesep])
            mkdir([pathavi,'Topo',filesep]);
        end
        if label.dataHemo==1 | label.dataHemo==2
        if ~isdir([pathavi,'Topo',filesep,'HbO',filesep])
            mkdir([pathavi,'Topo',filesep,'HbO',filesep]);
        end
        end
        if label.dataHemo==1 | label.dataHemo==3
        if ~isdir([pathavi,'Topo',filesep,'HbR',filesep])
            mkdir([pathavi,'Topo',filesep,'HbR',filesep])
        end
        end
        if 0 %label.dataHemo==3
        if ~isdir([pathavi,'Topo',filesep,'HbT',filesep])
            mkdir([pathavi,'Topo',filesep,'HbT',filesep])
        end
        end
    end
        
      
      guiHOMER = getappdata(0,'gui_SPMvideo');
      guiHelmet = getappdata(0,'guiHelmet');
        if isempty(guiHelmet)
            msgbox('Sorry, you should open the helmet .prj to enable this option')
            return
        end   
           
            %% GRAPH 1
            if flag_viewangle 
            numwindows = numel(a_face);
            deltay = 0.1; 
            xaxiszoom = 0.13;
            yaxiszoom = 0.13; 
            zaxiszoom = 0.13;
            speed = 1; 
            xzoom = 1.2;
            haut = 0.45;
            large = 1/(numwindows);
            ydistlabel = 0.2;
            col = 1/numwindows - deltay;
            ini_col = 0.05;
            dist = 0.3;
            haxes = [];
            ini_row = 0.05;         % image du bas
            row = 0.45;
            hfiguredisplay=figure;
            set(hfiguredisplay,'visible','off')
            set(hfiguredisplay,'unit', 'normalized','position',[0,0,1,1],'color',[1,1,1])
            hback = axes;
            set(hback ,'units','normalize','position',[ini_col,ini_row,large*numwindows-deltay,haut*2])
            set(hback,'xtick',[],'ytick',[],'ztick',[])
            resetview_axes(hback);
            xlim = get(hback,'clim');
            delete(hback);
            hback = axes;
            set(hback ,'units','normalize','position',[ini_col,ini_row,0.8,haut*2]);
            set(hback,'xtick',[],'ytick',[],'ztick',[], 'xcolor',[1,1,1],'ycolor',[1,1,1],'zcolor',[1,1,1]);
            caxis(xlim)
          
            if 0
            colorbar
            saveas(hfiguredisplay,[pathavi,'colorbar.tif'],'tif')
            imgall = imread([pathavi,'colorbar.tif'],'tif'); 
            X2d = sum(imgall,3);
            idx = find(765~=sum(X2d,1)./size(X2d,1));
            idy = find(765~=sum(X2d,2)./size(X2d,2));
            idytrim = (idy(1)):(idy(end));
            idxtrim = (idx(1)-40):(idx(end)+40);
            colorbartif = imgall(idytrim,idxtrim,:);
            imwrite(colorbartif,[pathavi,'colorbar.tif'],'tif');
            clear colorbartif imgall
            clear idx idy 
            end
            
            %Orientation of the topo
            for i=0:numwindows-1
                hsubaxes = axes;%HbO
                haxes=[haxes,hsubaxes];
                angle_cote = a_cote(i+1);
                angle_face = a_face(i+1);
                cp_x = cos(angle_cote*pi/180) * cos(angle_face*pi/180)*dist;
                cp_y = sin(angle_cote*pi/180)* cos(angle_face*pi/180)*dist;
                cp_z = sin(angle_face*pi/180)*dist;
                set(hsubaxes,'CameraTarget',[0,0,0]);
                set(hsubaxes,'CameraPosition',[cp_x,cp_y,cp_z]);
                set(hsubaxes,'CameraUpVector',[0,0,1]);
                set(hsubaxes,'UserData','Movies');
                set(hsubaxes,'units','normalize','position',[ini_col+col*i,ini_row,large,haut]);
                set(hsubaxes,'xtick',[],'ytick',[],'ztick',[])
                axis([-1*xaxiszoom,1*xaxiszoom, -1*yaxiszoom,1*yaxiszoom,-1*zaxiszoom,1*zaxiszoom])                
                hsubaxes = axes;%HbR
                haxes=[haxes,hsubaxes];
                set(hsubaxes,'CameraTarget',[0,0,0]);
                set(hsubaxes,'CameraPosition',[cp_x,cp_y,cp_z]);
                set(hsubaxes,'CameraUpVector',[0,0,1]);
                set(hsubaxes,'units','normalize','position',[ini_col+col*i,row,large,haut]);
                set(hsubaxes,'UserData','Movies');
                axis([-1*xaxiszoom,1*xaxiszoom, -1*yaxiszoom,1*yaxiszoom,-1*zaxiszoom,1*zaxiszoom])
                set(hsubaxes,'xtick',[],'ytick',[],'ztick',[])
                
            end
            end
            guiHelmet = getappdata(0,'guiHelmet');
            PrjStruct = getappdata(guiHelmet,'PrjStruct');
            
            %Affichage structure d1 , labeltime 
            %sauvegarde du fichier et d1(time, ch, HbO, HbR et HbT)
            
            for itime = 1:numel(labeltime)
                %Affichage titre
                if label.vcolortype == 1
                    hHbO=uicontrol('Style','text','string',[label.file],'units','normalize',...
                        'Position',[0.01,0.9-.45/2,0.04 0.05],'FontSize',12,'tag','edittimingHbO','BackgroundColor',[1,1,1],'HorizontalAlignment','center');
                    hTime=uicontrol('Style','text','string',['Time   ', labeltime{itime},'(s)'],'units','normalize',...
                        'Position',[0.01,0.85,0.1 0.1],'FontSize',12,'tag','edittimingHbO','BackgroundColor',[1,1,1],'HorizontalAlignment','center');

                    hHbR = uicontrol('Style','text','string',[''],'units','normalize',...
                        'Position',[0.01,0.45-.45/2,0.04 0.05],'FontSize',12,'tag','edittimingHbO','BackgroundColor',[1,1,1],'HorizontalAlignment','center');
                else
                     hHbO=uicontrol('Style','text','string',['HbO'],'units','normalize',...
                    'Position',[0.01,0.9-.45/2,0.04 0.05],'FontSize',12,'tag','edittimingHbO','BackgroundColor',[1,1,1],'HorizontalAlignment','center');
                    hTime=uicontrol('Style','text','string',['Time   ', labeltime{itime},'(s)'],'units','normalize',...
                    'Position',[0.01,0.85,0.1 0.1],'FontSize',12,'tag','edittimingHbO','BackgroundColor',[1,1,1],'HorizontalAlignment','center');
                    hHbR = uicontrol('Style','text','string',['HbR'],'units','normalize',...
                    'Position',[0.01,0.45-.45/2,0.04 0.05],'FontSize',12,'tag','edittimingHbO','BackgroundColor',[1,1,1],'HorizontalAlignment','center');
                end
                %%%%%%%Début HbO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if  label.dataHemo==1|label.dataHemo==2
                if label.vcolortype == 0;
                    if label.projectiontype==0
                        PrjStruct = display_MRIcolor(PrjStruct,PMI,d1(itime,:,1),skintype);
                    elseif label.projectiontype==1
                    	PrjStruct = display_MRIcolor_InverseWeightedDistance(PrjStruct,PMI,d1(itime,:,1),skintype);
                    end
                end
                if file~=0
                    oMRI = get_MRI_Data(PrjStruct);
                    if label.vcolortype == 1
                        vColor = d1(itime,:);
                        if skintype == 0
                            [VertexBuffer, IndexBuffer] = get_SkinMesh( oMRI );
                            vColor = d1(itime,1:size(VertexBuffer,1));
                            oMRI = set_SkinVcolor(oMRI,vColor');                           
                        elseif skintype == 1
                            [VertexBuffer, IndexBuffer] = get_CortexMeshLowRes( oMRI );         
                            vColor = d1(itime,1:size(VertexBuffer,1));
                            oMRI = set_CortexLowResVcolor(oMRI,vColor');
                        end
                        PrjStruct = set_MRI_Data(PrjStruct,oMRI);
                    else   
                        if skintype==0 %skin
                            vColor = get_SkinVcolor(oMRI);
                        elseif skintype==1 %Cortex
                            vColor = get_CortexLowResVcolor(oMRI);
                        end
                         savetopo([pathavi,'Topo/HbO/',file, 'HbOTime', labeltime{itime},'(s)','.img'],vColor,1);
                         TopoHbO.prj{1} = PMI{1}.prj_name;
                         TopoHbO.pp(1).pre = 'videoHbO';
                         TopoHbO.pp(1).p{itime} = [pathavi,'Topo/HbO/',file, 'HbOTime', labeltime{itime},'(s)','.vColor'];
                         TopoHbO.pp(1).srxtype(itime)=skintype;
                    end                    
                end
                setappdata(guiHelmet,'PrjStruct',PrjStruct);             
                idaxes=1;
                if flag_viewangle 
                    for i = 2:2:numel(haxes)
                        resetview_axes(haxes(i));    
                        oldpos = get(haxes(i),'position');
                        %Affichage full screen sans rien d'autre
                        fullpos = [0,0,1,1];
                        hobject = get(hfiguredisplay,'children'); 
                        set(hobject,'visible','off');
                        himageall = [];
                        for j = 1:numel(haxes)
                            if j~=i
                                himage = get(haxes(j),'children');
                                himageall = [himageall,himage];                        
                            end                       
                        end
                        if ~isempty(himageall)
                            set(himageall,'visible','off');
                        end
                        set(haxes(i),'visible','on');
                        set(haxes(i),'position',fullpos);
                        labelangle = [num2str(a_cote(idaxes)) ,'_', num2str(a_face(idaxes))];
                        idaxes = idaxes+1;
                        if ~isdir([pathavi,labelangle])
                            mkdir([pathavi,labelangle]);
                        end
                        if label.vcolortype == 0;
                        if ~isdir([pathavi,labelangle,filesep,'HbO'])
                            mkdir([pathavi,labelangle,filesep,'HbO'])
                        end
                        end
                        if label.vcolortype == 0;                            
                            saveas(hfiguredisplay,[pathavi,labelangle,filesep,'HbO',filesep,labelangle,file, 'HbOTime', labeltime{itime},'(s)','.tif'],'tif');
                        else
                            saveas(hfiguredisplay,[pathavi,labelangle,filesep,labelangle,file, 'Time', labeltime{itime},'(s)','.tif'],'tif');
                        end
                        %Affichage retour à la vue film avec toute les images;
                        set(haxes(i),'position',oldpos)
                        if ~isempty(himageall)
                            set(himageall,'visible','on');
                        end
                        set(hobject,'visible','on');               
                    end
                    if label.vcolortype == 1 %Enlever les axes initialement prévu pour HbR
                          for i = 1:2:numel(haxes)
                              set(haxes(i),'visible','off');
                          end
                    end
                end
            end
                %%%%%%%%%%%%FIN HbO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                %Début HbT
                 if 0 %label.dataHemo==3
                if label.vcolortype == 0
                
                    if label.projectiontype==0
                        PrjStruct = display_MRIcolor(PrjStruct,PMI,d1(itime,:,1),skintype);
                    elseif label.projectiontype==1
                    	PrjStruct = display_MRIcolor_InverseWeightedDistance(PrjStruct,PMI,d1(itime,:,1),skintype);
                    end
                    
                    if file~=0
                        oMRI = get_MRI_Data(PrjStruct);
                        if skintype==0 %skin
                            vColor = get_SkinVcolor(oMRI);    
                        elseif skintype==1 %Cortex
                            vColor = get_CortexLowResVcolor(oMRI);
                        end
                        filename = [pathavi,'Topo',filesep,'HbT',filesep,file,  'HbTTime', labeltime{itime},'(s)','.vcolor'];
                        savetopo(filename,vColor,1)
                        TopoHbT.prj{1} = PMI{1}.prj_name;
                        TopoHbT.pp(1).pre = 'videoHbT';
                        TopoHbT.pp(1).p{itime} = [pathavi,'Topo',filesep,'HbT',filesep,file, 'HbTTime', labeltime{itime},'(s)','.vColor'];
                        TopoHbT.pp(1).srxtype(itime)=skintype;
                    end                    
                    if flag_viewangle 
                    idaxes=1;
                    setappdata(guiHelmet,'PrjStruct',PrjStruct);
                    for i = 1:2:numel(haxes)
                        resetview_axes(haxes(i));
                        %save image .tif
                        oldpos = get(haxes(i),'position');
                        %Affichage full screen sans rien d'autre
                        fullpos = [0,0,1,1];
                        hobject = get(hfiguredisplay,'children'); 
                        set(hobject,'visible','off');
                        himageall = [];
                        for j = 1:numel(haxes)
                            if j~=i
                                himage = get(haxes(j),'children');
                                himageall = [himageall,himage];                        
                            end                       
                        end
                        if ~isempty(himageall)
                            set(himageall,'visible','off');
                        end
                        set(haxes(i),'visible','on');
                        set(haxes(i),'position',fullpos);
                        labelangle = [num2str(a_cote(idaxes)) ,'_', num2str(a_face(idaxes))];
                        idaxes = idaxes+1;
                        if ~isdir([pathavi,labelangle,filesep,'HbT'])
                            mkdir([pathavi,labelangle,filesep,'HbT'])
                        end
                        saveas(hfiguredisplay,[pathavi,labelangle,filesep,'HbT',filesep,labelangle,file, 'HbTTime', labeltime{itime},'(s)','.tif'],'tif');
                        %Affichage retour à la vue film avec toute les images;
                        set(haxes(i),'position',oldpos)
                        if ~isempty(himageall)
                            set(himageall,'visible','on');
                        end
                        set(hobject,'visible','on');   
                    end
                    end
                end
            end
                %%%%%%%%%Fin HbT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                %Début HbR
                if label.dataHemo==1 | label.dataHemo==3
                if label.vcolortype == 0
                      if label.projectiontype==0
                        PrjStruct = display_MRIcolor(PrjStruct,PMI,d1(itime,:,2),skintype);
                      elseif label.projectiontype==1
                    	PrjStruct = display_MRIcolor_InverseWeightedDistance(PrjStruct,PMI,d1(itime,:,2),skintype);
                    end
                    %PrjStruct = display_MRIcolor(PrjStruct,PMI,d1(itime,:,2),skintype);
                    if file~=0
                        oMRI = get_MRI_Data(PrjStruct);
                        if skintype==0 %skin
                            vColor = get_SkinVcolor(oMRI);
                        elseif skintype==1 %Cortex
                            vColor = get_CortexLowResVcolor(oMRI);
                        end
                        filename = [pathavi,'Topo',filesep,'HbR',filesep,file, 'HbRTime', labeltime{itime},'(s)','.vcolor'];
                        savetopo(filename,vColor,1)
                        TopoHbR.prj{1} = PMI{1}.prj_name;
                        TopoHbR.pp(1).pre = 'videoHbR';
                        TopoHbR.pp(1).p{itime} = [pathavi,'Topo',filesep,'HbR',filesep,file, 'HbRTime', labeltime{itime},'(s)','.vColor'];
                        TopoHbR.pp(1).srxtype(itime)=skintype;
                    end
                    setappdata(guiHelmet,'PrjStruct',PrjStruct);
                    idaxes=1;
                    if flag_viewangle
                    for i = 1:2:numel(haxes)
                        resetview_axes(haxes(i));
                        %save image .tif
                        oldpos = get(haxes(i),'position');
                        %Affichage full screen sans rien d'autre
                        fullpos = [0,0,1,1];
                        hobject = get(hfiguredisplay,'children'); 
                        set(hobject,'visible','off');
                        himageall = [];
                        for j = 1:numel(haxes)
                            if j~=i
                                himage = get(haxes(j),'children');
                                himageall = [himageall,himage];                        
                            end                       
                        end
                        if ~isempty(himageall)
                            set(himageall,'visible','off');
                        end
                        set(haxes(i),'visible','on');
                        set(haxes(i),'position',fullpos);
                        labelangle = [num2str(a_cote(idaxes)) ,'_', num2str(a_face(idaxes))];
                        idaxes = idaxes+1;
                        if ~isdir([pathavi,labelangle,filesep,'HbR'])
                            mkdir([pathavi,labelangle,filesep,'HbR']);
                        end
                        saveas(hfiguredisplay,[pathavi,labelangle,filesep,'HbR',filesep,labelangle,file, 'HbRTime', labeltime{itime},'(s)','.tif'],'tif');
                        %Affichage retour à la vue film avec toute les images;
                        set(haxes(i),'position',oldpos)
                        if ~isempty(himageall)
                            set(himageall,'visible','on');
                        end
                        set(hobject,'visible','on');   
                    end
                    end
                end
            end
                    %FIN HbR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  if flag_viewangle 
                    uiwait(hfiguredisplay,0.5)
                    if file~=0
                        %F = getframe(close);
                        figure(hfiguredisplay)
                       % F = getframe %(hfiguredisplay);
                       %old version
%                         F.cdata = screencapture(hfiguredisplay)
%                         F.colormap = [];                        
%                         mov = addframe(mov,F);
                        %new version
                    end 
                    v.frames(itime) = getframe(hfiguredisplay);
                    v.times(itime) = itime/ 1; %30; % display at 1fps
                    disp(['Frame: ', num2str( itime),' Time: ',labeltime{itime}]);
                    v.width=size(v.frames(1).cdata,2);
                    v.height=size(v.frames(1).cdata,1);
                    uiwait(gcf,speed)
                    delete(hHbO)
                    delete(hHbR)
                    delete(hTime) 
                  end                
            end
%             nameavi =[file,'.avi'];
%         pathavi = [label.pathout,filesep];
            nameavi = 'video.wmv';
            disp(['Create video file: ', fullfile([pathavi,nameavi])])         
            mmwrite([fullfile([pathavi,nameavi])],v);
            
            if label.vcolortype == 0; %sortie pour tout les modes sauf le .vcolor
                if label.dataHemo==1|label.dataHemo==2
                    [dir,fil,ext] = fileparts(TopoHbO.pp(1).p{1});
                    TOPO = TopoHbO;
                    save(fullfile(dir,[file,'HbOTOPO.mat']),'TOPO','-mat');
                    TOPOmat{1} = fullfile(dir,[file,'HbOTOPO.mat']);
                end
              if label.dataHemo==1|label.dataHemo==3
                [dir,fil,ext] = fileparts(TopoHbR.pp(1).p{1});
                TOPO = TopoHbR;
                save(fullfile(dir,[file,'HbRTOPO.mat']),'TOPO','-mat');
                 TOPOmat{2} = fullfile(dir,[file,'HbRTOPO.mat']);
              end
              if 0 % label.dataHemo==3
                [dir,fil,ext] = fileparts(TopoHbT.pp(1).p{1});
                TOPO = TopoHbT;
                save(fullfile(dir,[file,'HbTTOPO.mat']),'TOPO','-mat');
                TOPOmat{3} = fullfile(dir,[file,'HbTTOPO.mat']);
              end
          else
              TOPOmat = [];
          end
          
          
%         catch exception
%             disp(exception.message);
%             close(mov)
%         end
%         if flag_viewangle
%             if file~=0
%                 mov = close(mov);
%             end 
%             close(hfiguredisplay);
%             try
%                 [video, audio] = mmread([pathavi,nameavi]);
%                 filename = [pathavi,nameavi];
%                 mmwrite([filename(1:end-3),'wmv'],audio,video);
%             catch
%                 disp('Video format was not convert in wmv')
%             end
%         end
        %FIN FONCTION 
        warning on all
end
