function varargout = GUI_ERLCOH_EEG(varargin)
% GUI_ERLCOH_EEG M-file for GUI_ERLCOH_EEG.fig
%      GUI_ERLCOH_EEG, by itself, creates a new GUI_ERLCOH_EEG or raises the existing
%      singleton*.
%
%      H = GUI_ERLCOH_EEG returns the handle to a new GUI_ERLCOH_EEG or the handle to
%      the existing singleton*.
%
%      GUI_ERLCOH_EEG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_ERLCOH_EEG.M with the given input arguments.
%
%      GUI_ERLCOH_EEG('Property','Value',...) creates a new GUI_ERLCOH_EEG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_ERLCOH_EEG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_ERLCOH_EEG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_ERLCOH_EEG

% Last Modified by GUIDE v2.5 02-Oct-2019 10:33:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_ERLCOH_EEG_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_ERLCOH_EEG_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI_ERLCOH_EEG is made visible.
function GUI_ERLCOH_EEG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_ERLCOH_EEG (see VARARGIN)

% Choose default command line output for GUI_ERLCOH_EEG
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_ERLCOH_EEG wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_ERLCOH_EEG_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;








% --- Executes on button press in btn_openEEG.
function btn_openEEG_Callback(hObject, eventdata, handles)
% hObject    handle to btn_openEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path]=uigetfile('.dat','MultiSelect','on')
set(handles.edit_pathEEG,'string',[path])
set(handles.listbox_fileEEG,'string',[file])
set(handles.listbox_fileEEG,'value',1)

% --- Executes on button press in btn_ERLCOHcluster.
function btn_ERLCOHcluster_Callback(hObject, eventdata, handles)
% hObject    handle to btn_ERLCOHcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%load('C:\data\BCVareta\Barbados\Barbados\data\mycolormap_brain_basic_conn.mat');
%load('C:\data\BCVareta\Barbados\Barbados\data\Barbados_DATA_6k_15k_INFO.mat');

        mfunction = mfilename('fullpath');
        [pathfct,filefct,ext]= fileparts(mfunction);
        load([pathfct,'\mycolormap_brain_basic_conn.mat']);
        load([pathfct,'\Barbados_DATA_6k_15k_INFO.mat']);
       % filein = {get(handles.edit_fileEEG,'string')};
        elefile = get(handles.edit_elefile,'string');
        Fs = str2num(get(handles.edit_resample,'string'));
        
        parameter.Fs = Fs;
        parameter.layer = str2num(get(handles.edit_frequence,'string'));
        parameter.width = str2num(get(handles.edit_width,'string'));
        parameter.baseline = str2num(get(handles.edit_Baseline,'string'));
        parameter.ERLCOHcut = str2num(get(handles.edit_ERLCOHcut,'string'));
        parameter.maxcluster = str2num(get(handles.edit_MaxCluster,'string'));
        pathdata = get(handles.edit_pathEEG,'string');
        filedata = get(handles.listbox_fileEEG,'string');
        if ~iscell(filedata)
            filedata = {filedata}
        end
        for ifile=1:numel( filedata)
         filenametot =  filedata{ifile};
         filenametot =filenametot(1:end-4);
        indlist =findelevhdr([elefile],fullfile(pathdata,[filenametot,'.vhdr' ]));
        [nameele,x,y,z] = readelefile(elefile);
        fid = fopen(fullfile(pathdata,[filenametot,'.dat']),'r');
        d1= fread(fid,'float32','ieee-le');
        fclose(fid);
        info = read_vhdr(fullfile(pathdata,[filenametot,'.vhdr' ]));
        tmp= reshape(d1,info.DataPoints,info.NumberOfChannels); 
        tmp = tmp(:,indlist);
        %RESAMPLE         
        Fsini =  1/(info.SamplingInterval*1e-6);
        if isempty(Fs)
            y=tmp;
            parameter.Fs = Fsini;
            Fs =  Fsini;
        else
            y = resample(tmp, Fs,Fsini);
        end
        clear tmp d1
        %tmptrial = reshape(y,info.SegmentDataPoints,  info.DataPoints/info.SegmentDataPoints, numel(indlist));
        %time x trial x ele
        data_all = [];

        vmrk_path = fullfile(pathdata,[filenametot,'.vmrk' ]);
        %         handles.file_vmrk = handles.NIRS.Dt.fir.pp(end).p{idfile}; %used
        %         in save noise
        mrk_type_arr = cellstr('Bad Interval');
        mrks = [];
        ind = [];
        [ind_dur_ch_bad] = read_vmrk_find(vmrk_path,'Bad Interval'); %get BV marking 
        noise = ind_dur_ch2mat(ind_dur_ch_bad, info.DataPoints,info.NumberOfChannels);

        if isempty(Fs)
            noise= double(noise);           
        else
            noise= resample(double(noise), Fs,Fsini);
        end
        
        %RESAMPLE 
        noise = noise(:, indlist);
        %apply noise
        id = find(noise==1);
        y(id) = nan;
        clear id
        nbtrial = info.DataPoints/info.SegmentDataPoints;
        ytrial= reshape(y,size(y,1)/nbtrial,nbtrial,numel(indlist));
        noisetrial = reshape(noise,size(y,1)/nbtrial,nbtrial,numel(indlist));
 %       figure; imagesc(noise)        
 %          [Svv_channel,F,Nseg,PSD] = xspectrum(data,Fs,Fm,deltaf);
        nsample = size(y,1)/nbtrial
        t = (((1/Fs:1/Fs:(1/Fs*nsample)))*1000)+parameter.baseline
        idgood = 1;
        idtrial = [];
        
        %Only keep good trial 
        for itrial = 1:nbtrial
            if sum(noisetrial(:,itrial,:))<10
            data = squeeze(ytrial(:,itrial,:))';
             idgood = idgood +1;
             S(:,:,idgood ) = data';
             idtrial = [idtrial, itrial];
             data_all = [data_all,data];             
            end
        end
        if get(handles.popup_method,'value')==2
            ERP = repmat(nanmean(S,3),1,1,size(S,3));
            S = S - ERP;
            Args.method = 'induce';
        else
            Args.method = 'total';
        end
        %sample  x ele x trial
          %MY ERLCOH FUNCTION      
          for itrial=1:size(S,3)
            itrial          
            St = S(:,:,itrial);       
            TFR(:,:,:,itrial) = morletcwtEd(St,  parameter.layer ,parameter.Fs,  parameter.width);
            %time layer ele trial
          end     
          clear S St ERP noise ind_dur_ch_bad y ytrial noisetrial data_all data
%           figure;%subplot(2,1,2)
%           avgtfr  = mean(mean(abs(TFR),4),3)./(ones(size(TFR,1),1)*mean(mean(mean(abs(TFR(:,:,:)),4),3),1))
%           imagesc(t, parameter.layer, avgtfr');
%           axis xy
          
            %RUN COHERENCE NAN
    nele = numel( indlist);
   sampleremove= round(parameter.ERLCOHcut / 1000 * Fs);
   samplekeep = sampleremove:(numel(t)-sampleremove);
   idbaseline = 1:51;
    TF = TFR(samplekeep,:,:,:); %0 à 600
    size(TFR);
   Args.time =  t(samplekeep);
   Args.layer = parameter.layer;
   Args.ZoneList= nameele;
   Args.Fs = parameter.Fs;
    Args.width = parameter.width;
    Args.ntrial = size(TFR,4);
    WAV = abs(nanmean(TFR.*conj(TFR),4));
    iele = 1
            for iele =1:size(TFR,3)
                alltfX= permute(squeeze(TFR(:,:,iele,:)),[2 1 3]);
                P  = alltfX.*conj(alltfX); % power for wavelets (newtimef 1192)
                Pori  = nanmean(P, 3); 
                 mbase = nanmean(Pori(:,:),2);
                P = bsxfun(@rdivide, P, mbase);          
                P = nanmean(P(:,:, :), 3);
                 P = 10 * log10(P);
                ERSP(:,:,iele)=P;
            end
   Args.idbaseline= idbaseline;
   fileWAV =  [ pathdata, filesep, filenametot,'WAV','.mat'];
   save(fileWAV,'WAV','Args');
   fileERSP =  [ pathdata, filesep, filenametot,'ERSP','.mat'];
   ERSP=permute(ERSP,[2,1,3]);
   save(fileERSP,'ERSP','Args');
    clear ERSP
     % ERLCOH = zeros(numel( parameter.layer),numel(samplekeep),nele,nele);
      id = 1;
      if 1
    for ielex=2:nele
        ielex;
        ieley = 1;
        while ieley < ielex 
            alltfX=permute(TF(:,:,ielex,:),[2,1,4,3]);%layer x time x trial
            alltfY=permute(TF(:,:,ieley,:),[2,1,4,3]);                 
                if 0 %levelmad_ERLCOH %outlier here in the wavelet domain... madnan on erlcoh
                    tmp = abs(alltfX(:,:,:) .* conj(alltfY(:,:,:)));
                    vmad = repmat(squeeze(mad(tmp,1,3)),[1,1,ntrial]);
                    meanmad = repmat(squeeze(nanmean(tmp,3)),[1,1,ntrial]);
                    scoremad = (tmp -  meanmad)./vmad;
                    %figure;hist(scoremad(:))
                    idnan = find(scoremad>madtr);
                    alltfX(idnan) = nan;
                    alltfY(idnan) = nan;                
                end
                %ERLCOH(:,:,ielex,ieley) = nansum(alltfX(:,:,:) .* conj(alltfY(:,:,:)), 3) ./ sqrt(nansum(abs(alltfX(:,:,:)).^2,3) .* nansum(abs(alltfY(:,:,:)).^2,3));
                %ERLCOH(:,:,ielex,ieley) = sum(alltfX(:,:,:) .* conj(alltfY(:,:,:)), 3) ./ sqrt(sum(abs(alltfX(:,:,:)).^2,3) .* sum(abs(alltfY(:,:,:)).^2,3));
                %ERLCOH(:,:,ieley,ielex) = ERLCOH(:,:,ielex,ieley);
                MATCORR(:,:,id) =nansum(alltfX(:,:,:) .* conj(alltfY(:,:,:)), 3) ./ sqrt(nansum(abs(alltfX(:,:,:)).^2,3) .* nansum(abs(alltfY(:,:,:)).^2,3));
                id = id+1
                ieley = ieley + 1;
       
     end
    end
      1
    
    % find cluster in time and frequency from the ERLCOH !
    %X = reshape(MATCORR, numel(samplekeep)*numel( parameter.layer),(nele*nele-nele)/2);
        Args.ZoneList= nameele;
        Args.nbcluster = str2num(get(handles.edit_MaxCluster,'string'));   
        Args.dt = 1/Fs;
        Args.time =   t(samplekeep)
        Args.period = parameter.layer.^2;
        Args.layer = parameter.layer;
        Args.tRatio = 1;
        Args.fwin = 1:numel(parameter.layer);
        Args.coix = 1:numel(samplekeep);
        Args.pathout = pathdata;
        Args.filloutput = filenametot;
        if ~isempty( Args.nbcluster)
            C = clusterMATCORR(MATCORR,Args.nbcluster)
        else
            C = [];
        end
        fileCOH =  [Args.pathout,filesep, Args.filloutput,'COH','.mat']
        save(fileCOH,'MATCORR','C','Args') 
        clear MATCORR TFR TF WAV St S P Pori alltfX alltfY mbase
      end
      1
        end
%     MATCORR
%     size(ERLCOH)
%     [IDX,C,sumd,D] = kmeans(abs(X),4);
%     IDtF =reshape(IDX,numel( parameter.layer),numel(samplekeep))
%     figure
%     imagesc(t(samplekeep), parameter.layer,IDtF)
%     title('Kmeans Connectivity cluster in function of time and frequency')         
%      idcluster = find(IDX==2)
%      avgcluster = mean(X(idcluster,:));
%     avgcluster = reshape(avgcluster, nele,nele)
    


%     figure
%           imagesc( abs( avgcluster))
%           set(gca, 'XTick',1:length(elect_58_343.conv_ASA343),'YTick',1:length(elect_58_343.conv_ASA343),...
%             'XTickLabel',elect_58_343.label,'XTickLabelRotation',90,...
%             'YTickLabel',elect_58_343.label,'YTickLabelRotation',0);
%       end

function  C = clusterMATCORR(MATCORR,ncluster)
        %%%%%%%%%%%%% MATCORR %%%%%
        %Try to find several cluster of connectivity matrix in the time
        %frequency space
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        X = reshape(MATCORR, size(MATCORR,1)*size(MATCORR,2),size(MATCORR,3));
        clear MATCORR
        for icluster= 1:ncluster
            icluster
            [C{icluster}.IDX,C{icluster}.mean,sumd,D] = kmeans(abs(X),icluster);
            distsum(icluster) = sum(sumd);
            C{icluster}.sumd = sumd;
        end
%          figure
%         plot(distsum)
%         xlabel('Nb Cluster')
%         ylabel('Error')


function edit_MaxCluster_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MaxCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MaxCluster as text
%        str2double(get(hObject,'String')) returns contents of edit_MaxCluster as a double


% --- Executes during object creation, after setting all properties.
function edit_MaxCluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MaxCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_frequence_Callback(hObject, eventdata, handles)
% hObject    handle to edit_frequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_frequence as text
%        str2double(get(hObject,'String')) returns contents of edit_frequence as a double


% --- Executes during object creation, after setting all properties.
function edit_frequence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_frequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_elefile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_elefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_elefile as text
%        str2double(get(hObject,'String')) returns contents of edit_elefile as a double


% --- Executes during object creation, after setting all properties.
function edit_elefile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_elefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_width_Callback(hObject, eventdata, handles)
% hObject    handle to edit_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_width as text
%        str2double(get(hObject,'String')) returns contents of edit_width as a double


% --- Executes during object creation, after setting all properties.
function edit_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_openele.
function btn_openele_Callback(hObject, eventdata, handles)
% hObject    handle to btn_openele (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, path]= uigetfile('.ele')
set(handles.edit_elefile,'string',[path,file])



function edit_resample_Callback(hObject, eventdata, handles)
% hObject    handle to edit_resample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_resample as text
%        str2double(get(hObject,'String')) returns contents of edit_resample as a double


% --- Executes during object creation, after setting all properties.
function edit_resample_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_resample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_calculateBCVareta.
function radio_calculateBCVareta_Callback(hObject, eventdata, handles)
% hObject    handle to radio_calculateBCVareta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_calculateBCVareta



function edit_Baseline_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Baseline as text
%        str2double(get(hObject,'String')) returns contents of edit_Baseline as a double


% --- Executes during object creation, after setting all properties.
function edit_Baseline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ERLCOHcut_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ERLCOHcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ERLCOHcut as text
%        str2double(get(hObject,'String')) returns contents of edit_ERLCOHcut as a double


% --- Executes during object creation, after setting all properties.
function edit_ERLCOHcut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ERLCOHcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pathEEG_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pathEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pathEEG as text
%        str2double(get(hObject,'String')) returns contents of edit_pathEEG as a double


% --- Executes during object creation, after setting all properties.
function edit_pathEEG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pathEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_fileEEG.
function listbox_fileEEG_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_fileEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_fileEEG contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_fileEEG


% --- Executes during object creation, after setting all properties.
function listbox_fileEEG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_fileEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_COHCluster.
function btn_COHCluster_Callback(hObject, eventdata, handles)
% hObject    handle to btn_COHCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_COHClusterToMatrices


% --- Executes on selection change in popup_method.
function popup_method_Callback(hObject, eventdata, handles)
% hObject    handle to popup_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_method


% --- Executes during object creation, after setting all properties.
function popup_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_KmeanCluster.
function btn_KmeanCluster_Callback(hObject, eventdata, handles)
% hObject    handle to btn_KmeanCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function btn_KmeanCluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to btn_KmeanCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
GUI_ERLCOH_VIEW_Kmean