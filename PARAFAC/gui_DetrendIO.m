function varargout = gui_DetrendIO(varargin)
% GUI_DETRENDIO M-file for gui_DetrendIO.fig
%      GUI_DETRENDIO, by itself, creates a new GUI_DETRENDIO or raises the existing
%      singleton*.
%
%      H = GUI_DETRENDIO returns the handle to a new GUI_DETRENDIO or the handle to
%      the existing singleton*.
%
%      GUI_DETRENDIO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_DETRENDIO.M with the given input arguments.
%
%      GUI_DETRENDIO('Property','Value',...) creates a new GUI_DETRENDIO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_DetrendIO_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_DetrendIO_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_DetrendIO

% Last Modified by GUIDE v2.5 07-Sep-2018 14:30:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_DetrendIO_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_DetrendIO_OutputFcn, ...
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


% --- Executes just before gui_DetrendIO is made visible.
function gui_DetrendIO_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_DetrendIO (see VARARGIN)

% Choose default command line output for gui_DetrendIO
handles.output = hObject;

%CALCUL INTERACTIF PARAFAC SELECTIONNNER POUR LA FENETRE DE TEMPS DANS LE
%GUI
handles.output = hObject;
tstart = varargin{1};
tstop = varargin{2};
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
d = PMI{currentsub}.data(cf).HRF.AvgC;

indt = tstart(end):tstop(end);
intensnorm = d(indt,:);
 
% %Detrent DATA segment for centrering 
% X = 1:1:size(intensnorm,1);
% Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
% Mb2 =  intensnorm(1,:)'; %offset    
% A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
% spar = intensnorm- A;


listbad = find(PMI{currentsub}.data(cf).MeasListAct==0);
MeasListActplotLst = zeros(size(intensnorm,2),1);
MeasListActplotLst(PMI{1}.plotLst,:)=1;
if min(PMI{1}.plotLst-size(intensnorm,2)/2) >0
    MeasListActplotLst(PMI{1}.plotLst-size(intensnorm,2)/2,:)=1;
end
MeasListActplotLst((size(intensnorm,2)/2+1):end,1)=0; %forcer la fin a zeros

%listgood = find(double(PMI{currentsub}.data(cf).MeasListAct==1).*double(MeasListActplotLst==1));
listgood = PMI{currentsub}.plotLst ; % find(double(PMI{currentsub}.data(cf).MeasListAct==1).*double(MeasListActplotLst==1));
nc830 = numel(PMI{currentsub}.data(cf).MeasListAct)/2;
all830 = zeros(numel(PMI{currentsub}.data(cf).MeasListAct),1);
all830(listgood)=1;
tmp = reshape(all830,nc830 ,2);
t = sum(tmp');
listgood = find([t,t]);
d1 = intensnorm(:,listgood);

t = PMI{currentsub}.data(cf).HRF.tHRF;
PMI{1}.tmpDetrend.checksumd = sum(intensnorm(:)); 
PMI{1}.tmpDetrend.selected =1     
PMI{currentsub}.tmpDetrend.d1 = d1;                    %Data initial
PMI{currentsub}.tmpDetrend.listgood = listgood ; 
PMI{currentsub}.tmpDetrend.indt = [tstart(end),tstop(end)];%Time indice
PMI{currentsub}.tmpDetrend.equationtype = 'linear';
PMI{currentsub}.tmpDetrend.equation = 'NA';
PMI{currentsub}.tmpDetrend.Rsquare = 1;
set(guiHOMER,'UserData',PMI);
btn_LinearTrend_Callback(hObject, eventdata, handles)





% Update handles structure





% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_DetrendIO wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_DetrendIO_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_linearoption.
function popup_linearoption_Callback(hObject, eventdata, handles)
% hObject    handle to popup_linearoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_linearoption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_linearoption
if get(handles.linearoption,'value')==1
    set(handles.popup_displaydetrendtopo,'string',{'B1', 'B2','R2'})
end
    

% --- Executes during object creation, after setting all properties.
function popup_linearoption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_linearoption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_LinearTrend.
function btn_LinearTrend_Callback(hObject, eventdata, handles)
% hObject    handle to btn_LinearTrend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
d1 = PMI{currentsub}.tmpDetrend.d1;
X = PMI{currentsub}.tmpDetrend.indt(1):PMI{currentsub}.tmpDetrend.indt(2);

if  get(handles.popup_linearoption,'value')==1 %linear
    ball = zeros(size(d1,2),2);
    rall = zeros(size(d1,2),1);
    pall = zeros(size(d1,2),1);
    for ich=1:size(d1,2)
        [b,bint,r,rint,stats] = regress(d1(:,ich),[ones(size(X));X]');
         ball(ich,:) = b;
         rall(ich) = stats(1);
         pall(ich) = stats(3);
    end
PMI{currentsub}.tmpDetrend.equationtype = 'linear slope';
PMI{currentsub}.tmpDetrend.equation = ['b1','+' 'b2','*','X' ];
PMI{currentsub}.tmpDetrend.beta = ball;
PMI{currentsub}.tmpDetrend.Rsquare = rall;
PMI{currentsub}.tmpDetrend.p =pall;
PMI{currentsub}.tmpDetrend.display = 2
%Equation type Linear slope
%Display B1,B2, R2

end

set(guiHOMER,'UserData',PMI);
guidata(hObject, handles);

resetview(handles)

function resetview(handles)
%All normal 
newfigure  = 0

cla
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;

    if newfigure==0
    axes(handles.axes_trend);
    cla
    hold on
    else
        figure
        hold on
    end
    if get(handles.popup_linearoption,'value')==1
    t = PMI{currentsub}.data(cf).HRF.tHRF(PMI{currentsub}.tmpDetrend.indt(1):PMI{currentsub}.tmpDetrend.indt(2));
    d1 = PMI{currentsub}.tmpDetrend.d1 ; 

    
    X = PMI{currentsub}.tmpDetrend.indt(1):PMI{currentsub}.tmpDetrend.indt(2);
    ich = 1;
    for ich = 1:size(d1,2)
    Y = PMI{currentsub}.tmpDetrend.beta(ich,1) + PMI{currentsub}.tmpDetrend.beta(ich,2).*X;
        plot(t,d1(:,ich))
        plot(t,Y)    
    end
    
    set(handles.text_equation,'string', PMI{currentsub}.tmpDetrend.equation)
    set(handles.popup_displaydetrendtopo,'value',PMI{currentsub}.tmpDetrend.display)
   % set(handles.text_Rsquare,'string', Rsquare)

    end
% --- Executes on selection change in popup_displaydetrendtopo.
function popup_displaydetrendtopo_Callback(hObject, eventdata, handles)
% hObject    handle to popup_displaydetrendtopo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_displaydetrendtopo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_displaydetrendtopo
value = get(handles.popup_displaydetrendtopo,'value')
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
PMI{currentsub}.tmpDetrend.display = value;

set(guiHOMER,'UserData',PMI);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popup_displaydetrendtopo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_displaydetrendtopo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
