function varargout = gui_COH(varargin)
% GUI_COH MATLAB code for gui_COH.fig
%      GUI_COH, by itself, creates a new GUI_COH or raises the existing
%      singleton*.
%
%      H = GUI_COH returns the handle to a new GUI_COH or the handle to
%      the existing singleton*.
%
%      GUI_COH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_COH.M with the given input arguments.
%
%      GUI_COH('Property','Value',...) creates a new GUI_COH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_COH_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_COH_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_COH

% Last Modified by GUIDE v2.5 22-Oct-2024 10:09:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_COH_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_COH_OutputFcn, ...
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
end

% --- Executes just before gui_COH is made visible.
function gui_COH_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_COH (see VARARGIN)

% Choose default command line output for gui_COH
handles.output = hObject;

if isempty(varargin)
    resetview(handles,0)
else
handles.output = hObject;
tstart = varargin{1};
tstop = varargin{2};
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');
cf = PMI{currentsub}.currentFile;
d = PMI{currentsub}.data(cf).HRF.AvgC;

indt = tstart(end):tstop(end);
if tstart(end)==tstop(end)
    disp('Please verify time start and time stop.')
    msgbox('Please verify time start and time stop.')
end
intensnorm = d(indt,:);
  
%Detrent DATA segment for centrering 
X = 1:1:size(intensnorm,1);
Mb1 =  ((intensnorm(end,:)-intensnorm(1,:))./numel(X))';
Mb2 =  intensnorm(1,:)'; %offset    
A = reshape(X,numel(X),1)*reshape( Mb1,1,numel(Mb1)) +ones(numel(X),1)*reshape( Mb2,1,numel(Mb2));
spar = intensnorm - A;

listbad = find(PMI{currentsub}.data(cf).MeasListAct==0);
MeasListActplotLst = zeros(size(A,2),1);
MeasListActplotLst(PMI{1}.plotLst,:)=1;
if min(PMI{1}.plotLst-size(A,2)/2) >0
    MeasListActplotLst(PMI{1}.plotLst-size(A,2)/2,:)=1;
end
MeasListActplotLst((size(A,2)/2+1):end,1)=0; %forcer la fin a zeros

listgood = find(double(PMI{currentsub}.data(cf).MeasListAct==1).*double(MeasListActplotLst==1));



% nbcomp = str2num(get(handles.edit_nbComponent,'string'));
% t = PMI{currentsub}.data(cf).HRF.tHRF;
% spar = cat(3,spar(:,1:end/2),spar(:,end/2+1:end));
% spartmp = spar(:,listgood,:);
% PMI{1}.tmpPARAFAC.checksumd = sum(intensnorm(:)); 
% PMI{1}.tmpPARAFAC.selected =1;
% PMI{currentsub}.tmpPARAFAC.spar = spartmp;                    %Data initial
% PMI{currentsub}.tmpPARAFAC.listgood = listgood ; 
% PMI{currentsub}.tmpPARAFAC.indt = [tstart(end),tstop(end)];%Time indice
% PMI{currentsub}.tmpPARAFAC.dotted = zeros(nbcomp,1); 

set(guiHOMER,'UserData',PMI);
btn_COH_Callback(hObject, eventdata, handles, tstart,tstop)

% label = [];
% for i=1:nbcomp
%     label = [label,{['F',num2str(i)]}];
% end
% set(handles.popupmenu_Factor,'string',label)

end
% Update handles structure
guidata(hObject, handles);


end


% UIWAIT makes gui_COH wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_COH_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
end

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in btn_COH.
function btn_COH_Callback(hObject, eventdata, handles, tstart,tstop)
% hObject    handle to btn_COH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to btn_Parafac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
currentsub=1;
PMI = get(guiHOMER,'UserData');


d1 = PMI{1}.data(1).HRF.AvgC';
 
 fs = 1/(PMI{1}.data(1).HRF.tHRF(2)-PMI{1}.data(1).HRF.tHRF(1))

t = PMI{1}.data(1).HRF.tHRF;
idwindow = tstart:tstop; 
listHBO = PMI{1}.plotLst;
                St = d1(listHBO,idwindow); %prendre segment pour wavelet
           %     figure;plot(t (idwindow),St')
                if 0
                %from Thierry Beausoleil
                Args.dt = t(2)-t(1);
                Args.Pad = 1;
                Args.Dj = 0.0833;
                Args.J1 = 113;
                Args.Mother = 'Morlet';
                Args.S0 = 2*Args.dt;
                Args.tRatio = 1;
                TFR = zeros(Args.J1+1 ,size(St,2),size(St,1));
                sTFR = zeros(Args.J1+1 ,size(St,2),size(St,1));       %smoothwavelet and power
                %time x layer x ch
                for i=1:size(St,1)
                    i
                    Sig1Seg= St(i,:);
                    [TFR(:,:,i),period,scale,coix]  =wavelet(Sig1Seg,Args.dt ,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);%#ok
                    if i==1
                        sinv=1./(scale');
                    end
                  %  [sTFR(:,:,i)] = smoothwavelet(sinv(:,ones(1,size(TFR,2))).*abs(TFR(:,:,i)).^2,Args.dt,period,Args.Dj,scale);
                end
                Args.period = period;
                Args.scale = scale;
                Args.coix = coix;
                %CHECK COH 
                x =  St(1,:);
                     y =  St(2,:);
                     figure; subplot(2,1,2); plot(x)
                COH=wtc(x,y,Args)

 
                figure
                imagesc(t,Args.period, COH)
                axis i,j
                %MATLAB EXAMPLE
                end
                %figure
                for i = 1:numel(listHBO) 
                    for j = 1:numel(listHBO)
                        if ~(i==j)                   % i = listHBO(1)
                    % j  = listHBO(2)
                  %   axes(handles.axes_COH);
                            figure
                            x = St(i,:)'
                                 y= St(j,:)'
                                %padx=[ flipud(x);x;flipud(x)]
                               %  pady=[ flipud(y);y;flipud(y)]
                            wcoherence(x ,y,fs,'numscales',16)
                            
                            title(['COH ', num2str(listHBO(i)),' ' ,num2str(listHBO(j)), ' TIME:', num2str((PMI{1}.data(1).HRF.tHRF(tstart)))])
                              [WCOH,WCS,F,COI,WTX,WTY]  = wcoherence(St(i,:)',St(j,:)',fs,'numscales',16);
                              figure;imagesc(t(idwindow),F,WCOH);axis xy  
% theta = angle(WCS);
% plotPhaseVectors(AX,theta,tax,pax,tspace,pspace);
% 
% arrowpatch = [-1 0 0 1 0 0 -1; 0.1 0.1 0.5 0 -0.5 -0.1 -0.1]';
% for ii=numel(tgrid):-1:1
%     % Multiply each arrow by the rotation matrix for the given theta
%     rotarrow = arrowpatch*[cos(theta(ii)) sin(theta(ii));...
%         -sin(theta(ii)) cos(theta(ii))];
%     patch(tgrid(ii)+rotarrow(:,1)*dx,pgrid(ii)+rotarrow(:,2)*dy,[0 0 0],...
%         'edgecolor','none' ,'Parent',axhandle);
% end
                        
                        
                      end
              %angle(wcs)
              %plotPhaseVectors(AX,theta,tax,pax,tspace,pspace);
                    end
                end
%figure imagesc(wcoh




                %[wcoh,~,P,coi] = wcoherence(St(1,:)',St(2,:)',seconds(fs),'numscales',16);
             
            
                if 0
                            id = 1;
                for i=1:size(TFR,3)
                    j = 1;
                    while j<i
                        Gxy(:,1,id) = TFR(:,1,i).*conj(TFR(:,1,j));
                        id = id + 1;
                        j = j+1;
                    end
                end
                 i = 15; 
                    j = 20;
                     samplemat = 1:size(Gxy,1);
                    in1 =  TFR(:,:,i);
                    in1 = in1(:);
                 in2 =  TFR(: ,:,j);
                        in2 = in2(:);
                        COVC1C2  = nansum(in1.*conj(in2),1);
                        COVC1 =nansum(in1.*conj(in1),1);
                        COVC2 =nansum(in2.*conj(in2),1);
                        temp = abs(COVC1C2).^2 ./ (COVC1.*COVC2);
                figure; imagesc (temp )
                
                end


set(guiHOMER,'UserData',PMI);
guidata(hObject, handles);

resetview(handles,0)
end

function resetview(handles,newfigure)
%All normal 
% 
 axes(handles.axes_COH);
% cla
% guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
% currentsub=1;
% PMI = get(guiHOMER,'UserData');
% 
% 
%     FacTime = PMI{currentsub}.tmpPARAFAC.Factors{1};
%     FacSpatial = PMI{currentsub}.tmpPARAFAC.Factors{2};
%     FacHemo = PMI{currentsub}.tmpPARAFAC.Factors{3};
%     Nc = str2num(get(handles.edit_nbComponent,'string'));
%     colorlist = [0 0 1,
%         1 0 0,
%         0 1 0,
%         0 0 0,
%         1 0.5 0];
%     cf = PMI{currentsub}.currentFile;
%     if newfigure==0
%         axes(handles.axes_timecourse);
%         cla
%         hold on
%     else
%         figure
%         subplot(2,2,1)
%         hold on
%     end
%     t = PMI{currentsub}.data(cf).HRF.tHRF(PMI{currentsub}.tmpPARAFAC.indt(1):PMI{currentsub}.tmpPARAFAC.indt(2));
% 
%     for i=1:Nc
%          if PMI{1}.tmpPARAFAC.selected==i
% 
%             hline = plot(t, FacTime(:,i),'color',colorlist(i,:),'linewidth',4 );
%          else
%             hline = plot(t, FacTime(:,i),'color',colorlist(i,:) );
%          end
%          try
%          if PMI{currentsub}.tmpPARAFAC.dotted(i)==1
%              set(hline,'linestyle','--');
%          else
%               set(hline,'linestyle','-');
%          end
%          catch
%              1
%          end
% 
%     end
% 
%     if newfigure==0
%         axes(handles.axes_Topo);
%         cla
%         hold on
%     else
%         subplot(2,2,2)
%         hold on
%     end
%     for i=1:Nc
%         if PMI{1}.tmpPARAFAC.selected==i
%             plot(FacSpatial(:,i),'color',colorlist(i,:),'linewidth',4 )
%         else
%             plot(FacSpatial(:,i),'color',colorlist(i,:) )
%         end
%     end
%     if newfigure==0
%         axes(handles.axes_hemo);
%         cla
%         hold on
%      else
%         subplot(2,2,3)
%         hold on
%     end
%     for i=1:Nc
%         if PMI{1}.tmpPARAFAC.selected==i
%             plot(FacHemo(:,i),'color',colorlist(i,:),'linewidth',4 )
%         else
%             plot(FacHemo(:,i),'color',colorlist(i,:) )        
%         end
%     end
end
    %WAVELET  1D Wavelet transform with optional singificance testing
%
%   [WAVE,PERIOD,SCALE,COI] = wavelet(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
%
%   Computes the wavelet transform of the vector Y (length N),
%   with sampling rate DT.
%
%   By default, the Morlet wavelet (k0=6) is used.
%   The wavelet basis is normalized to have total energy=1 at all scales.
%
%
% INPUTS:
%
%    Y = the time series of length N.
%    DT = amount of time between each Y value, i.e. the sampling time.
%
% OUTPUTS:
%
%    WAVE is the WAVELET transform of Y. This is a complex array
%    of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
%    ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
%    The WAVELET power spectrum is ABS(WAVE)^2.
%    Its units are sigma^2 (the time series variance).
%
%
% OPTIONAL INPUTS:
%
% *** Note *** setting any of the following to -1 will cause the default
%               value to be used.
%
%    PAD = if set to 1 (default is 0), pad time series with enough zeroes to get
%         N up to the next higher power of 2. This prevents wraparound
%         from the end of the time series to the beginning, and also
%         speeds up the FFT's used to do the wavelet transform.
%         This will not eliminate all edge effects (see COI below).
%
%    DJ = the spacing between discrete scales. Default is 0.25.
%         A smaller # will give better scale resolution, but be slower to plot.
%
%    S0 = the smallest scale of the wavelet.  Default is 2*DT.
%
%    J1 = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
%        to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
%
%    MOTHER = the mother wavelet function.
%             The choices are 'MORLET', 'PAUL', or 'DOG'
%
%    PARAM = the mother wavelet parameter.
%            For 'MORLET' this is k0 (wavenumber), default is 6.
%            For 'PAUL' this is m (order), default is 4.
%            For 'DOG' this is m (m-th derivative), default is 2.
%
%
% OPTIONAL OUTPUTS:
%
%    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
%           to the SCALEs.
%
%    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J1
%            where J1+1 is the total # of scales.
%
%    COI = if specified, then return the Cone-of-Influence, which is a vector
%        of N points that contains the maximum period of useful information
%        at that particular time.
%        Periods greater than this are subject to edge effects.
%        This can be used to plot COI lines on a contour plot by doing:
%
%              contour(time,log(period),log(power))
%              plot(time,log(coi),'k')
%
%----------------------------------------------------------------------------
%   Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
%   University of Colorado, Program in Atmospheric and Oceanic Sciences.
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
%
% Notice: Please acknowledge the use of this program in any publications:
%   ``Wavelet software was provided by C. Torrence and G. Compo,
%     and is available at URL: http://paos.colorado.edu/research/wavelets/''.
%
% Notice: Please acknowledge the use of the above software in any publications:
%    ``Wavelet software was provided by C. Torrence and G. Compo,
%      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
%
% Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
%            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
%
% Please send a copy of such publications to either C. Torrence or G. Compo:
%  Dr. Christopher Torrence               Dr. Gilbert P. Compo
%  Advanced Study Program                 NOAA/CIRES Climate Diagnostics Center
%  National Center for Atmos. Research    Campus Box 216
%  P.O. Box 3000                          University of Colorado at Boulder
%  Boulder CO 80307--3000, USA.           Boulder CO 80309-0216, USA.
%  E-mail: torrence@ucar.edu              E-mail: gpc@cdc.noaa.gov
%----------------------------------------------------------------------------
function [wave,period,scale,coi] = ...
    wavelet(Y,dt,pad,dj,s0,J1,mother,param)

if (nargin < 8), param = -1; end
if (nargin < 7), mother = -1; end
if (nargin < 6), J1 = -1; end
if (nargin < 5), s0 = -1; end
if (nargin < 4), dj = -1; end
if (nargin < 3), pad = 0; end
if (nargin < 2)
    error('Must input a vector Y and sampling time DT')
end

n1 = length(Y);

if (s0 == -1), s0=2*dt; end
if (dj == -1), dj = 1./4.; end
if (J1 == -1), J1=fix((log(n1*dt/s0)/log(2))/dj); end
if (mother == -1), mother = 'MORLET'; end

%....construct time series to analyze, pad if necessary
x(1:n1) = Y - mean(Y);
if (pad == 1)
    base2 = fix(log(n1)/log(2) + 0.4999);   % power of 2 nearest to N
    x = [x,zeros(1,2^(base2+1)-n1)];
end
n = length(x);

%....construct wavenumber array used in transform [Eqn(5)]
k = [1:fix(n/2)];
k = k.*((2.*pi)/(n*dt));
k = [0., k, -k(fix((n-1)/2):-1:1)];

%....compute FFT of the (padded) time series
f = fft(x);    % [Eqn(3)]

%....construct SCALE array & empty PERIOD & WAVE arrays
scale = s0*2.^((0:J1)*dj);
period = scale;
wave = zeros(J1+1,n);  % define the wavelet array
wave = wave + i*wave;  % make it complex

% loop through all scales and compute transform
for a1 = 1:J1+1
    [daughter,fourier_factor,coi,dofmin]=wave_bases(mother,k,scale(a1),param);
    wave(a1,:) = ifft(f.*daughter);  % wavelet transform[Eqn(4)]
end

period = fourier_factor*scale;
coi = coi*dt*[1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];  % COI [Sec.3g]
wave = wave(:,1:n1);  % get rid of padding before returning

return
end
%WAVE_BASES  1D Wavelet functions Morlet, Paul, or DOG
%
%  [DAUGHTER,FOURIER_FACTOR,COI,DOFMIN] = ...
%      wave_bases(MOTHER,K,SCALE,PARAM);
%
%   Computes the wavelet function as a function of Fourier frequency,
%   used for the wavelet transform in Fourier space.
%   (This program is called automatically by WAVELET)
%
% INPUTS:
%
%    MOTHER = a string, equal to 'MORLET' or 'PAUL' or 'DOG'
%    K = a vector, the Fourier frequencies at which to calculate the wavelet
%    SCALE = a number, the wavelet scale
%    PARAM = the nondimensional parameter for the wavelet function
%
% OUTPUTS:
%
%    DAUGHTER = a vector, the wavelet function
%    FOURIER_FACTOR = the ratio of Fourier period to scale
%    COI = a number, the cone-of-influence size at the scale
%    DOFMIN = a number, degrees of freedom for each point in the wavelet power
%             (either 2 for Morlet and Paul, or 1 for the DOG)
%
%----------------------------------------------------------------------------
%   Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
%   University of Colorado, Program in Atmospheric and Oceanic Sciences.
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
%----------------------------------------------------------------------------
function [daughter,fourier_factor,coi,dofmin] = ...
    wave_bases(mother,k,scale,param)

mother = upper(mother);
n = length(k);

if (strcmp(mother,'MORLET'))  %-----------------------------------  Morlet
    if (param == -1), param = 6.; end
    k0 = param;
    expnt = -(scale.*k - k0).^2/2.*(k > 0.);
    norm = sqrt(scale*k(2))*(pi^(-0.25))*sqrt(n);    % total energy=N   [Eqn(7)]
    daughter = norm*exp(expnt);
    daughter = daughter.*(k > 0.);     % Heaviside step function
    fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); % Scale-->Fourier [Sec.3h]
    coi = fourier_factor/sqrt(2);                  % Cone-of-influence [Sec.3g]
    dofmin = 2;                                    % Degrees of freedom
elseif (strcmp(mother,'PAUL'))  %--------------------------------  Paul
    if (param == -1), param = 4.; end
    m = param;
    expnt = -(scale.*k).*(k > 0.);
    norm = sqrt(scale*k(2))*(2^m/sqrt(m*prod(2:(2*m-1))))*sqrt(n);
    daughter = norm*((scale.*k).^m).*exp(expnt);
    daughter = daughter.*(k > 0.);     % Heaviside step function
    fourier_factor = 4*pi/(2*m+1);
    coi = fourier_factor*sqrt(2);
    dofmin = 2;
elseif (strcmp(mother,'DOG'))  %--------------------------------  DOG
    if (param == -1), param = 2.; end
    m = param;
    expnt = -(scale.*k).^2 ./ 2.0;
    norm = sqrt(scale*k(2)/gamma(m+0.5))*sqrt(n);
    daughter = -norm*(i^m)*((scale.*k).^m).*exp(expnt);
    fourier_factor = 2*pi*sqrt(2./(2*m+1));
    coi = fourier_factor/sqrt(2);
    dofmin = 1;
else
    error('Mother must be one of MORLET,PAUL,DOG')
end

return
end

function [coef, pval] = corr(x, varargin)
%CORR Linear or rank correlation.
%   RHO = CORR(X) returns a P-by-P matrix containing the pairwise linear
%   correlation coefficient between each pair of columns in the N-by-P
%   matrix X.
%
%   RHO = CORR(X,Y,...) returns a P1-by-P2 matrix containing the pairwise
%   correlation coefficient between each pair of columns in the N-by-P1 and
%   N-by-P2 matrices X and Y.
%
%   [RHO,PVAL] = CORR(...) also returns PVAL, a matrix of p-values for
%   testing the hypothesis of no correlation against the alternative that
%   there is a non-zero correlation.  Each element of PVAL is the p-value
%   for the corresponding element of RHO.  If PVAL(i,j) is small, say less
%   than 0.05, then the correlation RHO(i,j) is significantly different
%   from zero.
%
%   [...] = CORR(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%   parameters and their values.  Valid parameters are the following:
%
%        Parameter  Value
%         'type'    'Pearson' (the default) to compute Pearson's linear
%                   correlation coefficient, 'Kendall' to compute Kendall's
%                   tau, or 'Spearman' to compute Spearman's rho.
%         'rows'    'all' (default) to use all rows regardless of missing
%                   values (NaNs), 'complete' to use only rows with no
%                   missing values, or 'pairwise' to compute RHO(i,j) using
%                   rows with no missing values in column i or j.
%         'tail'    The alternative hypothesis against which to compute
%                   p-values for testing the hypothesis of no correlation.
%                   Choices are:
%                      TAIL         Alternative Hypothesis
%                   ---------------------------------------------------
%                     'both'     correlation is not zero (the default)
%                     'right'    correlation is greater than zero
%                     'left'     correlation is less than zero
%
%   The 'pairwise' option for the 'rows' parameter can produce RHO that is
%   not positive semi-definite.  The 'complete' option always produces a
%   positive semi-definite RHO, but when data are missing, the estimates
%   may be based on fewer observations.
%
%   CORR computes p-values for Pearson's correlation using a Student's t
%   distribution for a transformation of the correlation.  This is exact
%   when X and Y are normal.  CORR computes p-values for Kendall's tau and
%   Spearman's rho using either the exact permutation distributions (for
%   small sample sizes), or large-sample approximations.
%
%   CORR computes p-values for the two-tailed test by doubling the more
%   significant of the two one-tailed p-values.
%
%   See also CORRCOEF, PARTIALCORR, TIEDRANK.

%   References:
%      [1] Gibbons, J.D. (1985) Nonparametric Statistical Inference,
%          2nd ed., M. Dekker.
%      [2] Hollander, M. and D.A. Wolfe (1973) Nonparametric Statistical
%          Methods, Wiley.
%      [3] Kendall, M.G. (1970) Rank Correlation Methods, Griffin.
%      [4] Best, D.J. and D.E. Roberts (1975) "Algorithm AS 89: The Upper
%          Tail Probabilities of Spearman's rho", Applied Statistics,
%          24:377-379.

%   Copyright 1984-2018 The MathWorks, Inc.


%   Spearman's rho is equivalent to the linear (Pearson's) correlation
%   between ranks.  Kendall's tau is equivalent to the linear correlation
%   between the "concordances" sign(x(i)-x(j))*sign(y(i)-y(j)), i<j, with
%   an adjustment for ties.  This is often referred to as tau-b.
%
%   Kendall's tau-b is identical to the standard tau (or tau-a) when there
%   are no ties.  However, tau-b includes an adjustment for ties in the
%   normalizing constant.

%   Spearman's rho and Kendall's tau are discrete-valued statistics, and
%   their distributions have positive probability at 1 and -1.  For small
%   sample sizes, CORR uses the exact permutation distributions, and thus,
%   the on-diagonal p-values from CORR(X,X) in those cases.

%   When there are ties in the data, the null distribution of Spearman's
%   rho and Kendall's tau may not be symmetric.  Computing a two-tailed
%   p-value in such cases is not well-defined.  CORR computes p-values for
%   the two-tailed test by doubling the smaller of the one-tailed p-values.

if nargin > 1
    [varargin{:}] = varargin{:};
end

if nargin < 1 || isempty(x)
    error(message('stats:corr:TooFewInputs'));
end
[n,p1] = size(x);

% Only x given, compute pairwise rank correlations
if (nargin < 2) || ischar(varargin{1})
    corrXX = true;
    y = x;
    p2 = p1;

% Both x and y given, compute the pairwise rank cross correlations
else
    y = varargin{1};
    varargin = varargin(2:end);
    if size(y,1) ~= n
        error(message('stats:corr:InputSizeMismatch'));
    end
    corrXX = false;
    p2 = size(y,2);
end

pnames = {'type'  'rows' 'tail'};
dflts  = {'p'     'a'    'both'};
try
[type,rows,tail] = internal.stats.parseArgs(pnames,dflts,varargin{:});
catch
    type = 'p';
    rows = 'a';
    tail = 'both';
end

% Validate the rows parameter.
rowsChoices = {'all' 'complete' 'pairwise'};
if ischar(rows)
    i = find(strncmpi(rows,rowsChoices,length(rows)));
    if isscalar(i)
        rows = rowsChoices{i}(1);
    end
end
switch rows
case 'a' % 'all'
    % Missing values are left in, so all cols have the same length.
case 'c' % 'complete'
    % Complete rows with missing values are removed, so all cols have the
    % same length.
    ok = ~any(isnan(x),2);
    if ~corrXX
        ok = ok & ~any(isnan(y),2);
    end
    n = sum(ok);
    if n < size(x,1)
        x = x(ok,:);
        if ~corrXX
            y = y(ok,:);
        end
    end
case 'p' % 'pairwise'
    % Missing values are removed pairwise, so each column pair may have a
    % different length, we'll have to see.
otherwise
    error(message('stats:corr:UnknownRows'));
end

% Can't do anything without at least two observations.
if n < 2
    outClass = superiorfloat(x,y);
    coef = NaN(p1,p2,outClass);
    pval = NaN(p1,p2,outClass);
    return
end

% Validate the type parameter, and set a handle to the correlation function.
typeChoices = {'pearson' 'kendall' 'spearman'};
if ischar(type)
    i = find(strncmpi(type,typeChoices,length(type)));
    if isscalar(i)
        type = typeChoices{i}(1);
    end
end
switch type
case 'p' % Pearson linear correlation
    corrFun = @corrPearson;
case 'k' % Kendall rank correlation
    corrFun = @corrKendall;
case 's' % Spearman rank correlation
    corrFun = @corrSpearman;
otherwise
    error(message('stats:corr:UnknownType'));
end

% Validate the tail parameter.
tailChoices = {'left','both','right'};
if ischar(tail) && (size(tail,1)==1)
    i = find(strncmpi(tail,tailChoices,length(tail)));
    if isempty(i)
        i = find(strncmpi(tail,{'lt','ne','gt'},length(tail)));
    end
    if isscalar(i)
        tail = tailChoices{i}(1);
    elseif isempty(i)
        error(message('stats:corr:UnknownTail'));
    end
else
    error(message('stats:corr:UnknownTail'));
end

% Some things we just can't do for complex data.
complexdata = ~(isreal(x) || all(imag(x(:))==0));
if ~complexdata && ~corrXX
    complexdata = ~(isreal(y) || all(imag(y(:))==0));
end
if complexdata
    if nargout > 1
        error(message('stats:corr:PValComplexData'));
    elseif type ~= 'p'
        error(message('stats:corr:RankCorComplexData'));
    end
end

if corrXX
    if nargout<2
        coef = corrFun(rows,tail,x);
    else
        [coef,pval] = corrFun(rows,tail,x);
    end
else
    if nargout<2
        coef = corrFun(rows,tail,x,y);
    else
        [coef,pval] = corrFun(rows,tail,x,y);
    end
end
end

%--------------------------------------------------------------------------

function [coef,pval] = corrPearson(rows,tail,x,y)
%CORRPEARSON Computation of Pearson correlation matrix
[n,p1] = size(x);
corrXX = (nargin < 4);
if corrXX
    p2 = p1;
    y = x;
else
    p2 = size(y,2);
end

% Do column-pairwise NaN removal only if there are missing values.
if rows == 'p'
    isnanX = isnan(x);
    pairwise = any(isnanX(:));
    if corrXX
        isnanY = isnanX;
    else
        isnanY = isnan(y);
        pairwise = pairwise || any(isnanY(:));
    end
else
    pairwise = false;
end

% Computation of linear correlation with column-pairwise NaN removal
if pairwise
    % Preallocate the outputs.
    outClass = superiorfloat(x,y);
    coef = zeros(p1,p2,outClass);
    nij = zeros(p1,p2,outClass);
    for i = 1:p1
        xi = x(:,i);
        isnanXi = isnanX(:,i);

        % Loop over all columns of y for corr(x,y), or over the lower triangle
        % without the diagonal for corr(x).
        j1 = p2*(1-corrXX) + (i-1)*corrXX; % 1:p2, or 1:(i-1)
        for j = 1:j1
            yj = y(:,j);
            ok = ~(isnanXi | isnanY(:,j));
            nij(i,j) = sum(ok);
            if nij(i,j) < n
                x0 = xi(ok); x0 = x0 - sum(x0)/nij(i,j);
                y0 = yj(ok); y0 = y0 - sum(y0)/nij(i,j);
            else
                x0 = xi - sum(xi)/n;
                y0 = yj - sum(yj)/n;
            end
            coef(i,j) = (x0'*y0) ./ (norm(x0) .* norm(y0));
        end
        % Compute diagonal elements for corr(x) to match the non-pairwise
        % case.  This is either 1, if the data in xi are not constant, or NaN
        % if they are.
        if corrXX
            ok = ~isnanXi;
            nij(i,i) = sum(ok);
            if nij(i,i) < n
                x0 = xi(ok); x0 = x0 - mean(x0);
            else
                x0 = xi - mean(xi);
            end
            coef(i,i) = (x0'*x0)./norm(x0)^2;
        end
    end

    % If this is autocorrelation, reflect the lower triangle into the upper.
    if corrXX
        coef = tril(coef) + tril(coef,-1)'; % leave the diagonal alone
    end

    n = nij;

% Computation of linear correlation, all elements at once.  NaNs not removed,
% or have already been removed in complete rows.
else
    x = bsxfun(@minus,x,sum(x,1)/n);  % Remove mean
    if corrXX
        coef = x' * x; % 1/(n-1) doesn't matter, renormalizing anyway
        d = sqrt(diag(coef)); % sqrt first to avoid under/overflow
        coef = bsxfun(@rdivide,coef,d); coef = bsxfun(@rdivide,coef,d'); % coef = coef ./ d*d';
    else
        y = bsxfun(@minus,y,sum(y,1)/n);  % Remove mean
        coef = x' * y; % 1/(n-1) doesn't matter, renormalizing anyway      
        dx = norm(x);    %dx = vecnorm(x,2,1);
        dy = norm(y);      %vecnorm(y,2,1);
        coef = bsxfun(@rdivide,coef,dx'); coef = bsxfun(@rdivide,coef,dy); % coef = coef ./ dx'*dy;
    end
end

% Limit off-diag elements to [-1,1], and put exact ones on the diagonal for
% autocorrelation.
t = find(abs(coef) > 1); coef(t) = coef(t)./abs(coef(t)); % preserves NaNs
if corrXX
    coef(1:p1+1:end) = sign(diag(coef)); % preserves NaNs
end

if nargout > 1
    if corrXX
        pval = zeros(p1);
        ltri = (tril(ones(size(pval)),-1) > 0); % lower triangle only, no diagonal
        if pairwise
            pval(ltri) = pvalPearson(tail, coef(ltri), n(ltri));
        else
            pval(ltri) = pvalPearson(tail, coef(ltri), n);
        end
        % Reflect the p-values from the lower triangle into the upper.
        pval = pval + pval';
        % The p-values along the diagonal are always exactly one (unless
        % they're NaN), because Pr{coef(i,i) as/more extreme than 1} == 1,
        % regardless of which tail(s) we're testing against.
        pval(1:p1+1:end) = sign(diag(coef)); % preserves NaNs on diag
    else
        pval = pvalPearson(tail, coef, n);
    end
end
end

%--------------------------------------------------------------------------

function [coef,pval] = corrKendall(rows,tail,x,y)
%CORRKENDALL Computation of Kendall correlation matrix.
[n,p1] = size(x);
corrXX = (nargin < 4);
if corrXX
    p2 = p1;
    y = x;
else
    p2 = size(y,2);
end

% Do column-pairwise NaN removal only if there are missing values.
if rows == 'p'
    isnanX = isnan(x);
    pairwise = any(isnanX(:));
    if corrXX
        isnanY = isnanX;
    else
        isnanY = isnan(y);
        pairwise = pairwise || any(isnanY(:));
    end
else
    pairwise = false;
end

% Compute the ranks once when not doing column-pairwise NaN removal.
if ~pairwise
    if rows == 'a'
        % For 'all', it's quicker to recognize that a column has a NaN in it
        % than to do the computation explicitly and end up with a NaN anyway.
        nancolX = any(isnan(x),1)';
        if corrXX
            nancol = bsxfun(@or,nancolX,nancolX');
        else
            nancolY = any(isnan(y),1);
            nancol = bsxfun(@or,nancolX,nancolY);
        end

    else % rows == 'c' | rows  == 'p'
        % For 'complete' or 'pairwise', there are no NaNs at this point
        nancol = false(p1,p2);
    end

    [xrank, xadj] = tiedrank(x,1);
    if corrXX
        yrank = xrank; yadj = xadj;
    else
        [yrank, yadj] = tiedrank(y,1);
    end
end

% Preallocate the outputs.
outClass = superiorfloat(x,y);
coef = zeros(p1,p2,outClass);
if nargout > 1
    pval = zeros(p1,p2,outClass);
    Kstat = zeros(p1,p2,outClass); % save the obs. stat. for p-value computation
    if corrXX
         needPVal = tril(true(p1),-1); % lower triangle only, no diagonal
    else
         needPVal = true(p1,p2);
    end
end

for i = 1:p1
    if ~pairwise
        xranki = xrank(:,i); xadji = xadj(:,i);
    end

    % Loop over all columns of y for corr(x,y), or over the lower triangle and
    % the diagonal for corr(x).
    j1 = p2*(1-corrXX) + i*corrXX; % 1:p2, or 1:i
    for j = 1:j1
        if pairwise % column-pairwise NaN removal
            ok = ~(isnanX(:,i) | isnanY(:,j));
            nij = sum(ok);
            anyRowsRemoved = (nij < n);
            if anyRowsRemoved
                [xranki, xadji] = tiedrank(x(ok,i),1);
                [yrankj, yadjj] = tiedrank(y(ok,j),1);
            else
                [xranki, xadji] = tiedrank(x(:,i),1);
                [yrankj, yadjj] = tiedrank(y(:,j),1);
            end
        else % no NaN removal, or NaNs already removed in complete rows
            % Quicker to check and bail out rather than compute the NaN below
            if nancol(i,j)
                coef(i,j) = NaN;
                if nargout > 1
                    pval(i,j) = NaN;
                    needPVal(i,j) = false;
                end
                continue
            end
            nij = n;
            yrankj = yrank(:,j); yadjj = yadj(:,j);
            anyRowsRemoved = false;
        end
        n2const = nij*(nij-1) / 2;
        ties = ((xadji(1)>0) || (yadjj(1)>0));

        K = 0;
        for k = 1:nij-1
            K = K + sum(sign(xranki(k)-xranki(k+1:nij)).*sign(yrankj(k)-yrankj(k+1:nij)));
        end
        coef(i,j) = K ./ sqrt((n2const - xadji(1)).*(n2const - yadjj(1)));

        % Clean up the diagonal for autocorrelation
        if (i == j) && corrXX
            if ~ties
                % Put an exact one on the diagonal only when there's no ties.
                % Kendall's tau may not be exactly one along the diagonal when
                % there are ties.
                coef(i,i) = sign(coef(i,i)); % preserves NaN
            end

            % Compute on-diag p-values for autocorrelation later

        elseif nargout > 1
            % If there are ties, or if there has been pairwise removal of
            % missing data, compute the p-value separately here.
            if ties
                % The tied case is sufficiently slower that we don't
                % want to do it if not necessary.
                %
                % When either of the data vectors is constant, stdK should be
                % zero, get that exactly.
                if (xadji(1) == n2const) || (yadjj(1) == n2const)
                    stdK = 0;
                else
                    stdK = sqrt(n2const*(2*nij+5)./9 ...
                            + xadji(1)*yadjj(1)./n2const ...
                            + xadji(2)*yadjj(2)./(18*n2const*(nij-2)) ...
                            - (xadji(3) + yadjj(3))./18);
                end
                pval(i,j) = pvalKendall(tail, K, stdK, xranki, yrankj);
                needPVal(i,j) = false; % this one's done
            elseif anyRowsRemoved
                stdK = sqrt(n2const*(2*nij+5)./9);
                pval(i,j) = pvalKendall(tail, K, stdK, nij);
                needPVal(i,j) = false; % this one's done
            else
                Kstat(i,j) = K;
            end
        end
    end
end
% Limit off-diag correlations to [-1,1].
t = find(abs(coef) > 1); coef(t) = coef(t)./abs(coef(t)); % preserves NaNs

% Calculate the remaining p-values, except not the on-diag elements for
% autocorrelation.  All cases with no ties and no removed missing values can
% be computed based on a single null distribution.
if nargout > 1 && any(needPVal(:))
    n2const = n*(n-1) ./ 2;
    stdK = sqrt(n2const * (2*n+5) ./ 9);
    pval(needPVal) = pvalKendall(tail,Kstat(needPVal),stdK,n);
end

% If this is autocorrelation, reflect the lower triangle into the upper.
if corrXX
    coef = tril(coef) + tril(coef,-1)'; % leave the diagonal alone
    if nargout > 1
        pval = pval + pval';
        % The p-values along the diagonal are always exactly one (unless
        % they're NaN) conditional on the pattern of ties, because
        % Pr{coef(i,i) as/more extreme than observed value} == 1, regardless
        % of which tail(s) we're testing against.
        pval(1:p1+1:end) = sign(diag(coef)); % preserves NaNs on diag
    end
end
end

%--------------------------------------------------------------------------

function [coef,pval] = corrSpearman(rows,tail,x,y)
%CORRSPEARMAN Computation of Spearman correlation matrix.
[n,p1] = size(x);
corrXX = (nargin < 4);
if corrXX
    p2 = p1;
    y = x;
else
    p2 = size(y,2);
end
outClass = superiorfloat(x,y);

% Do column-pairwise NaN removal only if there are missing values.
if rows == 'p'
    isnanX = isnan(x);
    pairwise = any(isnanX(:));
    if corrXX
        isnanY = isnanX;
    else
        isnanY = isnan(y);
        pairwise = pairwise || any(isnanY(:));
    end
else
    pairwise = false;
end

% Computation of rank correlation, with column-pairwise NaN removal.  Can't
% hand off to corrPearson as for 'all' or 'complete', because we need to
% compute the ranks separately for each pair of columns.
if pairwise
    % Preallocate the outputs.
    coef = zeros(p1,p2,outClass);
    if nargout > 1
        pval = zeros(p1,p2,outClass);
        Dstat = zeros(p1,p2,outClass); % save the obs. stat. for p-value computation
        if corrXX
             needPVal = tril(true(p1),-1); % lower triangle only, no diagonal
        else
             needPVal = true(p1,p2);
        end
    end

    for i = 1:p1
        % Loop over all columns of y for corr(x,y), or over the lower triangle
        % and the diagonal for corr(x).
        j1 = p2*(1-corrXX) + i*corrXX; % 1:p2, or 1:i
        for j = 1:j1
            ok = ~(isnanX(:,i) | isnanY(:,j));
            nij = sum(ok);
            anyRowsRemoved = (nij < n);
            if anyRowsRemoved
                [xranki, xadj] = tiedrank(x(ok,i),0);
                [yrankj, yadj] = tiedrank(y(ok,j),0);
            else
                [xranki, xadj] = tiedrank(x(:,i),0);
                [yrankj, yadj] = tiedrank(y(:,j),0);
            end
            n3const = (nij+1)*nij*(nij-1) ./ 3;
            ties = ((xadj>0) || (yadj>0));

            D = sum((xranki - yrankj).^2);
            meanD = (n3const - (xadj+yadj)./3) ./ 2;

            % When either of the data vectors is constant, stdD should be
            % zero, get that exactly.
            n3const2 = (nij+1)*nij*(nij-1)/2;
            if (xadj == n3const2) || (yadj == n3const2)
                stdD = 0;
            else
                stdD = sqrt((n3const./2 - xadj./3).*(n3const./2 - yadj./3)./(nij-1));
            end

            coef(i,j) = (meanD - D) ./ (sqrt(nij-1)*stdD);
            if (i == j) && corrXX
                % Put an exact one on the diagonal for autocorrelation.
                coef(i,i) = sign(coef(i,i)); % preserves NaNs

                % Compute on-diag p-values for autocorrelation later

            elseif nargout > 1
                % If there are ties, or if there has been pairwise removal
                % of missing data, we'll compute the p-value separately here.
                if (anyRowsRemoved || ties)
                    pval(i,j) = pvalSpearman(tail, D, meanD, stdD, xranki, yrankj);
                    needPVal(i,j) = false; % this one's done
                else
                    Dstat(i,j) = D;
                end
            end
        end
    end
    % Limit off-diag correlations to [-1,1].
    t = find(abs(coef) > 1); coef(t) = coef(t)./abs(coef(t)); % preserves NaNs

    % Calculate the remaining p-values, except not the on-diag elements for
    % autocorrelation.  All cases with no ties and no removed missing values
    % can be computed based on a single null distribution.
    if nargout > 1 && any(needPVal(:))
        n3const = (n+1)*n*(n-1) ./ 3;
        meanD = n3const./2;
        stdD = n3const ./ (2*sqrt(n-1));
        pval(needPVal) = pvalSpearman(tail,Dstat(needPVal),meanD,stdD,n);
    end

    % If this is autocorrelation, reflect the lower triangle into the upper.
    if corrXX
        coef = tril(coef) + tril(coef,-1)'; % leave the diagonal alone
        if nargout > 1
            pval = pval + pval';
            % The p-values along the diagonal are always exactly one (unless
            % they're NaN), because Pr{coef(i,i) as/more extreme than 1} == 1,
            % regardless of which tail(s) we're testing against.
            pval(1:p1+1:end) = sign(diag(coef)); % preserves NaNs on diag
        end
    end

% Vectorized computation of rank correlation.  No NaN removal, or NaNs already
% removed in complete rows
else
    [xrank, xadj] = tiedrank(x,0);
    ties = any(xadj>0);
    if corrXX
        yrank = xrank; yadj = xadj;
    else
        [yrank, yadj] = tiedrank(y,0);
        ties = ties || any(yadj>0);
    end
    n3const = (n+1)*n*(n-1) ./ 3;

    % Compute all elements at once, fastest when there are no ties, or
    % if we don't need p-values
    if nargout == 1 || ~ties
        if corrXX
            coef = corrPearson(rows,tail,xrank);
        else
            coef = corrPearson(rows,tail,xrank,yrank);
        end
        if nargout > 1 % p-values requested, but no ties
            meanD = n3const./2;
            stdD = n3const ./ (2*sqrt(n-1));
            D = meanD - round(coef*(sqrt(n-1)*stdD));
            if corrXX
                % Compute p-values for the lower triangle and reflect into the upper.
                pval = zeros(p1,outClass);
                ltri = (tril(ones(size(pval)),-1) > 0); % lower triangle, no diagonal
                pval(ltri) = pvalSpearman(tail,D(ltri),meanD,stdD,n);
                pval = pval + pval';
                % The p-values along the diagonal are always exactly one
                % (unless they're NaN), because Pr{coef(i,i) as/more extreme
                % than 1} == 1, regardless of which tail(s) we're testing
                % against.
                pval(1:p1+1:end) = sign(diag(coef)); % preserves NaNs on diag
            else
                pval = pvalSpearman(tail,D,meanD,stdD,n);
            end
        end

    % Compute one row at a time when we need p-values and there are ties
    else
        coef = zeros(p1,p2,outClass);
        pval = zeros(p1,p2,outClass);
        for i = 1:p1
            xranki = xrank(:,i);
            if corrXX
                j = 1:i; % lower triangle and diagonal
                yrankj = xrank(:,j);
                yadjj = yadj(j);
            else
                j = 1:p2;
                yrankj = yrank;
                yadjj = yadj;
            end
            D = sum(bsxfun(@minus,xranki,yrankj).^2); % sum((xranki - yrankj).^2);
            meanD = (n3const - (xadj(i)+yadjj)./3) ./ 2;
            stdD = sqrt((n3const./2 - xadj(i)./3)*(n3const./2 - yadjj./3)./(n-1));

            % When either of the data vectors is constant, stdD should be
            % zero, get that exactly.
            n3const2 = (n+1)*n*(n-1)/2;
            stdD((xadj(i) == n3const2) | (yadjj == n3const2)) = 0;

            coef(i,j) = (meanD - D) ./ (sqrt(n-1)*stdD);
            pval(i,j) = pvalSpearman(tail,D,meanD,stdD,xranki,yrankj);
        end

        % Limit off-diag correlations to [-1,1].
        t = find(abs(coef) > 1); coef(t) = coef(t)./abs(coef(t)); % preserves NaNs

        % If this is autocorrelation, reflect the lower triangle into the upper.
        if corrXX
            coef = tril(coef) + tril(coef,-1)'; % leave diagonal alone
            % The p-values along the diagonal are always exactly one (unless
            % they're NaN), because Pr{coef(i,i) as/more extreme than 1} == 1,
            % regardless of which tail(s) we're testing against.
            pval = pval + pval';
            pval(1:p1+1:end) = sign(diag(coef)); % preserve NaNs on diag
        end
    end
end
end

%--------------------------------------------------------------------------

function p = pvalPearson(tail, rho, n)
%PVALPEARSON Tail probability for Pearson's linear correlation.
t = rho.*sqrt((n-2)./(1-rho.^2)); % +/- Inf where rho == 1
switch tail
case 'b' % 'both or 'ne'
    p = 2*tcdf(-abs(t),n-2);
case 'r' % 'right' or 'gt'
    p = tcdf(-t,n-2);
case 'l' % 'left' or 'lt'
    p = tcdf(t,n-2);
end
end

%--------------------------------------------------------------------------

function p = pvalKendall(tail, K, stdK, arg1, arg2)
%PVALKENDALL Tail probability for Kendall's K statistic.

% Without ties, K is symmetric about zero, taking on values in
% -n(n-1)/2:2:n(n-1)/2.  With ties, it's still in that range, but not
% symmetric, and can take on adjacent integer values.

% K and stdK may be vectors when n is given (i.e., when there were no ties).
% When xrank and yrank are given (i.e., when there were ties), K and stdK must
% be scalars and xrank and yrank must both be a single column.

if nargin < 5 % pvalKendall(tail, K, stdK, n), no ties
    noties = true;
    n = arg1;
    exact = (n < 50);
else % pvalKendall(tail, K, stdK, xrank, yrank), ties in data
    noties = false;

    % If stdK is zero, at least one of the data vectors was constant.  The
    % correlation coef in these cases falls out of the calculations correctly
    % as NaN, but the exact p-value calculations below would regard Pr{K==0}
    % as 1.  Return a NaN p-value instead.
    K(stdK == 0) = NaN;

    xrank = arg1;
    yrank = arg2;
    n = length(xrank);
    exact = (n < 10);
end
nfact = factorial(n);
n2const = n*(n-1)/2;

if exact
    if noties
        % No ties, use recursion to get the cumulative distribution of
        % the number, C, of positive (xi-xj)*(yi-yj), i<j.
        %
        % K = #pos-#neg = C-Q, and C+Q = n(n-1)/2 => C = (K + n(n-1)/2)/2
        freq = [1 1];
        for i = 3:n
            freq = conv(freq,ones(1,i));
        end
        freq = [freq; zeros(1,n2const+1)]; freq = freq(1:end-1)';

    else
        % Ties, take permutations of the midranks.
        %
        % With ties, we could consider only distinguishable permutations
        % (those for which equal ranks (ties) are not simply interchanged),
        % but generating only those is a bit of work.  Generating all
        % permutations uses more memory, but gives the same result.
        xrank = repmat(xrank(:)',nfact,1);
        yrank = perms(yrank(:)');
        Kperm = zeros(nfact,1);
        for k = 1:n-1
            U = sign(repmat(xrank(:,k),1,n-k)-xrank(:,k+1:n));
            V = sign(repmat(yrank(:,k),1,n-k)-yrank(:,k+1:n));
            Kperm = Kperm + sum(U .* V, 2);
        end
        freq = histc(Kperm,-(n2const+.5):(n2const+.5)); freq = freq(1:end-1);
    end

    % Get the tail probabilities.  Reflect as necessary to get the correct
    % tail.
    switch tail
    case 'b' % 'both or 'ne'
        % Use twice the smaller of the tail area above and below the
        % observed value.
        tailProb = min(cumsum(freq), rcumsum(freq,nfact)) ./ nfact;
        tailProb = min(2*tailProb, 1); % don't count the center bin twice
    case 'r' % 'right' or 'gt'
        tailProb = rcumsum(freq,nfact) ./ nfact;
    case 'l' % 'left' or 'lt'
        tailProb = cumsum(freq) ./ nfact;
    end
    p = NaN(size(K),class(K));
    t = ~isnan(K(:));
    p(t) = tailProb(K(t) + n2const+1); % bins at integers, starting at -n2const

else
    switch tail
    case 'b' % 'both or 'ne'
        p = normcdf(-(abs(K)-1) ./ stdK);
        p = 2*p; p(p>1) = 1; % Don't count continuity correction at center twice
    case 'r' % 'right' or 'gt'
        p = normcdf(-(K-1) ./ stdK);
    case 'l' % 'left' or 'lt'
        p = normcdf((K+1) ./ stdK);
    end
end
end

%--------------------------------------------------------------------------

function p = pvalSpearman(tail, D, meanD, stdD, arg1, arg2)
%PVALSPEARMAN Tail probability for Spearman's D statistic.

% Without ties, D is symmetric about (n^3-n)/6, taking on values
% 0:2:(n^3-n)/3.  With ties, it's still in (but not on) that range, but
% not symmetric, and can take on odd and half-integer values.

% D, meanD, and stdD may be vectors when n is given (i.e., when there were no
% ties). When xrank and yrank are given (i.e., when there were ties), yrank
% may be a matrix, but xrank must be a single column, and D, meanD, stdD must
% have the same lengths as yrank has columns.

if nargin < 6 % pvalSpearman(tail, D, meanD, stdD, n), no ties
    noties = true;
    n = arg1;
else % pvalSpearman(tail, D, meanD, stdD, xrank, yrank), ties in data
    noties = false;

    % If stdD is zero, at least one of the data vectors was constant.  The
    % correlation coef in these cases falls out of the calculations correctly
    % as NaN, but the exact p-value calculations below would regard Pr{D==D0}
    % as 1.  Return a NaN p-value instead.
    D(stdD == 0) = NaN;

    xrank = arg1;
    yrank = arg2;
    n = length(xrank);
end
exact = (n < 10);
if exact
    p = NaN(size(D),class(D));
    if noties
        tailProb = spearmanExactSub(tail,n);
        t = ~isnan(D(:));
        p(t) = tailProb(2*D(t)+1); % bins at half integers, starting at zero
    else
        for j = 1:size(yrank,2)
            if ~isnan(D(j))
                tailProb = spearmanExactSub(tail,xrank,yrank(:,j));
                p(j) = tailProb(2*D(j)+1); % bins at half integers, starting at zero
            end
        end
    end

else
    if noties
        % Use AS89, an Edgeworth expansion for upper tail prob of D.
        n3const = (n^3 - n)/3;
        switch tail
        case 'b' % 'both or 'ne'
            p = AS89(max(D, n3const-D), n, n3const);
            p = 2*p; p(p>1) = 1; % Don't count continuity correction at center twice
        case 'r' % 'right' or 'gt'
            p = AS89(n3const - D, n, n3const);
        case 'l' % 'left' or 'lt'
            p = AS89(D, n, n3const);
        end
    else
        % Use a t approximation.
        r = (meanD - D) ./ (sqrt(n-1)*stdD);
        t = Inf*sign(r);
        ok = (abs(r) < 1);
        t(ok) = r(ok) .* sqrt((n-2)./(1-r(ok).^2));
        switch tail
        case 'b' % 'both or 'ne'
            p = 2*tcdf(-abs(t),n-2);
        case 'r' % 'right' or 'gt'
            p = tcdf(-t,n-2);
        case 'l' % 'left' or 'lt'
            p = tcdf(t,n-2);
        end
    end
end
end

%--------------------------------------------------------------------------

function tailProb = spearmanExactSub(tail,arg1,arg2)
if nargin < 3
    % No ties, take permutations of 1:n
    n = arg1;
    nfact = factorial(n);
    Dperm = sum((repmat(1:n,nfact,1) - perms(1:n)).^2, 2);
else
    % Ties, take permutations of the midranks.
    %
    % With ties, we could consider only distinguishable permutations
    % (those for which equal ranks (ties) are not simply interchanged),
    % but generating only those is a bit of work.  Generating all
    % permutations uses more memory, but gives the same result.
    xrank = arg1;
    yrank = arg2;
    n = length(xrank);
    nfact = factorial(n);
    Dperm = sum((repmat(xrank(:)',nfact,1) - perms(yrank(:)')).^2, 2);
end
n3const = (n^3 - n)/3;
freq = histc(Dperm,(-.25):.5:(n3const+.25)); freq = freq(1:end-1);

% Get the tail probabilities.  Reflect as necessary to get the correct
% tail: the left tail of D corresponds to right tail of rho.
switch tail
case 'b' % 'both or 'ne'
    % Use twice the smaller of the tail area above and below the
    % observed value.
    tailProb = min(cumsum(freq), rcumsum(freq,nfact)) ./ nfact;
    tailProb = min(2*tailProb, 1); % don't count the center bin twice
case 'r' % 'right' or 'gt'
    tailProb = cumsum(freq) ./ nfact;
case 'l' % 'left' or 'lt'
    tailProb = rcumsum(freq,nfact) ./ nfact;
end
end

%--------------------------------------------------------------------------

function p = AS89(D, n, n3const)
%AS89 Upper tail probability for Spearman's D statistic.
%   Edgeworth expansion for the upper tail probability of D in the no ties
%   case, with continuity correction, adapted from Applied Statistics
%   algorithm AS89.
c = [.2274 .2531 .1745 .0758 .1033 .3932 .0879 .0151 .0072 .0831 .0131 .00046];
x = (2*(D-1)./n3const - 1) * sqrt(n - 1);
y = x .* x;
u = x .* (c(1) + (c(2) + c(3)/n)/n + ...
    y .* (-c(4) + (c(5) + c(6)/n)/n - ...
    y .* (c(7) + c(8)/n - y .* (c(9) - c(10)/n + ...
    y .* (c(11) - c(12) * y)/n))/n))/n;
p = u ./ exp(.5*y) + 0.5 * erfc(x./sqrt(2));
p(p>1) = 1; p(p<0) = 0; % don't ignore NaNs
end

%--------------------------------------------------------------------------

function y = rcumsum(x,sumx)
%RCUMSUM Cumulative sum in reverse direction.
if nargin < 2
    sumx = sum(x);
end
y = repmat(sumx,size(x));
y(2:end) = sumx - cumsum(x(1:end-1));
end

function [z,mu,sigma] = zscore(x,flag,dim)
%ZSCORE Standardized z score.
%   Z = ZSCORE(X) returns a centered, scaled version of X, the same size as X.
%   For vector input X, Z is the vector of z-scores (X-MEAN(X)) ./ STD(X). For
%   matrix X, z-scores are computed using the mean and standard deviation
%   along each column of X.  For higher-dimensional arrays, z-scores are
%   computed using the mean and standard deviation along the first
%   non-singleton dimension.
%
%   The columns of Z have sample mean zero and sample standard deviation one
%   (unless a column of X is constant, in which case that column of Z is
%   constant at 0).
%
%   [Z,MU,SIGMA] = ZSCORE(X) also returns MEAN(X) in MU and STD(X) in SIGMA.
%
%   [...] = ZSCORE(X,1) normalizes X using STD(X,1), i.e., by computing the
%   standard deviation(s) using N rather than N-1, where N is the length of
%   the dimension along which ZSCORE works.  ZSCORE(X,0) is the same as
%   ZSCORE(X).
%
%   [...] = ZSCORE(X,FLAG,'all') standardizes X by working on all the
%   elements of X. Pass in FLAG==0 to use the default normalization by N-1,
%   or 1 to use N.
%
%   [...] = ZSCORE(X,FLAG,DIM) standardizes X by working along the dimension
%   DIM of X.
%
%   [...] = ZSCORE(X,FLAG,VECDIM) standardizes X by working along the all the
%   dimensions of X specified in VECDIM.
%
%   See also MEAN, STD.

%   Copyright 1993-2018 The MathWorks, Inc. 
% [] is a special case for std and mean, just handle it out here.
if isequal(x,[]), z = x; return; end

if nargin < 2
    flag = 0;
end
if nargin < 3
    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Compute X's mean and sd, and standardize it
mu = mean(x,dim);
sigma = std(x,flag,dim);
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z, sigma0);
end


function plotPhaseVectors(axhandle,theta,tax,pax,tspace,pspace)
if ~isempty(findobj(axhandle,'type','patch'))
    delete(findobj(axhandle, 'type', 'patch'));
end

[tgrid,pgrid]=meshgrid(tax,log2(pax));
theta = theta(1:pspace:size(theta,1),1:tspace:size(theta,2));

idx = find(~any(isnan([tgrid(:) pgrid(:) theta(:)]),2));

tgrid = tgrid(idx);
pgrid = pgrid(idx);
theta = theta(idx);

% Determine extent of phase arrows in plot
[dx,dy] = determinearrowextent(axhandle);
%

% Create the arrow patch object for plotting the phase
arrowpatch = [-1 0 0 1 0 0 -1; 0.1 0.1 0.5 0 -0.5 -0.1 -0.1]';

for ii=numel(tgrid):-1:1
    % Multiply each arrow by the rotation matrix for the given theta
    rotarrow = arrowpatch*[cos(theta(ii)) sin(theta(ii));...
        -sin(theta(ii)) cos(theta(ii))];
    patch(tgrid(ii)+rotarrow(:,1)*dx,pgrid(ii)+rotarrow(:,2)*dy,[0 0 0],...
        'edgecolor','none' ,'Parent',axhandle);
end
end
