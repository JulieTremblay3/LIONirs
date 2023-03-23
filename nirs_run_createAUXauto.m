function out = nirs_run_createAUXauto(job)

load(job.NIRSmat{1},'-mat');
NC = NIRS.Cf.H.C.N;
rDtp = NIRS.Dt.fir.pp(end).p;
fs = NIRS.Cf.dev.fs;
if isfield(NIRS.Dt,'AUX')
    if job.m_replaceaux % if one add aux to the existing one
        idAUX = numel(NIRS.Dt.AUX)+1; %new fields for AUX
    else  %remove old field and create new
        NIRS.Dt = rmfield(NIRS.Dt,'AUX');
        idAUX = 1; %new fields for AUX
        NIRS.Dt.AUX(idAUX).pp.p{size(rDtp,1),1}=[];
        NIRS.Dt.AUX(idAUX).pp.sync_timesec{size(rDtp,1),1}=[];
        NIRS.Dt.AUX(idAUX).label = [];
    end
else
    idAUX = 1; %new fields for AUX
    NIRS.Dt.AUX(idAUX).pp.p{size(rDtp,1),1}=[];
    NIRS.Dt.AUX(idAUX).pp.sync_timesec{size(rDtp,1),1}=[];
    NIRS.Dt.AUX(idAUX).label = [];
end
onsetall = [];
fileall = [];
for ifile =1:numel(rDtp)
    %[num, txt, raw] = xlsread(fileevent)
    filedata = rDtp{ifile};
    d1 = fopen_NIR(rDtp{ifile,1},NC);
    Ns = size(d1,2); clear d1;
    [EEG.data,EEG.infoBV,EEG.marker,EEG.ind_dur_ch]= fopen_EEG(filedata);
    idtrigger = find(strcmp('trigger',EEG.marker(:,1))& EEG.ind_dur_ch(:,1)>0);
    EEG.marker = EEG.marker(idtrigger,:);
    EEG.ind_dur_ch = EEG.ind_dur_ch(idtrigger,:);
    if isfield(job.c_createAUXauto,'b_HRFtriggeronset') % use trigger as HRF onset
        trigtofind = job.c_createAUXauto.b_HRFtriggeronset.e_HRFtrigger;
        aux5 = NIRS.Dt.fir.aux5{ifile};
        onset = []; %initialise in case no trig available
        dur = [];
        weight = [];
        for itrig = 1:numel(trigtofind) %multiple trig support
            id=  find(sum(aux5(:,1)==trigtofind(itrig),2));
            if ~isempty(id)
                onset =[onset;aux5(id,2).*1./fs];
                dur = [dur;repmat(job.c_createAUXauto.b_HRFtriggeronset.e_HRFduration,numel(onset),1)] ;
                weight =[weight;repmat(1,numel(onset),1)] ;
            end
        end
        eventV =zeros(size(EEG.data,1),1);
        fs = 1/(EEG.infoBV.SamplingInterval/1e6);
        t = 1/fs:1/fs:1/fs*size(EEG.data,1);
        idevent = [];
        for ievent = 1:numel(onset)
            idon = find(t>onset(ievent) & t<(onset(ievent)+dur(ievent)));
            eventV(idon ) = weight(ievent,1);
        end
        
        if ~isempty(dur)
            frametimes = 1/fs:1/fs:(dur(1)+30);
            param1=job.c_createAUXauto.b_HRFtriggeronset.e_TimetoPeak1;
            param2=job.c_createAUXauto.b_HRFtriggeronset.e_FWHM1;
            param3=job.c_createAUXauto.b_HRFtriggeronset.e_TimetoPeak2;
            param4=job.c_createAUXauto.b_HRFtriggeronset.e_FWHM2;
            param5=job.c_createAUXauto.b_HRFtriggeronset.e_DIP;
            %  [5.4 5.2 10.8 7.35 0.35]
            H = fmridesign(frametimes,0,[1 0],[],[param1 param2 param3 param4 param5 0]); %Default
            %figure;plot(H.X(:,1,1))
            Sc1 = conv2(H.X(:,1,1), eventV);
            Scntmp = Sc1(1:end-numel(H.X(:,1,1))+1);
            Scntmp = Scntmp/max(Scntmp);
            %      Regressor.xdata = t';
            %      Regressor.ydata = Scn;
            %      Onset.xdata = [onset,onset];
            %      Onset.ydata = [min(Scn),max(Scn) ];
            if numel(Scntmp)>=Ns
                Scn = Scntmp(1:Ns);
            else
                Scn= zeros(Ns,1);
                Scn(1:numel(Scntmp)) = Scntmp;
            end
        else
            Scn = zeros(Ns,1);
        end
        
        %figure;plot(eventV)
        filename = 'AUX';
        HRFcurve.data = Scn;
        HRFcurve.ind_dur_ch=EEG.ind_dur_ch;
        HRFcurve.marker=EEG.marker;
        HRFcurve.infoBV = EEG.infoBV;
        label = job.c_createAUXauto.b_HRFtriggeronset.e_HRFlabel;
        HRFcurve.infoBV.DataType = 'TIMEDOMAIN';
        HRFcurve.infoBV.dataformat = 'BINARY';
        HRFcurve.infoBV.DataOrientation='VECTORIZED';
        HRFcurve.infoBV.NumberOfChannels = '1';
        HRFcurve.infoBV.name_ele = {deblank(label)};
        HRFcurve.infoBV.coor_r = 1;
        HRFcurve.infoBV.coor_theta  = 90;
        HRFcurve.infoBV.coor_phi = 45;
        HRFcurve.infoBV.DataPoints = numel(Scn);
        [pathstr, name, ext] = fileparts(filedata);
        fileoutHRFcurve = fullfile(job.c_createAUXauto.b_HRFtriggeronset.e_AUXdir{1}, [label,name,'.dat']);
        disp(['File HRF AUX create: ', fileoutHRFcurve]);
        fwrite_EEG(fileoutHRFcurve,HRFcurve,1,HRFcurve.infoBV.DataPoints );
        NIRS.Dt.AUX(idAUX).pp(1,numel(NIRS.Dt.AUX(1).pp)).p{ifile,1}=fileoutHRFcurve;
        NIRS.Dt.AUX(idAUX).pp(1,numel(NIRS.Dt.AUX(1).pp)).sync_timesec{ifile,1}=0;
        NIRS.Dt.AUX(idAUX).label = label;
        onsetall = [onsetall;onset];
        fileall = [fileall; repmat(ifile,numel(onset),1)];
           e_HRFduration =   dur;
    elseif isfield(job.c_createAUXauto, 'b_HRFxlsonset') %use xls file to define trig onset
        fileevent = job.c_createAUXauto.b_HRFxlsonset.e_HRFxlsfiles{1} ;
        aux5 = NIRS.Dt.fir.aux5{ifile};        
        [~,~,ext] =fileparts(fileevent);
        if strcmp(ext,'.xlsx')|strcmp(ext,'.xls')
            [num, txt, raw] = xlsread(fileevent);
        elseif strcmp(ext,'.txt')   
            [num, txt, raw] = readtxtfile_asxlsread(fileevent);
        end
  
    %HRF convolution
    onset =num(:,1);
    dur = num(:,2);
    weight = num(:,3);
    eventV =zeros(size(EEG.data,1),1);
    fs = 1/(EEG.infoBV.SamplingInterval/1e6);
    t = 1/fs:1/fs:1/fs*size(EEG.data,1);
    idevent = [];
    for ievent = 1:numel(onset)
        idon = find(t>onset(ievent) & t<(onset(ievent)+dur(ievent)));      
        eventV(idon) = weight(ievent,1);
    end
    
        if ~isempty(dur)
            frametimes = 1/fs:1/fs:(dur(1)+30);
            param1=job.c_createAUXauto.b_HRFxlsonset.e_TimetoPeak1;
            param2=job.c_createAUXauto.b_HRFxlsonset.e_FWHM1;
            param3=job.c_createAUXauto.b_HRFxlsonset.e_TimetoPeak2;
            param4=job.c_createAUXauto.b_HRFxlsonset.e_FWHM2;
            param5=job.c_createAUXauto.b_HRFxlsonset.e_DIP;
            %  [5.4 5.2 10.8 7.35 0.35]
            H = fmridesign(frametimes,0,[1 0],[],[param1 param2 param3 param4 param5 0]); %Default
            Sc1 = conv2(H.X(:,1,1), eventV);
            Scntmp = Sc1(1:end-numel(H.X(:,1,1))+1);
            Scntmp = Scntmp/max(Scntmp);
            %      Regressor.xdata = t';
            %      Regressor.ydata = Scn;
            %      Onset.xdata = [onset,onset];
            %      Onset.ydata = [min(Scn),max(Scn) ];
            if numel(Scntmp)>=Ns
                Scn = Scntmp(1:Ns);
            else
                Scn= zeros(Ns,1);
                Scn(1:numel(Scntmp)) = Scntmp;
            end
        else
            Scn = zeros(Ns,1);
        end

        H = fmridesign(frametimes,0,[1 0],[],[param1 param2 param3 param4 param5 0]); %Default
        Sc1 = conv2(H.X(:,1,1), eventV);
        Scn = Sc1(1:end-numel(H.X(:,1,1))+1);

        Scn = Scn/max(Scn);
         %figure;plot(eventV)
        filename = 'AUX';
        HRFcurve.data = Scn;
        HRFcurve.ind_dur_ch=EEG.ind_dur_ch;
        HRFcurve.marker=EEG.marker;
        HRFcurve.infoBV = EEG.infoBV;
        label = job.c_createAUXauto.b_HRFxlsonset.e_HRFlabel;
        HRFcurve.infoBV.DataType = 'TIMEDOMAIN';
        HRFcurve.infoBV.dataformat = 'BINARY';
        HRFcurve.infoBV.DataOrientation='VECTORIZED';
        HRFcurve.infoBV.NumberOfChannels = '1';
        HRFcurve.infoBV.name_ele = {deblank(label)};
        HRFcurve.infoBV.coor_r = 1;
        HRFcurve.infoBV.coor_theta  = 90;
        HRFcurve.infoBV.coor_phi = 45;
        HRFcurve.infoBV.DataPoints = numel(Scn);
        [pathstr, name, ext] = fileparts(filedata);
        fileoutHRFcurve = fullfile(job.c_createAUXauto.b_HRFxlsonset.e_AUXdir{1}, [label,name,'.dat']);
        disp(['File HRF AUX create: ', fileoutHRFcurve]);
        fwrite_EEG(fileoutHRFcurve,HRFcurve,1,HRFcurve.infoBV.DataPoints );
        NIRS.Dt.AUX(idAUX).pp(1,numel(NIRS.Dt.AUX(1).pp)).p{ifile,1}=fileoutHRFcurve;
        NIRS.Dt.AUX(idAUX).pp(1,numel(NIRS.Dt.AUX(1).pp)).sync_timesec{ifile,1}=0;
        NIRS.Dt.AUX(idAUX).label = label;
        onsetall = [onsetall;onset];
        fileall = [fileall; repmat(ifile,numel(onset),1)];
        e_HRFduration = num(:,2);
    end
end
%create defaults xls to help building extract configuration.

[filepath,name,ext] = fileparts(job.NIRSmat{1});
filenamexls = fullfile(filepath,'ExtractHRF.xlsx');
try
    [num,txt,raw] = xlsread(filenamexls)
catch
    raw = [];
end
e_HRFduration
    if ~isfield(job.c_createAUXauto.b_HRFtriggeronset,'e_HRF_SDmodel')
        A = {'NIRS.mat folder','File','Trig', 'tStart','tStop','label','X0'};
        VAL = [repmat(job.NIRSmat,numel(onsetall),1),  num2cell(fileall), num2cell(onsetall), num2cell(onsetall-5),...
            num2cell(onsetall+repmat(e_HRFduration(1),numel(onsetall),1)++10), repmat({label},numel(onsetall),1) , repmat({label},numel(onsetall),1)  ];
        if isempty(raw);raw = [A;VAL];else; raw = [raw;VAL];end
    else isfield(job.c_createAUXauto.b_HRFtriggeronset,'e_HRF_SDmodel') 
        if isempty(job.c_createAUXauto.b_HRFtriggeronset.e_HRF_SDmodel{1}) %HRF only
             [pathdefault, ~,~]= fileparts(NIRS.Dt.fir.pp(1).p{1});
             disp(['Use default :',fullfile(pathdefault,'Global.zone'),' in the regression'])
             job.c_createAUXauto.b_HRFtriggeronset.e_HRF_SDmodel{1} = fullfile(pathdefault,'Global.zone');
        end
        if strcmp(job.c_createAUXauto.b_HRFtriggeronset.e_HRF_SDmodel, 'No')
            A = {'NIRS.mat folder','File','Trig', 'tStart','tStop','label','X0'};
            VAL = [repmat(job.NIRSmat,numel(onsetall),1),  num2cell(fileall), num2cell(onsetall), num2cell(onsetall-5),...
            num2cell(onsetall+repmat(e_HRFduration(1),numel(onsetall),1)+10), repmat({label},numel(onsetall),1) , repmat({label},numel(onsetall),1)  ];
            if isempty(raw);raw = [A;VAL];else; raw = [raw;VAL];end

        else
           	A = {'NIRS.mat folder','File','Trig', 'tStart','tStop','label','X0', 'X1'};
            nameSDfile =job.c_createAUXauto.b_HRFtriggeronset.e_HRF_SDmodel{1}; %'C:\data\FRESH\Analyzed\STAT_Motor\MotorSD.zone'
            VAL = [repmat(job.NIRSmat,numel(onsetall),1),  num2cell(fileall), num2cell(onsetall), num2cell(onsetall-5),...
            num2cell(onsetall+repmat(e_HRFduration(1),numel(onsetall),1)+10), repmat({label},numel(onsetall),1) , repmat({label},numel(onsetall),1),repmat({nameSDfile},numel(onsetall),1)  ];
            if isempty(raw);raw = [A;VAL];else; raw = [raw;VAL];end
        end
    end

try
    xlswrite(filenamexls,raw)
    disp(['Create default xls config to GLM ', label, ' extract to visualized open : winopen ', filenamexls ]);
catch
    writetxtfile(filenamexls,raw);
    disp(['Create default xls config to GLM ', label, ' extract to visualized open : winopen ', filenamexls ]);
end

save(job.NIRSmat{1},'NIRS');

out.NIRSmat = job.NIRSmat;


function X_cache=fmridesign(frametimes,slicetimes,events,S, ...
    hrf_parameters,shift)

%FMRIDESIGN makes a set of design matrices for fmrilm.
%
% X_CACHE = FMRIDESIGN( FRAME_TIMES [, SLICE_TIMES [, EVENTS , [S ,
%                       [, HRF_PARAMETERS [, SHIFT]]]]]] )
%
% FRAME_TIMES is a row vector of frame acquisition times in seconds.
% With just the frametimes, it gives the hemodynamic response function.
%
% SLICE_TIMES is a row vector of relative slice acquisition times,
% i.e. absolute acquisition time of a slice is FRAME_TIMES + SLICE_TIMES.
% Default is 0.
%
% EVENTS is a matrix whose rows are events and whose columns are:
% 1. id - an integer from 1:(number of events) to identify event type;
% 2. times - start of event, synchronised with frame and slice times;
% 3. durations (optional - default is 0) - duration of event;
% 4. heights (optional - default is 1) - height of response for event.
% For each event type, the response is a box function starting at the event
% times, with the specified durations and heights, convolved with the
% hemodynamic response function (see below). If the duration is zero, the
% response is the hemodynamic response function whose integral is
% the specified height - useful for `instantaneous' stimuli such as visual
% stimuli. The response is then subsampled at the appropriate frame and slice
% times to create a design matrix for each slice, whose columns correspond
% to the event id number. EVENT_TIMES=[] will ignore event times and just
% use the stimulus design matrix S (see next). Default is [1 0].
%
% S: Events can also be supplied by a stimulus design matrix,
% whose rows are the frames, and column are the event types. Events
% are created for each column, beginning at the frame time for each row
% of S, with a duration equal to the time to the next frame, and a height
% equal to the value of S for that row and column. Note that a
% constant term is not usually required, since it is removed by the
% polynomial trend terms provided N_POLY>=0. Default is [].
%
% HRF_PARAMETERS is a matrix whose rows are 5 parameters for the
% hemodynamic response function, one row for each event type and column
% of S (if there is just one row, this is repeated as necessary).
% The hrf is modeled as the difference of two gamma density functions
% (Glover, NeuroImage, 9:416-429). The components of HRF_PARAMETERS are:
% 1. PEAK1: time to the peak of the first gamma density;
% 2. FWHM1: approximate FWHM of the first gamma density;
% 3. PEAK2: time to the peak of the second gamma density;
% 4. FWHM2: approximate FWHM of the second gamma density;
% 5. DIP: coefficient of the second gamma density;
%    Final hrf is:   gamma1/max(gamma1)-DIP*gamma2/max(gamma2)
%    scaled so that its total integral is 1.
% If PEAK1=0 then there is no smoothing of that event type with the hrf.
% If PEAK1>0 but FWHM1=0 then the design is simply lagged by PEAK1.
% Default is: [5.4 5.2 10.8 7.35 0.35] chosen by Glover (1999) for
% an auditory stimulus.
% If HRF_PARAMETERS is a structure, then HRF_PARAMETERS.T is a matrix
% whose rows are the times in seconds of a user-supplied HRF, one
% row for each event type and column of S (if there is just one row,
% this is repeated as necessary).  Times must start at 0 but need not
% be equally spaced; spacing of 0.02s is recommended. HRF_PARAMETERS.H
% are the corresponding values of the HRF at those times.
%
% SHIFT is a matrix whose rows are the min and max shift for the hrf in
% seconds, one row for each event type and column of S (if there is just
% one row, this is repeated as necessary). Default is [-4.5 4.5]*FWHM1/5.2.
%
% X_CACHE.TR: TR, average time between frames (secs).
%
% X_CACHE.X: A cache of the design matrices stored as a 4D array.
% Dim 1: frames; Dim 2: response variables; Dim 3: 4 values, corresponding
% to the stimuli convolved with: hrf, derivative of hrf, first and second
% spectral basis functions over the range in SHIFT; Dim 4: slices.
%
% X_CACHE.W: A 3D array of coefficients of the basis functions in X_CACHE.X.
% Dim 1: frames; Dim 2: response variables; Dim 3: 5 values:
% coefficients of the hrf and its derivative, coefficients of the first and
% second spectral basis functions, shift values.

%############################################################################
% COPYRIGHT:   Copyright 2002 K.J. Worsley
%              Department of Mathematics and Statistics,
%              McConnell Brain Imaging Center,
%              Montreal Neurological Institute,
%              McGill University, Montreal, Quebec, Canada.
%              worsley@math.mcgill.ca, liao@math.mcgill.ca
%
%              Permission to use, copy, modify, and distribute this
%              software and its documentation for any purpose and without
%              fee is hereby granted, provided that the above copyright
%              notice appear in all copies. The author and McGill University
%              make no representations about the suitability of this
%              software for any purpose.  It is provided "as is" without
%              express or implied warranty.
%############################################################################

% Defaults:

if nargin < 2
    slicetimes=0;
end
if nargin < 3
    events=[1 0];
end
if nargin < 4
    S=[];
end
if nargin < 5
    hrf_parameters=[5.4 5.2 10.8 7.35 0.35];
end
if nargin < 6
    if isstruct(hrf_parameters)
        shift=[-4.5 4.5];
    else
        shift=[-4.5 4.5]*max(max(hrf_parameters(:,2))/5.2,1);
    end
end

n=length(frametimes);
numslices=length(slicetimes);

% Keep time points that are not excluded:

if ~isempty(events)
    numevents=size(events,1);
    eventid=events(:,1);
    numeventypes=max(eventid);
    eventime=events(:,2);
    if size(events,2)>=3
        duration=events(:,3);
    else
        duration=zeros(numevents,1);
    end
    if size(events,2)>=4
        height=events(:,4);
    else
        height=ones(numevents,1);
    end
    mineventime=min(eventime);
    maxeventime=max(eventime+duration);
else
    numeventypes=0;
    mineventime=Inf;
    maxeventime=-Inf;
end

if ~isempty(S)
    numcolS=size(S,2);
else
    numcolS=0;
end

% Set up response matrix:

dt=0.02;
startime=min(mineventime,min(frametimes)+min([slicetimes 0]));
finishtime=max(maxeventime,max(frametimes)+max([slicetimes 0]));
numtimes=ceil((finishtime-startime)/dt)+1;
numresponses=numeventypes+numcolS;
response=zeros(numtimes,numresponses);

if ~isempty(events)
    height=height./(1+(duration==0)*(dt-1));
    for k=1:numevents
        type=eventid(k);
        n1=ceil((eventime(k)-startime)/dt)+1;
        n2=ceil((eventime(k)+duration(k)-startime)/dt)+(duration(k)==0);
        if n2>=n1
            response(n1:n2,type)=response(n1:n2,type)+height(k)*ones(n2-n1+1,1);
        end
    end
end

if ~isempty(S)
    for j=1:numcolS
        for i=find(S(:,j)')
            n1=ceil((frametimes(i)-startime)/dt)+1;
            if i<n
                n2=ceil((frametimes(i+1)-startime)/dt);
            else
                n2=numtimes;
            end
            if n2>=n1
                response(n1:n2,numeventypes+j)= ...
                    response(n1:n2,numeventypes+j)+S(i,j)*ones(n2-n1+1,1);
            end
        end
    end
end

if isstruct(hrf_parameters)
    if size(hrf_parameters.T,1)==1
        hrf_parameters.T=repmat(hrf_parameters.T,numresponses,1);
    end
    if size(hrf_parameters.H,1)==1
        hrf_parameters.H=repmat(hrf_parameters.H,numresponses,1);
    end
else
    if size(hrf_parameters,1)==1
        hrf_parameters=repmat(hrf_parameters,numresponses,1);
    end
end
if size(shift,1)==1
    shift=repmat(shift,numresponses,1);
end

eventmatrix=zeros(numtimes,numresponses,4);
nd=41;
X_cache.W=zeros(nd,numresponses,5);

for k=1:numresponses
    Delta1=shift(k,1);
    Delta2=shift(k,2);
    if isstruct(hrf_parameters)
        numlags=ceil((max(hrf_parameters.T(k,:))+Delta2-Delta1)/dt)+1;
    else
        peak1=hrf_parameters(k,1);
        fwhm1=hrf_parameters(k,2);
        peak2=hrf_parameters(k,3);
        fwhm2=hrf_parameters(k,4);
        dip=hrf_parameters(k,5);
        numlags=ceil((max(peak1+3*fwhm1,peak2+3*fwhm2)+Delta2-Delta1)/dt)+1;
    end
    numlags=min(numlags,numtimes);
    time=(0:(numlags-1))'*dt;
    
    % Taylor:
    if isstruct(hrf_parameters)
        hrf=interp1(hrf_parameters.T(k,:),hrf_parameters.H(k,:),time,'spline',0);
        d_hrf=-gradient(hrf,dt);
    else
        tinv=(time>0)./(time+(time<=0));
        if peak1>0 & fwhm1>0
            alpha1=peak1^2/fwhm1^2*8*log(2);
            beta1=fwhm1^2/peak1/8/log(2);
            gamma1=(time/peak1).^alpha1.*exp(-(time-peak1)./beta1);
            d_gamma1=-(alpha1*tinv-1/beta1).*gamma1;
        else
            gamma1=min(abs(time-peak1))==abs(time-peak1);
            d_gamma1=zeros(numlags,1);
        end
        if peak2>0 & fwhm2>0
            alpha2=peak2^2/fwhm2^2*8*log(2);
            beta2=fwhm2^2/peak2/8/log(2);
            gamma2=(time/peak2).^alpha2.*exp(-(time-peak2)./beta2);
            d_gamma2=-(alpha2*tinv-1/beta2).*gamma2;
        else
            gamma2=min(abs(time-peak2))==abs(time-peak2);
            d_gamma2=zeros(numlags,1);
        end
        hrf=gamma1-dip*gamma2;
        d_hrf=d_gamma1-dip*d_gamma2;
    end
    HS=[hrf d_hrf]/sum(hrf);
    temp=conv2(response(:,k),HS);
    eventmatrix(:,k,1:2)=temp(1:numtimes,:);
%figure;plot(d_hrf)
    % Shifted hrfs:
    H=zeros(numlags,nd);
    delta=((1:nd)-1)/(nd-1)*(Delta2-Delta1)+Delta1;
    for id=1:nd
        if isstruct(hrf_parameters)
            t=time+Delta1-delta(id);
            hrf=interp1(hrf_parameters.T(k,:),hrf_parameters.H(k,:),t,'spline',0);
        else
            t=(time+Delta1-delta(id)).*((time+Delta1)>delta(id));
            if peak1>0 & fwhm1>0
                gamma1=(t/peak1).^alpha1.*exp(-(t-peak1)./beta1);
            else
                gamma1=min(abs(t-peak1))==abs(t-peak1);
            end
            if peak2>0 & fwhm2>0
                gamma2=(t/peak2).^alpha2.*exp(-(t-peak2)./beta2);
            else
                gamma2=min(abs(t-peak2))==abs(t-peak2);
            end
            hrf=gamma1-dip*gamma2;
        end
        H(:,id)=hrf/sum(hrf);
    end
    
    % Taylor coefs:
    origin=-round(Delta1/dt);
    HS0=[zeros(origin,2); HS(1:(numlags-origin),:)];
    WS=pinv(HS0)*H;
    X_cache.W(:,k,1:2)=WS';
    prcnt_var_taylor=sum(sum(H.*(HS0*WS)))/sum(sum(H.*H))*100;
    
    % svd:
    [U,SS,V]=svd(H,0);
    prcnt_var_spectral=(SS(1,1)^2+SS(2,2)^2)/sum(diag(SS).^2)*100;
    sumU=sum(U(:,1));
    US=U(:,1:2)/sumU;
    WS=V(:,1:2)*SS(1:2,1:2)*sumU;
    if delta*WS(:,2)<0
        US(:,2)=-US(:,2);
        WS(:,2)=-WS(:,2);
    end
    temp=conv2(response(:,k),US);
    eventmatrix(:,k,3:4)=temp((1:numtimes)-round(Delta1/dt),:);
    X_cache.W(:,k,3:4)=WS;
    X_cache.W(:,k,5)=delta';
    
    if ~all(WS(:,1)>0)
        fprintf(['Warning: use only for magnitudes, not delays \n first coef not positive for stimulus ' num2str(k)]);
    end
    cubic_coef=pinv([delta' delta'.^3])*(WS(:,2)./WS(:,1));
    if prod(cubic_coef)<0
        fprintf(['\nWarning: use only for magnitudes, not delays \n svd ratio not invertible for stimulus ' num2str(k)]);
    end
end

X_cache.X=zeros(n,numresponses,4,numslices);

for slice = 1:numslices
    subtime=ceil((frametimes+slicetimes(slice)-startime)/dt)+1;
    X_cache.X(:,:,:,slice)=eventmatrix(subtime,:,:);
end

X_cache.TR=(max(frametimes)-min(frametimes))/(length(frametimes)-1);

return
