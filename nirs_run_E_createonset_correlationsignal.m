function out = nirs_run_E_createonset_correlationsignal(job)
%Remode linear trend
%filename prefix
prefix = 'o'; %for "onset"
DelPreviousData  = job.DelPreviousData;
min_corr = str2num(job.e_correlationsignal_min);

for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
    %Load NIRS.mat information
    %     try
    NIRS = [];
    load(job.NIRSmat{filenb,1});
    [dir2,tmp,tmp] = fileparts(job.NIRSmat{filenb,1});
    
    %use last step of preprocessing
    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N;
    fs = NIRS.Cf.dev.fs;
    fprintf('%s\n',['Correlate the artefact: ', job.f_correlationsignal{1}, ' in the data']);
    
    %get the curve from the figure
    hfig = open(job.f_correlationsignal{1});
    haxes=get(hfig,'children');
    hlines=get(haxes,'children');
    for iline = 1:numel(hlines);
        artifact(:,iline)= get(hlines(iline),'ydata');
        xdata(:,iline)= get(hlines(iline),'xdata');
    end
    close(hfig)
    %check if data and artefac have the same sample rate, else just add
    %a warning but continue
    load(job.zonecorrelation{1},'-mat');
    channeluse = zone.plotLst{1};
    Fs_model =  1/(xdata(2,1)-xdata(1,1));
    if round(Fs_model)~=round(fs)
        disp('Warning: signal model and data do not have the same sample rate.')
    end
    
    for f=1:size(rDtp,1) %Loop over all files of a NIRS.mat
        %             try
        
        d = fopen_NIR(rDtp{f,1},NC);
        [dir1,fil1,~] = fileparts(rDtp{f});
        dfilt =d(:,:)';
        sizebloctocheck = size(dfilt,1);
        maxwindow = size(artifact,1)-1;
        if ~isempty(str2num(job.e_correlationsignal_LPF))%Lowpass
            LowPass_filter = str2num(job.e_correlationsignal_LPF);
            Wn = LowPass_filter*2/fs;
            filt_ord =4;
            [fb,fa]=butter(filt_ord,Wn,'low');
            idx = 1;
            dfilt = filtfilt(fb,fa, dfilt);
            disp(['Apply ',job.e_correlationsignal_LPF, ' Hz lowpass filter for the detection.'])
        end
        
        if ~isempty(str2num(job.e_correlationsignal_HPF))%Lowpass
            HighPass_filter = str2num(job.e_correlationsignal_HPF);
            Wn = HighPass_filter*2/fs;
            filt_ord =4;
            [fb,fa]=butter(filt_ord,Wn,'high');
            idx = 1;
            dfilt = filtfilt(fb,fa, dfilt);
            disp(['Apply ',job.e_correlationsignal_HPF, ' Hz highpass filter for the detection.'])
        end
        
 

        %on each channel the average artifact
        tic
        artmean= mean(artifact,2);
        nbmodel = size(artifact,1)-1;
        step = round(maxwindow/10);
        disp(['Apply correlation using step 10 % the size of the artifact ', num2str(step),' samples (', num2str(step*1/fs),' s)'] )
        %correlate the avg for each channel
        for idx=1:numel(channeluse)
            for isample = 1:step:(sizebloctocheck-maxwindow )
                [r(isample,idx)] = corr(dfilt((isample:isample+ maxwindow),channeluse(idx)), artmean);
                
            end
        end
        toc
        
        if isfield(NIRS.Cf.H.C,'label')
            labelall = [];
            for idx=1:numel(channeluse)
                labelall = [labelall, ' ' ,NIRS.Cf.H.C.label{(channeluse(idx))}];
            end
            disp( ['Correlate with' labelall,' at more then ', num2str(min_corr*100), '%'])
        end
        
        iwindow = 0;
        curvepeak =mean(r,2);                                            %mean coorelation over channel in the zone
        curvepeak = curvepeak.*double(curvepeak>min_corr) ;              %mask too small correlation to avoid to be use as local maximal
        
        %                 figure
        %                 time = 1/fs:1/fs:size(r,1)*1/fs;
        %                 plot(time,r,'r')
        %                  plot(r,'r')
        %                figure;plot(time,curvepeak)
        maxpeak = [1];
        
        while iwindow < (size(curvepeak,1)-round(maxwindow))
            try
                iwindow = iwindow + step;
                [val,id] = max(curvepeak(iwindow: iwindow+round(maxwindow)));
                if id==1 | val==0
                    1;
                else
                    idreal = id+iwindow;
                    if idreal~=maxpeak(end) | val==0
                        maxpeak =  [maxpeak,idreal];
                        iwindow =  iwindow+round(maxwindow);
                        curvepeak(iwindow: (iwindow+round(maxwindow)))=0;
                    end
                end
            catch
                
            end
            
        end
        
        maxpeak(1)=[];
        disp(['Found ', num2str(numel(maxpeak)),' events which correlate with the artifact and marked with trigger ', job.e_onset])
        %adjust the marke in the midel of the artifact time window
        maxpeak = maxpeak + round(maxwindow/2);
        disp(['Model artifact time duration: ', num2str(maxwindow*1/fs*1000), ' ms'])
        
        [dir1,fil1,ext1] = fileparts(rDtp{f});
        infilevmrk = fullfile(dir1,[fil1 '.vmrk']);
        try
            infilevhdr = fullfile(dir1,[fil1 '.vhdr']);
        catch
        end
        
        if ~exist(dir2,'dir'), mkdir(dir2); end
        outfile = fullfile(dir2,[prefix fil1 ext1]);
        outfilevmrk = fullfile(dir2,[prefix fil1 '.vmrk']);
        outfilevhdr = fullfile(dir2,[prefix fil1 '.vhdr']);
        
        
        fwrite_NIR(outfile,d);
        fprintf('%s\n',outfile);
        
        %write new .vmrk file
        try
            copyfile(infilevmrk,outfilevmrk);
        catch;end
        try
            try
                info = read_vhdr_brainvision((fullfile(dir1,[fil1,'.vhdr'])));
                ChannelLabels = info.label;
            catch
                ChannelLabels = ConvertmlIDsrs2label(NIRS);
            end
            SamplingInterval =floor(1000000/NIRS.Cf.dev.fs);
            nirs_boxy_write_vhdr(outfilevhdr,... %Output file
                outfile,... %DataFile
                outfilevmrk,... %MarkerFile,...
                'nirs_E_detrend',... %Function that created the header
                '',... %Channel Resolution
                '',... %Channel Units
                ChannelLabels,... %names given as a column of cells
                SamplingInterval,...
                size(d,2)); %SamplingInterval in microseconds
        catch
        end
        if DelPreviousData
            delete(rDtp{f,1});
            delete(infilevmrk)
            delete(infilevhdr)
            disp(['Delete previous .nir data file: ',rDtp{f,1}]);
        end
        %add outfile name to NIRS
        if f == 1
            NIRS.Dt.fir.pp(lst+1).pre = 'Set onset based on correlation';
            NIRS.Dt.fir.pp(lst+1).job = job;
        end
        NIRS.Dt.fir.pp(lst+1).p{f,1} = outfile;
        tmpaux5 = NIRS.Dt.fir.aux5{1,1};
        tmpaux5 = [tmpaux5;[ones(numel(maxpeak),1)*str2num(job.e_onset),maxpeak']];
        NIRS.Dt.fir.aux5{1,1} =  tmpaux5;
        
    end
    
    save(fullfile(dir2,'NIRS.mat'),'NIRS');
    job.NIRSmat{1} =fullfile(dir2,'NIRS.mat');
    
end
out.NIRSmat = job.NIRSmat;
function [num, den, z, p] = butter(n, Wn, varargin)
%BUTTER Butterworth digital and analog filter design.
%   [B,A] = BUTTER(N,Wn) designs an Nth order lowpass digital
%   Butterworth filter and returns the filter coefficients in length
%   N+1 vectors B (numerator) and A (denominator). The coefficients
%   are listed in descending powers of z. The cutoff frequency
%   Wn must be 0.0 < Wn < 1.0, with 1.0 corresponding to
%   half the sample rate.
%
%   If Wn is a two-element vector, Wn = [W1 W2], BUTTER returns an
%   order 2N bandpass filter with passband  W1 < W < W2.
%   [B,A] = BUTTER(N,Wn,'high') designs a highpass filter.
%   [B,A] = BUTTER(N,Wn,'low') designs a lowpass filter.
%   [B,A] = BUTTER(N,Wn,'stop') is a bandstop filter if Wn = [W1 W2].
%
%   When used with three left-hand arguments, as in
%   [Z,P,K] = BUTTER(...), the zeros and poles are returned in
%   length N column vectors Z and P, and the gain in scalar K.
%
%   When used with four left-hand arguments, as in
%   [A,B,C,D] = BUTTER(...), state-space matrices are returned.
%
%   BUTTER(N,Wn,'s'), BUTTER(N,Wn,'high','s') and BUTTER(N,Wn,'stop','s')
%   design analog Butterworth filters.  In this case, Wn is in [rad/s]
%   and it can be greater than 1.0.
%
%   See also BUTTORD, BESSELF, CHEBY1, CHEBY2, ELLIP, FREQZ, FILTER.

%   Author(s): J.N. Little, 1-14-87
%   	   J.N. Little, 1-14-88, revised
%   	   L. Shure, 4-29-88, revised
%   	   T. Krauss, 3-24-93, revised
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.8.4.4 $  $Date: 2007/12/14 15:03:53 $

%   References:
%     [1] T. W. Parks and C. S. Burrus, Digital Filter Design,
%         John Wiley & Sons, 1987, chapter 7, section 7.3.3.

[btype,analog,errStr] = iirchk(Wn,varargin{:});
if ~isempty(errStr), error(generatemsgid('SigErr'),errStr); end

if n>500
    error(generatemsgid('InvalidRange'),'Filter order too large.')
end

% step 1: get analog, pre-warped frequencies
if ~analog,
    fs = 2;
    u = 2*fs*tan(pi*Wn/fs);
else
    u = Wn;
end

Bw=[];
% step 2: convert to low-pass prototype estimate
if btype == 1	% lowpass
    Wn = u;
elseif btype == 2	% bandpass
    Bw = u(2) - u(1);
    Wn = sqrt(u(1)*u(2));	% center frequency
elseif btype == 3	% highpass
    Wn = u;
elseif btype == 4	% bandstop
    Bw = u(2) - u(1);
    Wn = sqrt(u(1)*u(2));	% center frequency
end

% step 3: Get N-th order Butterworth analog lowpass prototype
[z,p,k] = buttap(n);

% Transform to state-space
[a,b,c,d] = zp2ss(z,p,k);

% step 4: Transform to lowpass, bandpass, highpass, or bandstop of desired Wn
if btype == 1		% Lowpass
    [a,b,c,d] = lp2lp(a,b,c,d,Wn);
    
elseif btype == 2	% Bandpass
    [a,b,c,d] = lp2bp(a,b,c,d,Wn,Bw);
    
elseif btype == 3	% Highpass
    [a,b,c,d] = lp2hp(a,b,c,d,Wn);
    
elseif btype == 4	% Bandstop
    [a,b,c,d] = lp2bs(a,b,c,d,Wn,Bw);
end

% step 5: Use Bilinear transformation to find discrete equivalent:
if ~analog,
    [a,b,c,d] = bilinear(a,b,c,d,fs);
end

if nargout == 4
    num = a;
    den = b;
    z = c;
    p = d;
else	% nargout <= 3
    % Transform to zero-pole-gain and polynomial forms:
    if nargout == 3
        [z,p,k] = ss2zp(a,b,c,d,1); %#ok
        z = buttzeros(btype,n,Wn,analog);
        num = z;
        den = p;
        z = k;
    else % nargout <= 2
        den = poly(a);
        num = buttnum(btype,n,Wn,Bw,analog,den);
        % num = poly(a-b*c)+(d-1)*den;
        
    end
end

%---------------------------------
function b = buttnum(btype,n,Wn,Bw,analog,den)
% This internal function returns more exact numerator vectors
% for the num/den case.
% Wn input is two element band edge vector
if analog
    switch btype
        case 1  % lowpass
            b = [zeros(1,n) n^(-n)];
            b = real( b*polyval(den,-j*0)/polyval(b,-j*0) );
        case 2  % bandpass
            b = [zeros(1,n) Bw^n zeros(1,n)];
            b = real( b*polyval(den,-j*Wn)/polyval(b,-j*Wn) );
        case 3  % highpass
            b = [1 zeros(1,n)];
            b = real( b*den(1)/b(1) );
        case 4  % bandstop
            r = j*Wn*((-1).^(0:2*n-1)');
            b = poly(r);
            b = real( b*polyval(den,-j*0)/polyval(b,-j*0) );
    end
else
    Wn = 2*atan2(Wn,4);
    switch btype
        case 1  % lowpass
            r = -ones(n,1);
            w = 0;
        case 2  % bandpass
            r = [ones(n,1); -ones(n,1)];
            w = Wn;
        case 3  % highpass
            r = ones(n,1);
            w = pi;
        case 4  % bandstop
            r = exp(j*Wn*( (-1).^(0:2*n-1)' ));
            w = 0;
    end
    b = poly(r);
    % now normalize so |H(w)| == 1:
    kern = exp(-j*w*(0:length(b)-1));
    b = real(b*(kern*den(:))/(kern*b(:)));
end

function z = buttzeros(btype,n,Wn,analog)
% This internal function returns more exact zeros.
% Wn input is two element band edge vector
if analog
    % for lowpass and bandpass, don't include zeros at +Inf or -Inf
    switch btype
        case 1  % lowpass
            z = zeros(0,1);
        case 2  % bandpass
            z = zeros(n,1);
        case 3  % highpass
            z = zeros(n,1);
        case 4  % bandstop
            z = j*Wn*((-1).^(0:2*n-1)');
    end
else
    Wn = 2*atan2(Wn,4);
    switch btype
        case 1  % lowpass
            z = -ones(n,1);
        case 2  % bandpass
            z = [ones(n,1); -ones(n,1)];
        case 3  % highpass
            z = ones(n,1);
        case 4  % bandstop
            z = exp(j*Wn*( (-1).^(0:2*n-1)' ));
    end
end
function [btype,analog,errStr] = iirchk(Wn,varargin)
%IIRCHK  Parameter checking for BUTTER, CHEBY1, CHEBY2, and ELLIP.
%   [btype,analog,errStr] = iirchk(Wn,varargin) returns the
%   filter type btype (1=lowpass, 2=bandpss, 3=highpass, 4=bandstop)
%   and analog flag analog (0=digital, 1=analog) given the edge
%   frequency Wn (either a one or two element vector) and the
%   optional arguments in varargin.  The variable arguments are
%   either empty, a one element cell, or a two element cell.
%
%   errStr is empty if no errors are detected; otherwise it contains
%   the error message.  If errStr is not empty, btype and analog
%   are invalid.

%   Copyright 1988-2002 The MathWorks, Inc.
% $Revision: 1.7 $

errStr = '';

% Define defaults:
analog = 0; % 0=digital, 1=analog
btype = 1;  % 1=lowpass, 2=bandpss, 3=highpass, 4=bandstop

if length(Wn)==1
    btype = 1;
elseif length(Wn)==2
    btype = 2;
else
    errStr = 'Wn must be a one or two element vector.';
    return
end

if length(varargin)>2
    errStr = 'Too many input arguments.';
    return
end

% Interpret and strip off trailing 's' or 'z' argument:
if length(varargin)>0
    switch lower(varargin{end})
        case 's'
            analog = 1;
            varargin(end) = [];
        case 'z'
            analog = 0;
            varargin(end) = [];
        otherwise
            if length(varargin) > 1
                errStr = 'Analog flag must be either ''z'' or ''s''.';
                return
            end
    end
end

% Check for correct Wn limits
if ~analog
    if any(Wn<=0) | any(Wn>=1)
        errStr = 'The cutoff frequencies must be within the interval of (0,1).';
        return
    end
else
    if any(Wn<=0)
        errStr = 'The cutoff frequencies must be greater than zero.';
        return
    end
end

% At this point, varargin will either be empty, or contain a single
% band type flag.

if length(varargin)==1   % Interpret filter type argument:
    switch lower(varargin{1})
        case 'low'
            btype = 1;
        case 'bandpass'
            btype = 2;
        case 'high'
            btype = 3;
        case 'stop'
            btype = 4;
        otherwise
            if nargin == 2
                errStr = ['Option string must be one of ''high'', ''stop'',' ...
                    ' ''low'', ''bandpass'', ''z'' or ''s''.'];
            else  % nargin == 3
                errStr = ['Filter type must be one of ''high'', ''stop'',' ...
                    ' ''low'', or ''bandpass''.'];
            end
            return
    end
    switch btype
        case 1
            if length(Wn)~=1
                errStr = 'For the ''low'' filter option, Wn must have 1 element.';
                return
            end
        case 2
            if length(Wn)~=2
                errStr = 'For the ''bandpass'' filter option, Wn must have 2 elements.';
                return
            end
        case 3
            if length(Wn)~=1
                errStr = 'For the ''high'' filter option, Wn must have 1 element.';
                return
            end
        case 4
            if length(Wn)~=2
                errStr = 'For the ''stop'' filter option, Wn must have 2 elements.';
                return
            end
    end
end
function [z,p,k] = buttap(n)
%BUTTAP Butterworth analog lowpass filter prototype.
%   [Z,P,K] = BUTTAP(N) returns the zeros, poles, and gain
%   for an N-th order normalized prototype Butterworth analog
%   lowpass filter.  The resulting filter has N poles around
%   the unit circle in the left half plane, and no zeros.
%
%   See also BUTTER, CHEB1AP, CHEB2AP, ELLIPAP.

%   Author(s): J.N. Little and J.O. Smith, 1-14-87
%   	   L. Shure, 1-13-88, revised
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.6 $  $Date: 2002/04/15 01:09:52 $

% Poles are on the unit circle in the left-half plane.
z = [];
p = exp(i*(pi*(1:2:n-1)/(2*n) + pi/2));
p = [p; conj(p)];
p = p(:);
if rem(n,2)==1   % n is odd
    p = [p; -1];
end
k = real(prod(-p));

function y = filtfilt(b,a,x)
%FILTFILT Zero-phase forward and reverse digital filtering.
%   Y = FILTFILT(B, A, X) filters the data in vector X with the filter described
%   by vectors A and B to create the filtered data Y.  The filter is described
%   by the difference equation:
%
%     y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                      - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%
%
%   After filtering in the forward direction, the filtered sequence is then
%   reversed and run back through the filter; Y is the time reverse of the
%   output of the second filtering operation.  The result has precisely zero
%   phase distortion and magnitude modified by the square of the filter's
%   magnitude response.  Care is taken to minimize startup and ending
%   transients by matching initial conditions.
%
%   The length of the input x must be more than three times
%   the filter order, defined as max(length(b)-1,length(a)-1).
%
%   Note that FILTFILT should not be used with differentiator and Hilbert FIR
%   filters, since the operation of these filters depends heavily on their
%   phase response.
%
%   See also FILTER.

%   References:
%     [1] Sanjit K. Mitra, Digital Signal Processing, 2nd ed., McGraw-Hill, 2001
%     [2] Fredrik Gustafsson, Determining the initial states in forward-backward
%         filtering, IEEE Transactions on Signal Processing, pp. 988--992, April 1996,
%         Volume 44, Issue 4

%   Author(s): L. Shure, 5-17-88
%   revised by T. Krauss, 1-21-94
%   Initial Conditions: Fredrik Gustafsson
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.7.4.3 $  $Date: 2007/12/14 15:04:33 $

error(nargchk(3,3,nargin,'struct'))
if (isempty(b) || isempty(a) || isempty(x))
    y = [];
    return
end

[m,n] = size(x);
if (n>1) && (m>1)
    y = x;
    for i=1:n  % loop over columns
        y(:,i) = filtfilt(b,a,x(:,i));
    end
    return
    % error(generatemsgid('InvalidDimensions'),'Only works for vector input.')
end
if m==1
    x = x(:);   % convert row to column
end
len = size(x,1);   % length of input
b = b(:).';
a = a(:).';
nb = length(b);
na = length(a);
nfilt = max(nb,na);

nfact = 3*(nfilt-1);  % length of edge transients

if (len<=nfact),    % input data too short!
    error(generatemsgid('InvalidDimensions'),'Data must have length more than 3 times filter order.');
end

% set up filter's initial conditions to remove dc offset problems at the
% beginning and end of the sequence
if nb < nfilt, b(nfilt)=0; end   % zero-pad if necessary
if na < nfilt, a(nfilt)=0; end
% use sparse matrix to solve system of linear equations for initial conditions
% zi are the steady-state states of the filter b(z)/a(z) in the state-space
% implementation of the 'filter' command.
rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2];
cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1];
data = [1+a(2) a(3:nfilt) ones(1,nfilt-2)  -ones(1,nfilt-2)];
sp = sparse(rows,cols,data);
zi = sp \ ( b(2:nfilt).' - a(2:nfilt).'*b(1) );
% non-sparse:
% zi = ( eye(nfilt-1) - [-a(2:nfilt).' [eye(nfilt-2); zeros(1,nfilt-2)]] ) \ ...
%      ( b(2:nfilt).' - a(2:nfilt).'*b(1) );

% Extrapolate beginning and end of data sequence using a "reflection
% method".  Slopes of original and extrapolated sequences match at
% the end points.
% This reduces end effects.
y = [2*x(1)-x((nfact+1):-1:2);x;2*x(len)-x((len-1):-1:len-nfact)];

% filter, reverse data, filter again, and reverse data again
y = filter(b,a,y,zi*y(1));
y = y(length(y):-1:1);
y = filter(b,a,y,zi*y(1));
y = y(length(y):-1:1);

% remove extrapolated pieces of y
y([1:nfact len+nfact+(1:nfact)]) = [];

if m == 1
    y = y.';   % convert back to row if necessary
end
function [at,bt,ct,dt] = lp2lp(a,b,c,d,wo)
%LP2LP Lowpass to lowpass analog filter transformation.
%   [NUMT,DENT] = LP2LP(NUM,DEN,Wo) transforms the lowpass filter
%   prototype NUM(s)/DEN(s) with unity cutoff frequency of 1 rad/sec
%   to a lowpass filter with cutoff frequency Wo (rad/sec).
%   [AT,BT,CT,DT] = LP2LP(A,B,C,D,Wo) does the same when the
%   filter is described in state-space form.
%
%   See also BILINEAR, IMPINVAR, LP2BP, LP2BS and LP2HP

%   Author(s): J.N. Little and G.F. Franklin, 8-4-87
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.8.4.2 $  $Date: 2004/12/26 22:16:17 $

if nargin == 3		% Transfer function case
    % handle column vector inputs: convert to rows
    if size(a,2) == 1
        a = a(:).';
    end
    if size(b,2) == 1
        b = b(:).';
    end
    % Transform to state-space
    wo = c;
    [a,b,c,d] = tf2ss(a,b);
end

error(abcdchk(a,b,c,d));

% Transform lowpass to lowpass
at = wo*a;
bt = wo*b;
ct = c;
dt = d;

if nargin == 3		% Transfer function case
    % Transform back to transfer function
    [z,k] = tzero(at,bt,ct,dt);
    num = k * poly(z);
    den = poly(at);
    at = num;
    bt = den;
end
function [zd, pd, kd, dd] = bilinear(z, p, k, fs, fp, fp1)
%BILINEAR Bilinear transformation with optional frequency prewarping.
%   [Zd,Pd,Kd] = BILINEAR(Z,P,K,Fs) converts the s-domain transfer
%   function specified by Z, P, and K to a z-transform discrete
%   equivalent obtained from the bilinear transformation:
%
%      H(z) = H(s) |
%                  | s = 2*Fs*(z-1)/(z+1)
%
%   where column vectors Z and P specify the zeros and poles, scalar
%   K specifies the gain, and Fs is the sample frequency in Hz.
%
%   [NUMd,DENd] = BILINEAR(NUM,DEN,Fs), where NUM and DEN are
%   row vectors containing numerator and denominator transfer
%   function coefficients, NUM(s)/DEN(s), in descending powers of
%   s, transforms to z-transform coefficients NUMd(z)/DENd(z).
%
%   [Ad,Bd,Cd,Dd] = BILINEAR(A,B,C,D,Fs) is a state-space version.
%
%   Each of the above three forms of BILINEAR accepts an optional
%   additional input argument that specifies prewarping.
%
%   For example, [Zd,Pd,Kd] = BILINEAR(Z,P,K,Fs,Fp) applies prewarping
%   before the bilinear transformation so that the frequency responses
%   before and after mapping match exactly at frequency point Fp
%   (match point Fp is specified in Hz).
%
%   See also IMPINVAR.

%   Author(s): J.N. Little, 4-28-87
%   	   J.N. Little, 5-5-87, revised
%   Copyright 1988-2006 The MathWorks, Inc.
%   $Revision: 1.7.4.4 $  $Date: 2007/12/14 15:03:48 $

%   Gene Franklin, Stanford Univ., motivated the state-space
%   approach to the bilinear transformation.

[mn,nn] = size(z);
[md,nd] = size(p);

if (nd == 1 && nn < 2) && nargout ~= 4	% In zero-pole-gain form
    if mn > md
        error(generatemsgid('InvalidRange'),'Numerator cannot be higher order than denominator.')
    end
    if nargin == 5		% Prewarp
        fp = 2*pi*fp;
        fs = fp/tan(fp/fs/2);
    else
        fs = 2*fs;
    end
    z = z(isfinite(z));	 % Strip infinities from zeros
    pd = (1+p/fs)./(1-p/fs); % Do bilinear transformation
    zd = (1+z/fs)./(1-z/fs);
    % real(kd) or just kd?
    kd = (k*prod(fs-z)./prod(fs-p));
    zd = [zd;-ones(length(pd)-length(zd),1)];  % Add extra zeros at -1
    
elseif (md == 1 && mn == 1) || nargout == 4 %
    if nargout == 4		% State-space case
        a = z; b = p; c = k; d = fs; fs = fp;
        error(abcdchk(a,b,c,d));
        if nargin == 6			% Prewarp
            fp = fp1;		% Decode arguments
            fp = 2*pi*fp;
            fs = fp/tan(fp/fs/2)/2;
        end
    else			% Transfer function case
        if nn > nd
            error(generatemsgid('InvalidRange'),'Numerator cannot be higher order than denominator.')
        end
        num = z; den = p;		% Decode arguments
        if nargin == 4			% Prewarp
            fp = fs; fs = k;	% Decode arguments
            fp = 2*pi*fp;
            fs = fp/tan(fp/fs/2)/2;
        else
            fs = k;			% Decode arguments
        end
        % Put num(s)/den(s) in state-space canonical form.
        [a,b,c,d] = tf2ss(num,den);
    end
    % Now do state-space version of bilinear transformation:
    t = 1/fs;
    r = sqrt(t);
    t1 = eye(size(a)) + a*t/2;
    t2 = eye(size(a)) - a*t/2;
    ad = t2\t1;
    bd = t/r*(t2\b);
    cd = r*c/t2;
    dd = c/t2*b*t/2 + d;
    if nargout == 4
        zd = ad; pd = bd; kd = cd;
    else
        % Convert back to transfer function form:
        p = poly(ad);
        zd = poly(ad-bd*cd)+(dd-1)*p;
        pd = p;
    end
else
    error(generatemsgid('SignalErr'),'First two arguments must have the same orientation.')
end
function [at,bt,ct,dt] = lp2hp(a,b,c,d,wo)
%LP2HP Lowpass to highpass analog filter transformation.
%   [NUMT,DENT] = LP2HP(NUM,DEN,Wo) transforms the lowpass filter
%   prototype NUM(s)/DEN(s) with unity cutoff frequency to a
%   highpass filter with cutoff frequency Wo.
%   [AT,BT,CT,DT] = LP2HP(A,B,C,D,Wo) does the same when the
%   filter is described in state-space form.
%
%   See also BILINEAR, IMPINVAR, LP2BP, LP2BS and LP2LP

%   Author(s): J.N. Little and G.F. Franklin, 8-4-87
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.7.4.2 $  $Date: 2004/12/26 22:16:16 $

if nargin == 3		% Transfer function case
    % handle column vector inputs: convert to rows
    if size(a,2) == 1
        a = a(:).';
    end
    if size(b,2) == 1
        b = b(:).';
    end
    % Transform to state-space
    wo = c;
    [a,b,c,d] = tf2ss(a,b);
end

error(abcdchk(a,b,c,d));

% Transform lowpass to highpass
at =  wo*inv(a);
bt = -wo*(a\b);
ct = c/a;
dt = d - c/a*b;

if nargin == 3		% Transfer function case
    % Transform back to transfer function
    [z,k] = tzero(at,bt,ct,dt);
    num = k * poly(z);
    den = poly(at);
    at = num;
    bt = den;
end