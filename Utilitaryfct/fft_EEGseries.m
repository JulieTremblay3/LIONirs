function [Cf, freq] = fft_EEGseries(data,SampFreq)

%Function for computing FFT for an EEG segment
%
% Input
% data   -> EEG segment: matrix No. of channels by No. of time points
% SampFreq -> Sampling Frequency in Hz
%
% Output
% Cf    -> Scaled Fourier coefficients: matrix No. of frequencies by No. of channels
% freq  -> Vector with frequency values. If freq(1)=0, Cf(1) is a real
%          value showing the DC level of the signal


%% Preliminaries
%check inputs?
[Nc,Nt] = size(data);
L = Nt; %2.^nextpow2(Nt);
data = [data zeros(Nc,L-Nt)];
deltaT = 1./SampFreq;    % time step
T=L.*deltaT;             % total time

%% Computing FFT
%data=detrend(data')';   % Removing linear trending.
Cf = fft(data')./L;      % Fast Fourier Transform for each channel (whole segment)

% without DC
%freq=(1:L/2)./T;         % Frequency grid
%Cf(1,:) = [];            % Removing DC, zero frequency
%Cf = 2*Cf(1:L/2,:);      % Keeping just one half of the frequency spectra and correcting for the energy (mult by 2)

% with DC
freq=(0:floor(L/2))./T;                  % Frequency grid
Cf = [Cf(1,:); 2*Cf(2:floor(L/2)+1,:)];    % Keeping just one half of the frequency spectra and correcting for the energy (mult by 2)

if 1==0
%% Test Example
clear Cf ffteegB
Nc=1;
Nt = 1000;
SampFreq = 100;
deltaT = 1./SampFreq;
t = (0:Nt-1)*deltaT; % time grid

%simulate signal
dc = .5; s0 = dc.*ones(Nc,Nt);    % DC level
s1 = zeros(Nc,Nt); A1 = 12; f1 = 10 ; Phi1 = 0; % volts, Hz and rad
s2 = zeros(Nc,Nt); A2 = 5; f2 = 23 ; Phi2 = pi/3; % volts, Hz and rad
s3 = zeros(Nc,Nt); A3 = 20; f3 = 4 ; Phi3 = pi; % volts, Hz and rad
range1 = [1:Nt];%[1:round(Nt/3)];%[1:Nt/2];%
range2 = [round(Nt/3)+1:2*round(Nt/3)];%[1:2*round(Nt/3)];%[Nt/2+1:Nt];%[1:Nt];%
range3 = [2*round(Nt/3)+1:Nt];%[1:Nt];%[1:2*round(Nt/3)];%[1:Nt/2];%
s1(:,range1) = A1.*sin(2*pi*f1.*t(range1) + Phi1); S = s0 + s1;
s2(:,range2) = A2.*sin(2*pi*f2.*t(range2) + Phi2); S = s0 + s2;
s3(:,range3) = A3.*sin(2*pi*f3.*t(range3) + Phi3); S = s0 + s3;
S = s0 + s1 + s2 + s3;
%correction for amplitude of spectrum
c1 = length(range1)./Nt; c2 = length(range2)./Nt; c3 = length(range3)./Nt;
AmpCorrected = [A1.*c1 A2.*c2 A3.*c3];
figure;
plot(t, S); xlabel('Time');title('Clean simulated Signal')
% S1=detrend(S')';
% p=polyfit(t,S,1);
% figure;
% plot(t, S1); xlabel('Time');title('Clean simulated Signal')

%add noise
SNR = 50;
noise = 0;%mvnrnd(zeros(Nc,1),max(S(:))./SNR.*eye(Nc),Nt)';
data = S + noise;
%data = S; 
% figure;
% plot(t, data); xlabel('Time');title('Noisy simulated Signal')

%%%using segments
Nwin=1; % use a divisor of Nt
data=reshape(data,Nc,Nt/Nwin,[]);

%computing FFT
tic
for i=1:Nwin
    [Cf(:,:,i), freq] = fft_EEGseries(data(:,:,i),SampFreq);
end
ted=toc;
Power1 = squeeze(Nt.*abs(Cf).^2);
for i=1:size(Cf,3); Sii(:,i) = Nt.*(squeeze(Cf(:,1,i)).*squeeze(Cf(:,1,i))'.')./(2*pi);end
Power2 = Sii;
figure;title('Power Estimates: blue=Nt*abs(x)^2; green=Nt*x*x''/2pi');for i=1:size(Power1,2);plot(freq,[Power1(:,i) Power2(:,i)]);pause;end
Ampest = squeeze(abs(Cf));
figure; title('Amplitude!');for i=1:size(Ampest,2);plot(freq,Ampest(:,i));pause;end

tic
[ffteegB, real_freqinterval] = fft_calcm(reshape(data,Nc,Nt), Nt/Nwin, deltaT,'none',1);
%[ffteegB, real_freqinterval] = fft_calc(reshape(data,Nc,Nt), Nt/Nwin, deltaT);
tb=toc;
%freq_arrayB = real_freqinterval(1):real_freqinterval(2):real_freqinterval(3);
freq_arrayB = 0:real_freqinterval(2):real_freqinterval(3);

%checking(optional)
for i=1:Nwin
    Cfp=Cf(:,:,i);
    ffteegBp=ffteegB(i,:);
    
    % Amplitude Spectrum
    figure
    plot(freq,abs(Cfp),'-o');  % Spectrum
    xlabel('Frequency');title(['Amplitude Spectrum. Segment' num2str(i)])
    hold on; plot(freq_arrayB,abs(ffteegBp),'-gs')
    plot([0,f1,f2,f3],[dc,A1 A2 A3],'*r')
    plot([0,f1,f2,f3],[dc AmpCorrected],'dk')
    % Phase Spectrum
    figure
    phase = angle(Cf);
    %phase = unwrap(angle(Cf));
    plot(freq,phase);  % Phases
    xlabel('Frequency');title('Phase Spectrum')
    hold on; plot(freq_arrayB,angle(ffteegBp),'-gs')
    plot([f1,f2,f3],[Phi1 Phi2 Phi3],'*r')
end

end
