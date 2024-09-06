function TFR = morletcwtEd(S, freqVec, Fs, width);

%%%Help or morletcwt
% INPUT
%    S         Signal: Array which first dimension is time in seconds.
%    freqVec   Frequency grid in Hz
%    Fs        sampling rate in Hertz
%    width     Wavelet parameter (see Tallon-Baudry, 1997)
% OUTPUT
%    TFR       Wavelets coefficients matrix
%    
% Authors: Fumi, Bornot, Eduardo

% dimensions of signal matrix S
size_S = size(S);
% number of time instants (latencies)
nlat = size_S(1);
% number of frequencies
nfreq = length(freqVec);

S = reshape(S, nlat, prod(size_S(2:end)));
TFR = zeros([nlat nfreq size(S, 2)]); 

S = detrend(S);
dt = 1/Fs;

for j = (1:nfreq)
%     disp(nfreq-j)
    
    f = freqVec(j);
    %sf = f/width;
    st = width/(2*pi*f);
    A = 1/sqrt(st*sqrt(2*pi)); % Eduardo
    
%     t=-3.5*st:dt:3.5*st;
    t=(-(nlat./2)+1)*dt:dt:(nlat./2)*dt;
    m = A*exp(-t.^2/(2*st^2)).*exp(i*2*pi*f.*t);
    y = conv2(S, m');
%     y = (2*y/Fs).^2;
    y = y(ceil(length(m)/2):size(y, 1)-floor(length(m)/2), :);
    %y = y(floor(length(m)/2)+1:size(y, 1)-floor(length(m)/2), :);
    TFR(:, j, :) = y; 
end

TFR = reshape(TFR, [nlat nfreq size_S(2:end)]);

return;