function [ind_dur_ch_new] = detect_subintervals_dur(ind_dur_ch)
% Process the duration of each subinterval and detects subintervals of
% sufficient duration.
% Inputs :
% Outputs :

a = 1;

for marker = 1:size(ind_dur_ch,1)-1
    ch = ind_dur_ch(marker,3);
    if ind_dur_ch(marker+1,3) == ch
        dur = ind_dur_ch(marker+1,1) - (ind_dur_ch(marker,1)+ind_dur_ch(marker,2));
        ind_dur_ch_new(a,1) = ind_dur_ch(marker,1)+ind_dur_ch(marker,2);
        ind_dur_ch_new(a,2) = dur;
        ind_dur_ch_new(a,3) = ch;
        ind_dur_ch_new(a,4) = 2; %For long_sub (will be nullified)
        a = a+1;
    end
end

try
    ind_dur_ch_new = [ind_dur_ch; ind_dur_ch_new];
catch 
    ind_dur_ch_new = ind_dur_ch;
end