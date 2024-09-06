function SD = nirs_boxy_extend_channel_names_lengths(SD)
%ChannelLabels and ChannelLengths (in centimeters
n_m=size(SD.ChnNames,1);
SD.ChannelLabels{2*n_m,1}='';
SD.ChannelLengths{2*n_m,1}='';
for Cidx = 1:n_m
    %For wavelength 830 nm
    SD.ChannelLabels{Cidx,1} = SD.ChnNames{Cidx,1};
    SD.ChannelLengths{Cidx,1} = SD.ChnDist(Cidx,1);
    %For wavelength 690 nm
    SD.ChannelLabels{Cidx+n_m,1} = SD.ChnNames{Cidx,1};
    SD.ChannelLengths{Cidx+n_m,1} = SD.ChnDist(Cidx,1);
end