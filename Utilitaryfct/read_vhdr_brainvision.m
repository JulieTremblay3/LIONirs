function info = read_vhdr_brainvision(file)
fid = fopen(file);
ilabel = 1;
while 1
    fline = fgetl(fid);
    if ~ischar(fline),   break,   end
    [rem,tok]=strtok(fline,'=');
    if numel(rem)<2
    elseif  strcmp(rem,'NumberOfChannels')
        info.NumberOfChannels = str2num(tok(2:end));
    elseif strcmp(rem,'DataType')
        info.DataType = tok(2:end);
    elseif strcmp(rem,'DataFormat')
        info.dataformat = (tok(2:end));
    elseif strcmp(rem,'DataOrientation')
        info.DataOrientation = (tok(2:end)); 
    elseif strcmp(rem,'BinaryFormat')
        info.BinaryFormat = tok(2:end);
    elseif strcmp(rem,'DataPoints')
        info.DataPoints = str2num(tok(2:end));
    elseif strcmp(rem,'SamplingInterval')
        info.SamplingInterval=str2num(tok(2:end));
    elseif strcmp(rem,'Layers')
        info.Layers=str2num(tok(2:end));
    elseif strcmp(rem,'SegmentDataPoints')
        info.SegmentDataPoints = str2num(tok(2:end));
    elseif strcmp(rem,'SegmentDataPointsPre')
        info.SegmentDataPointsPre = str2num(tok(2:end));  
    elseif strcmp(rem,'SegmentDataPointsPost')
        info.SegmentDataPointsPost = str2num(tok(2:end));  
    elseif strcmp(rem(1:2),'Ch') %& strcmp(tok(1:5),'=SEEG'
        [rem,tmp] =strtok(tok(2:end),',');
        info.label{ilabel,1} = rem; %(7:12);        
        ilabel = ilabel+1;
    end
end
fclose(fid);
