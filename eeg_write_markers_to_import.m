function eeg_write_markers_to_import(OutputFile,DataFile,Markers, SamplingRate)
%info required:
%DataFile
%New Segment and its date
%Other markers (optional)
%Markers are a column of cells, which is ordered along the field
%Position in data points

fid = fopen(OutputFile,'wt'); %write


%fprintf(fid, '%s%s%s%s%s%s%s', 'Sampling rate: ', num2str(SamplingRate), 'Hz, SamplingInterval: ', num2str(1/SamplingRate*1000), 'ms' )
fprintf(fid,'%s\n','Sampling rate: 320Hz, SamplingInterval: 3.125ms')
%fprintf(fid, '\n%s\n', 'Type, Description, Position, Lenght, Channel' )
fprintf(fid,'%s\n','Type, Description, Position, Length, Channel')

for Midx=1:size(Markers,1)
fprintf(fid,'%s%s%s%s%s%s%s%s%s',... %
Markers{Midx,1}.Type,',',...
Markers{Midx,1}.Description,',',int2str(Markers{Midx,1}.Position),',',...
int2str(Markers{Midx,1}.Size),',',Markers{Midx,1}.ChNumber);
%'Mk1=New Segment,,1,1,0,20100121160612567087');
if Midx<size(Markers,1)
    fprintf(fid,'\n')
end
end
fclose(fid);
