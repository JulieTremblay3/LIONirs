%___________________________________________________________________
% Copyright (C) 2019 LION Lab, Centre de recherche CHU Sainte-Justine 
% www.lionlab.umontreal.ca
%___________________________________________________________________
function nirs_boxy_write_markers(OutputFile,DataFile,Markers)
%info required:
%DataFile
%New Segment and its date
%Other markers (optional)
%Markers are a column of cells, which is ordered along the field
%Position in data points

fid = fopen(OutputFile,'wt'); %write
%first 6 lines
fprintf(fid,'%s\n\n%s\n%s\n%s%s\n\n',... %     
'NIRS Data Marker File',...
'[Common Infos]','Codepage=UTF-8',...
'DataFile=',DataFile);
%markers
fprintf(fid,'%s\n%s\n%s\n%s\n%s\n',... % 
'[Marker Infos]',...
'; Each entry: Mk<Marker number>=<Type>,<Description>,<Position in data points>,',...
'; <Size in data points>, <Channel number (0 = marker is related to all channels)>',...
'; Fields are delimited by commas, some fields might be omitted (empty).',...
'; Commas in type or description text are coded as "\1".');

for Midx=1:size(Markers,1)
fprintf(fid,'%s%s%s%s%s%s%s%s%s%s%s%s\n',... %
'Mk',int2str(Midx),'=',Markers{Midx,1}.Type,',',...
Markers{Midx,1}.Description,',',int2str(Markers{Midx,1}.Position),',',...
int2str(Markers{Midx,1}.Size),',',int2str(Markers{Midx,1}.ChNumber));
%'Mk1=New Segment,,1,1,0,20100121160612567087');
end
fclose(fid);
end