function write_vmrk_all(varargin)
% Write all events in the vmrk file.
% Inputs : filepath, ind_dur_ch, description of artefact and a N x 3
% matrix of the data points (indices), duration and channels for each of
% the N markers, label cell array type of artefact and description 
% ex.: write_vmrk(filepath,ind_dur_ch,label);

file = varargin{1};
ind_dur_ch = varargin{2};
label = varargin{3};
fid = fopen(file,'w+'); %write from begining
fprintf(fid, '%s\n\n','NIRS Data Marker File');
fprintf(fid, '%s\n','[Common Infos]');
fprintf(fid, '%s\n','Codepage=UTF-8');
fprintf(fid, '%s%s\n\n','DataFile=',[file(1:end-4),'nir']);
fprintf(fid, '%s\n','[Marker Infos]');
fprintf(fid, '%s\n','; Each entry: Mk<Marker number>=<Type>,<Description>,<Position in data points>,');
fprintf(fid, '%s\n',' ; <Size in data points>, <Channel number (0 = marker is related to all channels)>');
fprintf(fid, '%s\n',' ; Fields are delimited by commas, some fields might be omitted (empty).');
fprintf(fid, '%s\n',' ; Commas in type or description text are coded as "\1".');
for i = 1:size(label,1)
 fprintf(fid,'%s%d%s%d%s%d%s%d\n','Mk',i,['=',label{i,1},',',label{i,2},','],ind_dur_ch(i,1),',',ind_dur_ch(i,2),',',ind_dur_ch(i,3));
end
fclose(fid);




