function fwrite_EEG(filename,EEG,pstart,pstop )
%Write Segmented EEG file
%filename
%EEG structure
%EEG.data data 
%EEG.marker marker label 
%EEG.infoBV parameter of reading 
%EEG.ind_dur_ch indice duration channel marker 
%pstart in sample 
%pstop in sample  
[dirEEG,filEEG,extEEG]= fileparts(filename);   
if ~isdir(dirEEG)
    mkdir(dirEEG)
    disp(['Create dir: ',dirEEG])
end
outfileEEG_vmrk= fullfile(dirEEG,[filEEG '.vmrk']);
outfileEEG_vhdr = fullfile(dirEEG,[filEEG '.vhdr']);
    try
            %open a NIR data file and save the data            
            fid = fopen(filename,'w');
            %protection against change in vectorized of multiplex
            if ischar(EEG.infoBV.NumberOfChannels)
                Nch = str2num(EEG.infoBV.NumberOfChannels);
            else
                Nch = EEG.infoBV.NumberOfChannels;
            end
            if size(EEG.data,2)== Nch
                EEG.infoBV.DataOrientation = 'VECTORIZED';               
                d = EEG.data(pstart:pstop,:);
                EEG.infoBV.DataPoints  = size(d,1);
                d = d(:);
            elseif size(EEG.data,1)== Nch
                EEG.infoBV.DataOrientation = 'MULTIPLEXED';
                d = EEG.data(:,pstart:pstop);
                EEG.infoBV.DataPoints  = size(d,2);
                d = d(:);
            end
            fwrite(fid,d,'float32',0,'ieee-le');      
            fclose(fid);   
            if ~isfield(EEG.infoBV,'coor_r')
              EEG.infoBV.coor_r = zeros(EEG.infoBV.NumberOfChannels,1);end 
            if ~isfield(EEG.infoBV,'coor_theta')
              EEG.infoBV.coor_theta = zeros(EEG.infoBV.NumberOfChannels,1);end
            if ~isfield(EEG.infoBV,'coor_phi')
              EEG.infoBV.coor_phi = zeros(EEG.infoBV.NumberOfChannels,1);end
                
          
            fwrite_BVvhdr(outfileEEG_vhdr,EEG.infoBV);
            idseg = find(EEG.ind_dur_ch(:,1)>pstart&EEG.ind_dur_ch(:,1)<pstop);
            ind_dur_ch = EEG.ind_dur_ch(idseg,:);          
            ind_dur_ch(:,1) = EEG.ind_dur_ch(idseg,1)-pstart*ones(size(EEG.ind_dur_ch(idseg,1))); %remove pretime in the new file
            marker = EEG.marker(idseg,:);
            fwrite_BVvmrk(outfileEEG_vmrk,ind_dur_ch,marker,0);
    catch
        disp('Failed to write file'); 
    end
end

function nbMK=fwrite_BVvmrk(file,ind_dur_ch,type,nbMK)
% Write events in the vmrk binary file.
% Inputs : filepath, type of artefact, description of artefact and a N x 3
% matrix of the data points (indices), duration and channels for each of
% the N markers.
% ex.: write_BVvmrk(file,ind_dur_ch, type, mk_nb);
% nbMK = nombre de marker


if nbMK>1 %append just add new marker
    fid = fopen(file,'a');
else
    fid = fopen(file,'w');
    fprintf(fid,'%s\r\n','[Marker Infos]');
    fprintf(fid,'%s\r\n','; Each entry: Mk<Marker number>=<Type>,<Description>,<Position in data points>,');
    fprintf(fid,'%s\r\n','; <Size in data points>, <Channel number (0 = marker is related to all channels)>,');
    fprintf(fid,'%s\r\n','; <Date (YYYYMMDDhhmmssuuuuuu)>');
    fprintf(fid,'%s\r\n','; Fields are delimited by commas, some fields might be omited (empty).');
    fprintf(fid,'%s\r\n','; Commas in type or description text are coded as "\1".');
    nbMK=1;
end

if strcmp(type,'Time 0') 
    desc = ' ';
   fprintf(fid, '%s%s%s%s%s%s%s%s%s%s%s%s\r\n',...
            'Mk',num2str(nbMK),'=',type,',',...
            desc,',',int2str(ind_dur_ch(1)),',',...
            int2str(ind_dur_ch(2)),',',int2str(ind_dur_ch(3)));
        nbMK = nbMK + 1;
else
    for Midx=1:size(ind_dur_ch,1)   %For Step_detection
        desc = num2str(1); 
        %tmp = type{Midx,2}
        tmp = [type{Midx,1},',',type{Midx,2}];
%         fprintf(fid, '%s%s%s%s%s%s%s%s%s%s%s%s\r\n',...
%             'Mk',int2str(nbMK),'=',tmp,',',...
%             desc,',',int2str(ind_dur_ch(Midx,1)),',',...
%             int2str(ind_dur_ch(Midx,2)),',',int2str(ind_dur_ch(Midx,3)));
%             nbMK = nbMK+1;
        fprintf(fid, '%s%s%s%s%s%s%s%s%s%s\r\n',...
            'Mk',int2str(nbMK),'=',tmp,',',int2str(ind_dur_ch(Midx,1)),',',...
            int2str(ind_dur_ch(Midx,2)),',',int2str(ind_dur_ch(Midx,3)));
            nbMK = nbMK+1;
    end
end    
fclose(fid);
%SamplingInterval=51200
end

function fwrite_BVvhdr(file,infoBV)
fid = fopen(file,'w');   
fprintf(fid,'%s\r\n','Brain Vision Data Exchange Header File Version 1.0');
fprintf(fid,'%s\r\n','MATLAB FILE');
fprintf(fid,'\r\n%s\r\n','[Common Infos]');
fprintf(fid,'%s\r\n','Codepage=UTF-8');
fprintf(fid,'%s\r\n',['; File create from : ',file ]);

[pathstr, name, ext] = fileparts(file);
fprintf(fid,'%s\r',['DataFile=',[name, '.dat']]);
[pathstr, name, ext] = fileparts(file);

fprintf(fid,'%s\r\n',['MarkerFile=',[name, '.vmrk']]);
fprintf(fid,'%s%s\r\n','DataFormat=',infoBV.DataFormat);
fprintf(fid,'%s\r\n','; Data orientation: VECTORIZED=ch1,pt1, ch1,pt2..., MULTIPLEXED=ch1,pt1, ch2,pt1 ...');
fprintf(fid,'%s%s\r\n','DataOrientation=', infoBV.DataOrientation);
fprintf(fid,'%s%s\r\n','DataType=', infoBV.DataType);

fprintf(fid,'%s\r\n',['NumberOfChannels=',num2str(infoBV.NumberOfChannels)]);
fprintf(fid,'%s\r\n',['DataPoints=',num2str(infoBV.DataPoints)]);
fprintf(fid,'%s\r\n','; Sampling interval in microseconds if time domain (convert to Hertz:');
fprintf(fid,'%s\r\n','; 1000000 / SamplingInterval) or in Hertz if frequency domain:');
fprintf(fid,'%s\r\n',['SamplingInterval=',num2str(infoBV.SamplingInterval)]);
if isfield(infoBV,'SegmentDataPoints')
fprintf(fid,'%s\r\n','SegmentationType=MARKERBASED')
fprintf(fid,'%s\r\n',['SegmentDataPoints=',num2str(infoBV.SegmentDataPoints)])
end
fprintf(fid,'%s\r\n','[BINARY infos]');
fprintf(fid,'%s\r\n','BinaryFormat=IEEE_FLOAT_32');
fprintf(fid,'%s\r\n','[Channel Infos]');
fprintf(fid,'%s\r\n','; Each entry: Ch<Channel number>=<Name>,<Reference channel name>,');
fprintf(fid,'%s\r\n','; <Resolution in "Unit">,<Unit>, Future extensions...');
fprintf(fid,'%s\r\n','; Fields are delimited by commas, some fields might be omited (empty).');
fprintf(fid,'%s\r\n','; Commas in channel names are coded as "\1".');
 
for ich=1:numel(infoBV.name_ele) 
    fprintf(fid,'%s\r\n',['Ch',num2str(ich),'=',infoBV.name_ele{ich},',,,µV']);
end
fprintf(fid,'%s\r\n','[Coordinates]');
fprintf(fid,'%s\r\n','; Each entry: Ch<Channel number>=<Radius>,<Theta>,<Phi>');
for ich=1:numel(infoBV.name_ele)
    fprintf(fid,'%s%f,%f,%f\r\n',['Ch',num2str(ich),'='],infoBV.coor_r(ich),infoBV.coor_theta(ich),infoBV.coor_phi(ich));
end
fclose(fid);

end