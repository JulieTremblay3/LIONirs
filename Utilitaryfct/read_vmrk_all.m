function [label,ind_dur_ch] = read_vmrk_all(varargin)
% Read the .vmrk file for all marker type.
% Inputs : filepath and type of artefact
% Outputs : N x 3 matrix of the data points (indices), duration and
% channels for each of the N markers.

file = varargin{1};
fid = fopen(file,'r');
if fid==-1
    disp(['Verify file: ',file, ' could not be open.'])
end
label = [];
ind_dur_ch = [];
ind = [];
dur = [];
ch = [];
N=1;
i=4;
labellist = [];
%Get the markers infos

fline = fgetl(fid);

    
while ischar(fline)
    fline = fgetl(fid);
    if ischar(fline) && ~strcmp(fline,'');
        if strcmp(fline(1:2),'Mk')
            %Add all line in the table
            [tok, rem] = strtok(fline,'=');
            [label{N,1}, rem] = strtok(rem(2:end),',');
             sep=strfind(rem,',');      
              if (sep(2)-sep(1))>1 %Description n'est pas un espace vide
                [label{N,2}, rem]=strtok(rem,',');
              else
                label{N,2} = '' ;    
              end   
              [ind, rem]=strtok(rem,',');
              [dur, rem]=strtok(rem,',');
              [ch , rem]=strtok(rem,',');
              ind_dur_ch(N,:) = [str2num(ind), str2num(dur), str2num(ch)];
              N = N+1;
        end
    end
end


try
    fclose(fid);
catch
    disp('Error closing .vmrk file');
end