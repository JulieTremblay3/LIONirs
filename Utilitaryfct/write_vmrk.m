function write_vmrk(varargin)
% Write events in the vmrk ascii file.
% Inputs : filepath, type of artefact, description of artefact and a N x 3
% matrix of the data points (indices), duration and channels for each of
% the N markers.
% ex.: write_vmrk(filepath,type,description,ind_dur_ch);
% append in the current file

file = varargin{1};
type = varargin{2};
desc = varargin{3};
ind_dur_ch = varargin{4};

fid = fopen(file,'r+');
mk_nb = num2str(1);

%Get the last marker number (passer par la fin?...)
fline = fgetl(fid);
while ischar(fline)
    fline = fgetl(fid);
    if ischar(fline) && ~strcmp(fline,'')
        if strcmp(fline(1:2),'Mk')
            for i=4:7
                if strcmp(fline(i),'=')
                    mk_nb = fline(3:i-1);
                    break
                end 
            end
       end
    end
end

if size(desc,1) == 1
    for Midx=1:size(ind_dur_ch,1)   %For Step_detection
        fprintf(fid, '%s%s%s%s%s%s%s%s%s%s%s%s\n',...
            'Mk',int2str(str2num(mk_nb)+Midx),'=',type,',',...
            desc,',',int2str(ind_dur_ch(Midx,1)),',',...
            int2str(ind_dur_ch(Midx,2)),',',int2str(ind_dur_ch(Midx,3)));
    end
else
    for Midx=1:size(ind_dur_ch,1) %For ReadBOXY
    fprintf(fid, '%s%s%s%s%s%s%s%s%s%s%s%s\n',...
        'Mk',int2str(str2num(mk_nb)+Midx),'=',type,',',...
        desc{Midx},',',int2str(ind_dur_ch(Midx,1)),',',...
        int2str(ind_dur_ch(Midx,2)),',',int2str(ind_dur_ch(Midx,3)));
    end
end
    

try 
    fclose(fid);
catch
    disp('Error closing .vmrk file');
end


