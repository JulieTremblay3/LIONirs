function [ind_dur_ch] = read_vmrk_find(varargin)
% Read the .vmrk file to find the specified type of marker.
% Inputs : filepath and type of artefact
% Outputs : N x 3 matrix of the data points (indices), duration and 
% channels for each of the N markers.

file = varargin{1};
type = varargin{2};

if ~iscellstr(type) %In order to always work with a cell-array
    temp_type = type;
    type =  [];
    type{1} = temp_type;
end

fid = fopen(file,'r');

ind_dur_ch = [];
ind = [];
dur = [];
ch = [];
N=1;
i=4;

%Get the markers infos
fline = fgetl(fid);
try
    while ischar(fline)
        fline = fgetl(fid);
        if ischar(fline) && ~strcmp(fline,'')
            if strcmp(fline(1:2),'Mk')
                %Find if this line has one of the marker types
                mrk_ind = [];
                for b = 1:size(type,1) %Loop over all specified types
                    mrk_ind = [mrk_ind, strfind(fline,char(type{b}))]; %In the case of a cell array type
                end
                [token,remain]  = strtok(type{b},','); %check if coma in the label, take token to start identification after first comma              

                if ~isempty(mrk_ind)
                    i = mrk_ind+length(token)+1; % 1 indice after first comma
                    while ~strcmp(fline(i),',') %To skip desc
                        i = i+1;
                    end
                    i = i+1;
                    while ~strcmp(fline(i),',') %Get the data point
                        ind = [ind, fline(i)];
                        i = i+1;
                    end
                    i = i+1;
                    while ~strcmp(fline(i),',') %Get the duration
                        dur = [dur, fline(i)];
                        i = i+1;
                    end
                    i = i+1;
                    for j=i:length(fline) %Get the channel
                        ch = [ch, fline(j)];
                    end
                    ind_dur_ch(N,:) = [str2num(ind), str2num(dur), str2num(ch)];
                    N = N+1;
                    ind = [];
                    dur = [];
                    ch = [];
                    i=4;
                end
            end
        end
    end
catch
    disp('Failed to get markers info.');
end

try 
    fclose(fid);
catch
    disp('Error closing .vmrk file');
end