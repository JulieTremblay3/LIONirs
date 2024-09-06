function d = fopen_NIR(varargin)
%Convention: .nir files written in multiplexed format (channels x time),
%as a long vector, and require specification of number of channels
%while .nirs files are written in matrix format (time x channels)
%output of fopen_NIR will be in multiplexed format
try
    location = varargin{1};
    if length(varargin) >= 2
        NC = varargin{2};
    end
    [dummy,dummy2,ext1] = fileparts(deblank(location));
    if strcmp(ext1,'.nirs')
        load(location,'-mat');
        d = d';
    else
        %open a NIR data file and return the data 
        fid = fopen(location,'r');
        d = fread(fid,'float32',0,'ieee-le');
        %data is channels in rows by time points in columns
        d = reshape(d,NC,[]);
        fclose(fid);             
    end
catch
    disp(['Failed to open file :',location ]);
end
