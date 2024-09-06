function fwrite_NIR(varargin)
%Convention: .nir files written in multiplexed format (channels x time),
%as a long vector, and require specification of number of channels
%while .nirs files are written in matrix format (time x channels)
try
    location = varargin{1};
    d = varargin{2};

    [dummy,dummy2,ext1] = fileparts(deblank(location));
    if strcmp(ext1,'.nirs')
        d = d';
        save(location,'d','-mat');
    else
        %open a NIR data file and save the data 
        fid = fopen(location,'w');
        d = d(:);
        fwrite(fid,d,'float32',0,'ieee-le');
        fclose(fid);             
    end
catch
    disp('Failed to write file'); 
end