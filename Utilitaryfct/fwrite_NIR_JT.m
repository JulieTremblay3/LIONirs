function fwrite_NIR(varargin)
%Convention: .nir files written in multiplexed format (channels x time),
%as a long vector, and require specification of number of channels
%while .nirs files are written in matrix format (time x channels)
try
    location = varargin{1};
    d = varargin{2};
    try
        NIRS = varargin{3};
        HomerOn = 1;
        if numel(varargin)>3
            noise = varargin{4};
        else
            noise =[]
        end
    catch
        HomerOn=0;
    end
    [~, ~,ext1] = fileparts(location);
    if strcmp(ext1,'.nirs')
        d = d';
        if HomerOn
            SD = [];
            SD.Lambda = NIRS.Cf.dev.wl;
            SD.SrcPos = NIRS.Cf.H.S.r.o.mm.p';
            SD.DetPos = NIRS.Cf.H.D.r.o.mm.p';
            if numel(NIRS.Cf.H.C.id(2,:)) == size(d,2)
                ml = [NIRS.Cf.H.C.id(2,:)',NIRS.Cf.H.C.id(3,:)', ones(numel(NIRS.Cf.H.C.id(2,:)),1),NIRS.Cf.H.C.wl']; %srs,det,1,wav
                t = 1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:size(d,1)*1/NIRS.Cf.dev.fs;
                s = zeros(numel(t),1);
                aux = zeros(numel(t),6);
                prj_name = 's';
            end
                   
            if isempty(noise)
                save(location,'d','SD','ml','t','s','aux','-mat');
            else
                save(location,'d','SD','ml','t','s','aux','noise','-mat');
            end
        else
            save(location,'d','-mat');
        end
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