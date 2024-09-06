function d = fopen_NIR_removeNAN(varargin)
%Convention: .nir files written in multiplexed format (channels x time),
%as a long vector, and require specification of number of channels
%while .nirs files are written in matrix format (time x channels)
%output of fopen_NIR will be in multiplexed format
%This version will look for NaN segment and they will completly remove
%those it for the correlation purpose the resting part are dc ajust for
%continuity. 
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
        
% %remove old bad steps
% list_bad_step =  [];
% for i=1:size(label,1)
%     if strmatch(label{i,1},'bad_step')
%         list_bad_step = [list_bad_step; i];
%     end
% end
        %data is channels in rows by time points in columns
     %   [label,ind_dur_ch]=read_vmrk_all(last_file_vmrk);
        d = reshape(d,NC,[]);
        fclose(fid);           
        [label,ind_dur_ch] = read_vmrk_all([location(1:end-3),'vmrk']);
        list_bad_step = [];
        for i=1:size(label,1)
            if strmatch(label{i,1},'bad_step')
                list_bad_step = [list_bad_step; i];
            end
        end
        matok = zeros(size(d,2),1);
        for i = 1:numel(list_bad_step)
            tstart = ind_dur_ch(list_bad_step(i),1);
            tstop =  ind_dur_ch(list_bad_step(i),1)+ind_dur_ch(list_bad_step(i),2);
            matok(tstart:tstop,1)=1;
        end
%         figure
%         plot(matok)
        
        [ind_dur_ch2] = mat2d2ind_dur_ch(matok');
        dini = d;
        for i =1:size(ind_dur_ch2,1) %ajust
            pstart = ind_dur_ch2(i,1);
            pstop = ind_dur_ch2(i,1)+ind_dur_ch2(i,2);
            pend = size(d,2); 
            d(:,pstop:end) = d(:,pstop:end) - (d(:,pstop)-d(:,pstart))*ones(1,size(d(:,pstop:end),2));
        
        end
%         figure
%         hold on
%         plot(d(1,:),'r','displayname','corr')
%         plot(dini(1,:),'b','displayname','ini')
%         plot(matok)
%         idbad = find(matok);
%         dnew = d;
%         dnew(:,idbad ) = [];
%         plot(dnew(1,:),'g','displayname','corr remove bad')
        
    end
catch
    disp('Failed to open file');
end
