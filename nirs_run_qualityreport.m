function out = nirs_run_qualityreport(job)
%report pourcentage of good channel
idrow = 1
xlsout{idrow,1}='NIRS.mat';
xlsout{idrow,2}='File'
xlsout{idrow,3}='Duration (s)';
xlsout{idrow,4}='Percentage rejected'
xlsout{idrow,5}='Percentage corrected'
xlsout{idrow,6}='Maximal duration corrected'

for filenb=1:size(job.NIRSmat,1) %JT
    idrow = idrow + 1;
    NIRS = [];
    load(job.NIRSmat{filenb,1});

    lst = length(NIRS.Dt.fir.pp);
    rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
    NC = NIRS.Cf.H.C.N;
    fs = NIRS.Cf.dev.fs;
    fprintf('%s\n','File processed');
    for f=1:size(rDtp,1) %Loop over all files of a NIRS.mat   
         d = fopen_NIR(rDtp{f,1},NC);
         [dir1,fil1,~] = fileparts(rDtp{f});
         vmrk_path = fullfile(dir1,[fil1 '.vmrk']);
         [ind_dur_ch] = read_vmrk_find(vmrk_path,'bad_step');
        %load noise
        noise = logical(zeros(size(d')));
        mrk_type_arr = cellstr('bad_step');
     mrks = [];
    ind = [];
    dur = [];
    [ind_dur_ch] = read_vmrk_find(vmrk_path,mrk_type_arr);
    if ~isempty(ind_dur_ch)
        maxpoint  = ind_dur_ch(:,1)+ind_dur_ch(:,2);
        badind = find(maxpoint>size(noise,1));
        if ~isempty(badind)
            disp(['Warning file ' vmrk_path ' marker : ' num2str(badind') ' are out of range'])
            ind_dur_ch(badind,2)=size(noise,1)- ind_dur_ch(badind,1);
        end
        for Idx = 1:size(noise,2)
            mrks = find(ind_dur_ch(:,3)==Idx);
            ind = ind_dur_ch(mrks,1);
            indf = ind + ind_dur_ch(mrks,2) - 1;
            if ~isempty(ind)
                try
                    for i = 1:numel(ind)
                        noise(ind(i):indf(i),Idx) = 1;
                    end
                catch
                    msgbox('Noise reading problem')
                end
            end
        end
    end
  
     xlsout{idrow,1} = job.NIRSmat{filenb,1} ;
     [pathstr,name,ext] =  fileparts(job.NIRSmat{filenb,1});
     xlsout{idrow,2} = fil1;
     xlsout{idrow,3} = num2str( 1/fs*size(d,2));
     xlsout{idrow,4} = num2str( sum(noise(:))/numel(noise(:))*100);
     load(fullfile(pathstr,'CorrectionApply.mat'));
     for icorr = 1:numel
        
  
        if PARCORR(icorr).file ==f;
            if strcmp(PARCORR(icorr).type,'PCA') | strcmp(PARCORR(icorr).type,'PARAFAC')
                indtsizemax
                
            end
        end
        
    end
end
xlswrite('test.xls',xlsout)
out.NIRSmat = job.NIRSmat;