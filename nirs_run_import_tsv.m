function out = nirs_run_import_tsv(job)
load(job.NIRSmat{1,1});
fileevt = job.f_filetsv{1};
fs = NIRS.Cf.dev.fs; %sample frequency nirs device;
for ifile = 1:numel(NIRS.Dt.fir.pp(end).p)
    try
        if isfield('NIRS.Dt.fir','comments')
          comments = NIRS.Dt.fir.comments{1};    
        else
            comments = [];
        end
        T = readtable(fileevt,'FileType', 'text') ; 
        for i=1:size(T,1)
            sample = round(T.start(i)/1000*fs);
            duration = round( (T.xEnd(i)-T.start(i))/1000*fs);
            comments = [comments; [T.text(i), {sample}, {duration}]];    
        end
    catch
        disp(['TSV file ', fileevt, 'could not be read'])
    end
    NIRS.Dt.fir.comments{1} =  comments;
    % if job.m_trigmode==1 %ajouter les évènement à aux5
    %     aux5 = NIRS.Dt.fir.aux5{ifile};
    %     aux5 = [aux5;evt];
    %     NIRS.Dt.fir.aux5{ifile}=aux5;
    %     clear aux5;
    % else %remplacer aux5 par le trig manuel
    %     aux5 = [];
    %     aux5 = [aux5;evt];
    %     NIRS.Dt.fir.aux5{ifile}=aux5;
    %     clear aux5;
    % end
end
save(job.NIRSmat{1,1},'NIRS');
out.NIRSmat = job.NIRSmat;