function out = nirs_run_E_manualtriglsl(job)
load(job.NIRSmat{1,1});
for ifile = 1:numel(NIRS.Dt.fir.pp(end).p)
    try
        fileevt = job.m_trigfilelsl{ifile};
        evtLSL = load(fileevt);
        evt = evtLSL(:,[3,2]);
        disp(['LSL trigger file ', fileevt, ' is read with ' num2str(size(evt,1)), ' event'])
    catch
        disp(['LSL trigger file ', fileevt, 'could not be read'])
    end
    if job.m_trigmode==1 %ajouter les évènement à aux5
        aux5 = NIRS.Dt.fir.aux5{ifile};
        aux5 = [aux5;evt];
        NIRS.Dt.fir.aux5{ifile}=aux5;
        clear aux5;
    else %remplacer aux5 par le trig manuel
        aux5 = [aux5;evt];
        NIRS.Dt.fir.aux5{ifile}=aux5;
        clear aux5;
    end
end
save(job.NIRSmat{1,1},'NIRS');
out.NIRSmat = job.NIRSmat;