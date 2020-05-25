function out = nirs_run_E_manualtrig(job)

load(job.NIRSmat{1,1});

for itrigfile=1:numel(job.m_trigfile)
    run(job.m_trigfile{itrigfile});
    totaltrig  = 0;
    for ifile = 1:numel(NIRS.Dt.fir.pp(end).p)      
        if job.m_trigmode==1 %ajouter les évènement à aux5
            trigfile = [];
            for itrigifile = 1:numel(timingfile)      
            if ~isempty(timingfile{itrigifile}) 
                if isfield(NIRS.Dt.fir,'aux5')      
                   if numel(NIRS.Dt.fir.aux5)>=ifile
                        aux5 = NIRS.Dt.fir.aux5{ifile};
                   else
                        aux5 = [];
                   end
                else
                    aux5 = [];
                end
                timesec = timingfile{ifile};
                sample = round(NIRS.Cf.dev.fs*timesec);
                if size(sample,1)<size(sample,2)
                    sample = sample';
                end
                if ~isempty(sample)
                    sample = [ones(numel(sample),1)*Trigvalue,sample]; 
                    aux5 = [aux5;sample];
                end
                NIRS.Dt.fir.aux5{ifile}=aux5;
                clear aux5;
                end
            end
        else %remplacer aux5 par le trig manuel
            trigfile = [];         
                %Initialise AUX5
                if itrigfile==1
                    aux5 = [];
                    NIRS.Dt.fir.aux5{ifile} = [];
                end 
                %AJOUTE AU VALEURS 
                if isfield(NIRS.Dt.fir,'aux5')            
                    aux5 = NIRS.Dt.fir.aux5{ifile};
                else
                    aux5 = [];
                end             
                timesec = timingfile{ifile};
                sample = round(NIRS.Cf.dev.fs*timesec);           
                if size(sample,1)<size(sample,2) %sample doit être une column
                    sample = sample';
                end
                sample = [ones(numel(sample),1).*Trigvalue,sample]; 
                aux5 = [aux5;sample];
                NIRS.Dt.fir.aux5{ifile}=aux5;
                clear aux5;           
        end    
    end
end
save(job.NIRSmat{1,1},'NIRS');
out.NIRSmat = job.NIRSmat;