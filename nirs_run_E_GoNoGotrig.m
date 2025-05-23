function out = nirs_run_E_GoNoGotrig(job)
%Check The good and false response to the gonogo taks
% Check the pourcentage d'omission absence de 3 pour les blocs go


%Ajouter des trigs pour voir facilement le pourcentage de Go réussi entre
%100 et 200, 200 étant Go réussi à 100 %

% Ajouter des trigs pour voir facilement le pourcentage de NoGo réussi
% entre
% 200 et 300


for filenb=1:size(job.NIRSmat,1) %Loop over all subjects
        NIRS = [];
        load(job.NIRSmat{filenb,1});
        [dir2,tmp,tmp] = fileparts(job.NIRSmat{filenb,1});

        %use last step of preprocessing
        lst = length(NIRS.Dt.fir.pp);
        rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
        NC = NIRS.Cf.H.C.N;
        fs = NIRS.Cf.dev.fs;
% 
offsetNIRSEEG = NIRS.Dt.EEG.pp(end).sync_timesec{1}
fileEEG =  NIRS.Dt.EEG.pp(end).p{1}
 [~,EEG.infoBV,EEG.marker,EEG.ind_dur_ch]=fopen_EEG(fileEEG);

TimeEEGs = EEG.infoBV.SamplingInterval*EEG.ind_dur_ch(:,1)/1000000-offsetNIRSEEG;
aux5 = NIRS.Dt.fir.aux5{1}; 
TimeNIRSs =aux5(:,2).*1./NIRS.Cf.dev.fs;


%S11 bloc de GO
% voir le nombre de bonne réponse
% S  3 près affichage de S  2
disp('Bloc S11 tache GO')
disp(['Voir le nombre de bonne réponse event S 2 trig de réponse S  3 présent ' ...
    'Ajouter un trig entre 100 et 200 correspondant au pourcentage  '])
disp('Bloc S10 tache NOGO')
disp(['Voir le nombre de bonne réponse NOGO event S 1 trig sans réponse ou' ...
    'd''erreur  4 réponse à un no go erreur du participant'])

iGo = []; 
iNoGo= [];
iend = [];
for im=1:size(EEG.marker,1)% S 11
    if TimeEEGs(im)>=min(TimeNIRSs)& TimeEEGs(im)<=max(TimeNIRSs)
    if strcmp(EEG.marker(im,2), 'S 11')
        iGo = [iGo,im];
    end
    if strcmp(EEG.marker(im,2), 'S 10')
         iNoGo= [iNoGo,im];    
    end
     if strcmp(EEG.marker(im,2), 'S 12')
         iend= [iend,im];    
     end
    end
end

iallbloc = sort([iGo, iNoGo, iend(end)]);
%ajout d'un bloc de fin de 754-


%vérifier que ses temps font partie du blocs NIRS
%iallblocshort = iallbloc(find(TimeEEGs(iallbloc)>TimeNIRSs(1) & TimeEEGs(iallbloc)<TimeNIRSs(end)));

 
%Check le pourcentage de bonne réponse go et nogo
for ibloc = 1:(numel(iallbloc)-1);
    Reptot = 0;
    Repok = 0;
    NoGotot = 0;
    NoGook = 0;
    for istim=iallbloc(ibloc):iallbloc(ibloc+1);
        if strcmp(EEG.marker(istim,2),'S  2'); %stim  go
            Reptot = Reptot+1 ;
            if strcmp(EEG.marker(istim+1,2),'S  3'); %bonne réponse pour un go
            Repok = Repok + 1; 
            end
        end
        if strcmp(EEG.marker(istim,2),'S  1'); %stim  nogo
            NoGotot =  NoGotot+1;
            if strcmp(EEG.marker(istim+1,2),'S  4')              
            else
                NoGook = NoGook +1;
            end
        end 
    end
    iallReptot(ibloc)=Reptot;
    iallRepok(ibloc)= Repok;
    iallNoGotot(ibloc)= NoGotot;
    iallNoGook(ibloc)=  NoGook; 
end
pourcentageGoOK= round(iallRepok./iallReptot*100');
pourcentageNoGoOK =round(iallNoGook./iallNoGotot*100');

%Ajouter des trigs pour voir facilement le pourcentage de Go réussi entre
%100 et 200, 200 étant Go réussi à 100 %
aux5score =[];
for ibloc = 1:numel(iallNoGook);
   aux5score = [aux5score; 100+pourcentageGoOK(ibloc), round((TimeEEGs(iallbloc(ibloc)) - 2) *NIRS.Cf.dev.fs) ];
end


aux5scorenogo =[];
for ibloc = 1:numel(iallNoGook) 
    if ~isnan(pourcentageNoGoOK(ibloc))
    aux5scorenogo = [aux5scorenogo; 200+pourcentageNoGoOK(ibloc), round((TimeEEGs(iallbloc(ibloc)) - 4) *NIRS.Cf.dev.fs) ];

    end
end
[filepath,name,ext] = fileparts(job.NIRSmat{filenb,1})

id = strfind(filepath, '\');
filepath(id) = [];
id = strfind(filepath, 'Analyses_Kass_2025');
filename = filepath(id+18:end);
try
fid = fopen(fullfile(job.output_path,[filename,'.txt']),'w')
disp(['Create as report the file', fullfile(job.output_path,[filename,'.txt']),...
      'first column percentage GoOk, second column percentage NoGoOK'])



 

for ibloc = 1:numel(pourcentageGoOK)
    fprintf(fid, '%d\t%d\r ', pourcentageGoOK(ibloc), pourcentageNoGoOK(ibloc) );
% Ajouter des trigs pour voir facilement le pourcentage de NoGo réussi
% entre
% 200 et 300
end
fclose(fid);
catch
end

            NIRS.Dt.fir.aux5{1}=[aux5;aux5score;aux5scorenogo];
           save(fullfile(dir2,'NIRS.mat'),'NIRS');
           job.NIRSmat{1} =fullfile(dir2,'NIRS.mat');

end
out.NIRSmat = job.NIRSmat;
