%pathNonMot = 'H:\Dima\Lecture francais\Video\Video_NonMot\'
%pathoutNonMot = 'H:\Dima\Lecture francais\Video\Video GyrusAngulaire\NonMot\'
%pathoutMotmoinsNonMot = 'H:\Dima\Lecture francais\Video\Video GyrusAngulaire\MotmoinsNonMot\'
% MOT Gyrus Angulaire
namesubject = {'LFP05';
    'LFP07';
    'LFP14';
    'LFP15';
    'LFP18';
    'LFP19'}
% path = 'H:\Dima\Lecture francais\Video\Video_Mot\'
% pathout = 'H:\Dima\Lecture francais\Video\Video GyrusAngulaire\Mot\'
% NON MOT Gyrus Angulaire
% path = 'H:\Dima\Lecture francais\Video\Video_NonMot\'
% pathout = 'H:\Dima\Lecture francais\Video\Video GyrusAngulaire\NonMot\'

%MOT PARIETAL
namesubject = {'LFP07';
    'LFP08';
    'LFP13';
    'LFP14';
    'LFP15';
    'LFP17';
    'LFP18'}
% path = 'H:\Dima\Lecture francais\Video\Video_Mot\'
% pathout = 'H:\Dima\Lecture francais\Video\Video Parietal\Mot\'
% NON MOT Gyrus Angulaire
path = 'H:\Dima\Lecture francais\Video\Video_NonMot\'
pathout = 'H:\Dima\Lecture francais\Video\Video Parietal\NonMot\'

%MOT Gyrus Fusiforme
namesubject = {'LFP05';
    'LFP12';
    'LFP14';
    'LFP15';
    'LFP18';
    'LFP19';}
% path = 'H:\Dima\Lecture francais\Video\Video_Mot\'
% pathout = 'H:\Dima\Lecture francais\Video\Video GyrusFusiforme\Mot\'
% NON MOT Gyrus Angulaire
% path = 'H:\Dima\Lecture francais\Video\Video_NonMot\'
% pathout = 'H:\Dima\Lecture francais\Video\Video GyrusFusiforme\NonMot\'

for itime = 1:40
    %     all = [];
    try
    allHbO = [];
    labeltime = sprintf('%03.0f',itime)
    for iname=1:numel(namesubject);
        load([path,namesubject{iname,1},'.aviHbOTime   ', labeltime,'(s).vcolor'],'-mat');
        allHbO = [allHbO,vColor];
    end
    vColor = mean(allHbO')';
    [nameout]= [pathout,'HbOAvg',labeltime,'(s).vColor']
    save([nameout],'vColor','-mat')
    catch 
        'HbO not averaged'
    end
    try 
    allHbR = [];
    for iname=1:numel(namesubject);
        load([path,namesubject{iname,1},'.aviHbRTime   ', labeltime,'(s).vcolor'],'-mat');
        allHbR = [allHbR,vColor];
    end
    vColor = mean(allHbR')';
    [nameout]= [pathout,'HbRAvg',labeltime,'(s).vColor']
    save([nameout],'vColor','-mat')
    catch
        'HbR not averaged'
    end
    try 
    allHbT = [];
    for iname=1:numel(namesubject);
        load([path,namesubject{iname,1},'.aviHbRTime   ', labeltime,'(s).vcolor'],'-mat');
        allHbT = [allHbT,vColor];
    end
    vColor = mean(allHbT')';
    [nameout]= [pathout,'HbTAvg',labeltime,'(s).vColor']
    save([nameout],'vColor','-mat')
    catch 
        'HbT not averaged'
    end
end


