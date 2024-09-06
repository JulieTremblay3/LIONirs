%GyrusAngulaire
pathMot = 'H:\Dima\Lecture francais\Video\Video GyrusAngulaire\Mot\'
pathNonMot = 'H:\Dima\Lecture francais\Video\Video GyrusAngulaire\NonMot\'
pathout = 'H:\Dima\Lecture francais\Video\Video GyrusAngulaire\MotmoinsNonMot\'
%PARIETAL
pathMot = 'H:\Dima\Lecture francais\Video\Video Parietal\Mot\'
pathNonMot = 'H:\Dima\Lecture francais\Video\Video Parietal\NonMot\'
pathout = 'H:\Dima\Lecture francais\Video\Video Parietal\MotmoinsNonMot\'
%Video GyrusFusiforme
pathMot = 'H:\Dima\Lecture francais\Video\Video GyrusFusiforme\Mot\'
pathNonMot = 'H:\Dima\Lecture francais\Video\Video GyrusFusiforme\NonMot\'
pathout = 'H:\Dima\Lecture francais\Video\Video GyrusFusiforme\MotmoinsNonMot\'
for itime = 1:40
    all = [];
    for iname=1:numel(namesubject);
            Mot = load([pathMot,'HbOAVG', num2str(itime),'(s).vcolor'],'-mat');
            NonMot = load([pathNonMot,'HbOAVG', num2str(itime),'(s).vcolor'],'-mat');
    end
    vColor = Mot.vColor - NonMot.vColor;
    [nameout]= [pathout,'HbOAvg',num2str(itime),'(s).vColor']
    save([nameout],'vColor','-mat')
    
    for iname=1:numel(namesubject);
           Mot = load([pathMot,'HbRAVG', num2str(itime),'(s).vcolor'],'-mat');
           NonMot = load([pathNonMot,'HbRAVG', num2str(itime),'(s).vcolor'],'-mat');
    end
    vColor = Mot.vColor - NonMot.vColor;
    [nameout]= [pathout,'HbRAvg',num2str(itime),'(s).vColor']
    save([nameout],'vColor','-mat')
end


