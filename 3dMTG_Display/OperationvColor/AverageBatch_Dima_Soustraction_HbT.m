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

%Video HBT
pathMot = 'H:\Dima\Lecture francais\Video_Mot_Hbt\'
pathNonMot = 'H:\Dima\Lecture francais\Video_NonMot_HbT\'
pathout = 'H:\Dima\Lecture francais\VideoHbT\'

for itime = 1:25
    all = []; 
    Mot = load([pathMot,'HbTAVG', num2str(itime),'(s).vcolor'],'-mat');
    NonMot = load([pathNonMot,'HbTAVG', num2str(itime),'(s).vcolor'],'-mat');
    vColor = Mot.vColor - NonMot.vColor;
    [nameout]= [pathout,'HbTAvg',num2str(itime),'(s).vColor']
    save([nameout],'vColor','-mat')
    
end


