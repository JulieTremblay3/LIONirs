% VideoMain
%Gyrus Angulaire
% clear all, close all
% path = 'H:\Dima\Lecture francais\Video\Video GyrusAngulaire\Mot\'
% file = 'Mot.avi'
% texttoadd = 'Word'
% path = 'H:\Dima\Lecture francais\Video\Video GyrusAngulaire\MotmoinsNonMot\'
% file = 'MotmoinsNonMot.avi'
% texttoadd = 'Word - Non Word'
% path = 'H:\Dima\Lecture francais\Video\Video GyrusAngulaire\NonMot\'
% file = 'NonMot.avi'
% texttoadd = 'Non Word'

% GyrusFusiforme
% path = 'H:\Dima\Lecture francais\Video\Video GyrusFusiforme\Mot\'
% file = 'Mot.avi'
% texttoadd = 'Word'
% path = 'H:\Dima\Lecture francais\Video\Video GyrusFusiforme\MotmoinsNonMot\'
% file = 'MotmoinsNonMot.avi'
% texttoadd = 'Word - Non Word'
% path = 'H:\Dima\Lecture francais\Video\Video GyrusFusiforme\NonMot\'
% file = 'NonMot.avi'
% texttoadd = 'Non Word'

% Video Parietal
% path = 'H:\Dima\Lecture francais\Video\Video Parietal\Mot\'
% file = 'Mot.avi'
% texttoadd = 'Word'
% path = 'H:\Dima\Lecture francais\Video\Video Parietal\MotmoinsNonMot\'
% file = 'MotmoinsNonMot.avi'
% texttoadd = 'Word - Non Word'
% path = 'H:\Dima\Lecture francais\Video\Video Parietal\NonMot\'
% file = 'NonMot.avi'
% texttoadd = 'Non Word'
[file,path] =uigetfile('.avi')
%  path = 'C:\data\LFC01_SM\Analyse1_pc00\T vs Baseline -3_-0\'
% file = 'HbO.avi'

filelegende = ['C:\data\LFC01_SM\','Legende-15_15.TIF']
texttoadd = 'Adults'
filein = [path,file];
fileout =[path,'r1',file];
dx = 1:500;
dy = 165:900;
dframe = 1:11;
fileinfo = aviinfo(filein);
dframe = 2:fileinfo.NumFrames;
load([path,'time.mat'],'-mat')
VideoRescale(filein,fileout,dx,dy,dframe,time(2:end))
filein = [path,'r1',file];
fileout = [path,'c',file];
VideoColorBar(filein,fileout,filelegende,dframe)

filein = [path,'c',file];
fileout = [path,'y',file];
posx = 100;
posy = 20;
VideoText(filein,fileout,texttoadd,posx,posy,dframe)

[video, audio] = mmread(fileout)
fileout = [path,'final',file,'.wmv']
mmwrite([fileout],audio,video);
cla
name = 'C:\data\LFC01_SM\Analyse01_pc00_spiketask\finalvideo_allspk.avi.wmv'
figure
[video, audio]=mmread(name,17)
imagesc(video.frames.cdata)
set(gca,'visible','off')