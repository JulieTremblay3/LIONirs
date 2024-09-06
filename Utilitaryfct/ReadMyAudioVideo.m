function [video,audio]=ReadMyAudioVideo(filenamevideo,timeid)
%inspire form mmread toolbox
%use matlab VideoReader to use codec install on the machine
%32 bit matlab use 32 bit codec
%64 bit matlab use 64 bit codec

tstart = timeid(1);
tstop = timeid(end);
v=VideoReader(filenamevideo,'CurrentTime',tstart);
nbframe = round(( tstop -  tstart)*v.FrameRate);
video.width= v.Width;
video.height= v.Height;
video.nrFrameTotal= 1;
for iframe = 1:nbframe
    vidFrame = readFrame(v);
    video.frames(1,iframe).cdata = vidFrame;
    video.frames(1,iframe).colormap = [];
end
video.rate=v.FrameRate;
video.totalDuration=v.Duration;
video.times= tstart:1/v.FrameRate: tstop;

%READ AUDIO by video frame 
try
    info = audioinfo(filenamevideo);
    if tstart ==0
        tstart =1;
    end
    [y,Fs] = audioread(filenamevideo,round(([tstart tstop ])*info.SampleRate));
    audio.bits = info.BitsPerSample;
    audio.rate = info.SampleRate;
    audio.data = y;
    audio.totalDuration = info.Duration;
    audio.nrChannels = 2;
    audio.times = 1/video.rate : 1/video.rate : (1/video.rate*nbframe);
catch
    audio = [];
end
