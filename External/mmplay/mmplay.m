% function mmplay(...options...)
% 
% mmplay is a replacement for the movie command as well as plays movies 
% read by mmread.  It can play both audio and video simultaneously.  This 
% uses the Windows DirectX infrastructure to generate the movie, so other 
% OSs are out of luck.  The video will appear in the current figure (or pop
% one up if it doesn't exist) unless 'fullscreen' is specified.  The video 
% will cover the entire window, including the toolbar and menubar.  If this
% bothers you turn off the menubar with:
% set(gcf,'MenuBar','none');
% 
% INPUT:
%   video structure:  The video structure matches the output of mmread.  At
%     a minimum it must have 4 fields "frames", "times", "height" and
%     "width".  The "frames" field must be a struct array with a field
%     "cdata" that contains the raw frame data encoded as height by width
%     by color(3) as UINT8s.  The "times" field contains the time stamps of
%     the data stored in frames.  "times" and "frames.cdata" must be the
%     same length.
% 
%   audio structure:  The audio structure matches the output of mmread.  At
%     a minimum it must have 3 fields "data", "rate" and "times".  The
%     "data" field is a matrix nrSamples by nrChannels (the same format was
%     wavread/wavplay).  The field "rate" is the sampling rate of the data
%     in Hz, eg. 44100.  The field "frames" is used to specify the time
%     that audio should start, the rest of the time is extrapolated based
%     upon the "rate" and the nrSamples.
% 
%   'fullscreen':   play in fullscreen mode.
% 
% EXAMPLES:
%   [video,audio] = mmread('myfile');
%   mmplay(video,audio);
% 
%   %play fullscreen
%   mmplay(video,'fullscreen');
% 
%   % play just a subsample (the time range 10 to 20s) of a movie using mmread
%   [video, audio] = mmread('your file',[],[10 20]);
% 
%   subtract 10 seconds off the time stamps so that the audio and video
%   will start at the beginning of the movie.
%   video.times = video.times - 10;
%   audio.times = audio.times - 10;
% 
%   mmplay(audio,video); % play the movie...
% 
%   % make a movie ourselves...
%   f=figure;
%   x = -pi:.1:pi;
%   set(gca,'xlim',[-length(x) length(x)],'ylim',[-length(x) length(x)]);
%   for i=1:length(x)
%       patch(sin(x)*i,cos(x)*i,[abs(cos(x(i))) 0 0]);
%       v.frames(i) = getframe(gca);
%       v.times(i) = i/30; % display at 30fps
%   end
%   v.width=size(v.frames(1).cdata,2);
%   v.height=size(v.frames(1).cdata,1);
%   close(f);
%   mmplay(v);
% 
%   %as a replacement for the 'movie' command
%   change: movie(M)
%   to:
%   v.frames = M;
%   v.times = (1:length(M))/12; %movie defaults to 12fps
%   v.width = size(M(1).cdata,2);
%   v.height = size(M(1).cdata,1);
%   mmplay(v);
