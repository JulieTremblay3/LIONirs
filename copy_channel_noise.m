function [noise] = copy_channel_noise(noise)
% Takes noise matrix (time x channels) and copy noise markers for one
% wavelength to the other of the same channel.

NC = size(noise,2);
noise(:,1:NC/2) = noise(:,1:NC/2) | noise(:,NC/2+1:end);
noise(:,NC/2+1:end) = noise(:,1:NC/2);

