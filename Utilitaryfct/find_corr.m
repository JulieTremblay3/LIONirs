% function [corr_ch] = find_corr(mean_diff,stepmat)
% % For a specified time interval in the data series, find_corr looks for
% % correlated channels.
% % Inputs : Matrix of data points (or of mean difference), matrix for markers
% % 
% 
% %Correlation between channels for steps that have
% %just been detected
% steps = find(ind_dur_ch(:,3)==Idx);
% for i = 1:numel(steps)   %For all found steps of the current channel
%     interval = [ind_dur_ch(steps(i),1) ind_dur_ch(steps(i),1)+ind_dur_ch(steps(i),2) ind_dur_ch(steps(i),4)]; % interval = [ind_start ind_end type]
%     d1 = d(Idx, interval(1):interval(2));
%     for Ich = 1:NC   %Run through all the channels to find correlation for the selected interval
%         if NIRS.Cf.H.C.ok(Idx,f) ~= 0;   %Select a valid channel
%             d2 = d(Ich, interval(1):interval(2));
%             if corr(d1,d2) >= 0.95
%                 ind_dur_ch(a,1) = interval(1);
%                 ind_dur_ch(a,2) = interval(2);
%                 ind_dur_ch(a,3) = Ich;
%                 ind_dur_ch(a,4) = interval(3);
%             end
%         end
%     end
% end
%                         