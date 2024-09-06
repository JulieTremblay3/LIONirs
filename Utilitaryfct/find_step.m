% function [stepmat] = find_step(d,id_step_dur)
% % Compare channels with an ideal step data serie to mark high correlation.
% % Inputs : Data matrix, duration of the ideal step (in data points EN SECONDES??!!)
% % Outputs : 0 & 1 stepmat matrix of bad and good intervals
% 
% id_step_dur = round(id_step_dur);
% if ~mod(id_step_dur,2) %If even, make it odd
%     id_step_dur = id_step_dur+1;
% end
% 
% x = ones(id_step_dur,1);
% x(round(id_step_dur/2)+1:end) = 2;
% x(round(id_step_dur/2)) = 1.5;
% 
% NC = size(d,1);
% win = floor(numel(x)/2);
% samp_length = size(d,2);
% 
% for Idx = 1:NC
%     d1 = d(Idx,:)';
% %     for t = 1:win
% %         rho(Idx,t) = corr(d1(1:t*2),x);
% %     end
%     for t = 1+win:samp_length-win
%         rho(Idx,t) = corr(d1(t-win:t+win),x);
%     end
% %     for t = size(d,2)-win+1:samp_length
% %         rho(Idx,t) = corr(d1(samp_length-(samp_length-t)*2:samp_length),d2(samp_length-(samp_length-t)*2:samp_length));
% %     end
%     
%     stepmat(Idx,:) = rho(Idx,:) >= 0.8;
% end
% 
% 
% 
% 
% 
% 
% 
% 
% figure;
% imagesc(rho);figure(gcf);
% colorbar
% title('Correlation between channels and ideal step');
% %Fenêtre mobile sur toutes les données
% % if corr(x,d(fenêtre)) > 0.95
% %   STEP!!!
% 
% for Idx = 1:10
%     figure;
%     hold on
%     subplot(3,1,1); 
%     hold on
%         plot(d1)
% %         plot_triggers(d1, trig_aux)
%     hold off
%     axis([0  samp_length  min(d1) max(d1)])
%     legend(['Channel ',num2str(Idx)]);
%     title(['Correlation for channel ', num2str(Idx)]);
%     subplot(3,1,2)
%     hold on
%         plot(x,'r')
% %         plot_triggers(x, trig_aux)
%     hold off
%     % axis([0  samp_length  min(x) max(x)])
%     legend('Ideal step');
%     subplot(3,1,3)
%     hold on
%         plot(rho(Idx,:),'k')
% %         plot_triggers(rho(Idx,:), trig_aux)
%     hold off
%     axis([0  samp_length  -1.1  1.1])
%     legend('Rho');
%     hold off
% end

