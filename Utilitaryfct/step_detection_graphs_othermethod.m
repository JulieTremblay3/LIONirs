%%%GRAPHIQUES DE RUN_STEP_DETECTION %%%%%%%%%%%%%%



% %%%%%%%%%%%%%%%%%%%%%% GRAPHS %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% d2 = d';
% % d2(:,4) = d2(:,4)-mean(d2(:,4));
% % 
% % figure; 
% % hold on
% % plot(d2(:,4))
% % plot(50*mean_diff(4,:), 'g')
% % hold off
% 
% steps_ind = find(ind_dur_ch(:,4)==0); %Beginning of steps
% badpointsmean = mean_diff(ind_dur_ch(steps_ind,3),ind_dur_ch(steps_ind,1))';
% 
% figure; 
% hold on
% % plot(mean_diff(1,:))
% % plot(mean_diff(2,:), 'r')
% % plot(mean_diff(3,:), 'k')
% plot(mean_diff(4,:), '.')
% plot(ind_dur_ch(steps_ind,1),badpointsmean,'x', 'MarkerSize', 15,'MarkerEdgeColor','r')
% hold off
% 
% %%% Display of found steps on the raw data
% t = 1:size(d,2);
% x = (min(d(4,:))-1)*ones(size(d,2),1);
% % x(ind_dur_ch(steps_ind,1)) = max(d(4,:));
% for i = 1:size(steps_ind)
%     x(ind_dur_ch(steps_ind(i),1):ind_dur_ch(steps_ind(i),1)+ind_dur_ch(steps_ind(i),2)) = max(d(4,:));
% end
% 
% figure; 
% hold on
% plot(d(4,:), '.')
% plot(t,x,'r')
% legend('ch 1','ch 2','ch 3','ch 4')
% hold off
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%%%AUTRE MÉTHODE DE CALCUL DE DIFF DE MOYENNE



%%%% OTHER METHOD
% for Idx = 1:NC
%     waitbar(Idx/NC,hwaitbar,'Looking for steps in raw data');
%     for samp_pt = step+1:size(d,2)-step
%         mean_diff(Idx,samp_pt-step) = mean(d(Idx,samp_pt:samp_pt+step)) - mean(d(Idx,samp_pt-step:samp_pt));
%         if abs(mean_diff(Idx,samp_pt-step)) > threshold
%                if samp_pt > step+1 && abs(mean_diff(Idx,samp_pt-step-1)) > threshold
%                    ind_dur_ch(a-1,2) = ind_dur_ch(a-1,2)+1;
% %                    step_dur(Idx,cl-1) = step_dur(Idx,cl-1)+1;
%                else
%                    ind_dur_ch(a,1) = samp_pt;
%                    ind_dur_ch(a,2) = 1;
%                    ind_dur_ch(a,3) = Idx;
% %                    step_ind(Idx,cl) = samp_pt;
% %                    step_dur(Idx,cl) = 1;
%                    a = a+1;
%                end
%                %NIRS.Dt.fir.pp(Idx).bpi{f,1} = 1;
%                %NIRS.Dt.fir.pp(Idx).bpd{f,1} = 1;
%         end
%     end
% end


