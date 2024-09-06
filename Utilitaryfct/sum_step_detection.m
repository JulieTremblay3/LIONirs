function [stepmat] = sum_step_detection(d,thr_ind,corr_thr)
% Sums all channels and uses difference of means to detect steps in it.
% Correlation between the sum and each channel is then performed to mark
% bad intervals in each channel.
% Inputs : matrix of data, threshold factor for step detection, correlation
% threshold
% Outputs : data-size matrix of 0 and 1 (no step or step)

step = 15;
samp_length = size(d,2);

ch_sum = sum(d,1);

meantot = [];
mean_diff = [];
stepmat = zeros(size(d,1),size(d,2));

%Moving window average
 for i = 1:step %first point all same value
    meantot(i) = mean(d(Idx,1:step));
end
for i = 1+step:samp_length-step %window between
    meantot(i) = mean(d(Idx,i-step:i+step));
end
for i = samp_length-step:(samp_length+1) %end of the window all same value
   % meantot(i) = mean(d(Idx,i:size(d,2)));
   meantot(i) = mean(d(Idx,samp_length-step:samp_length)); 
end 
mean_diff = diff(meantot); 

 %Definition du threshold                         
zsc = zscore(abs(mean_diff));
indok = find(abs(zsc)<3);
valdata = median(abs(mean_diff(indok)));
threshold = valdata*thr_ind;                  
                            

%Find data above threshold
stepmat_sum = (abs(mean_diff) > threshold); %All positive elements are steps

% %Correlation with an ideal step NOT USED DUE TO DETECTION OF PULSE
% dcorr = ch_sum';
% for t = 1+win:samp_length-win
%     rho(t) = corr(dcorr(t-win:t+win),x);
% end
% rho(samp_length-win+1:samp_length-1) = 0; %Padding
% stepmat_sum = stepmat_sum + (abs(rho) >= 0.85);

% Looking for channels that are correlated to the steps of ch_sum
indbad = find(stepmat_sum);
if ~isempty(indbad)
    inddiff = diff(indbad);
    inddiff(end+1) = 0;
    step_dur = 0;
    for i = 1:numel(inddiff) %For all steps just found
        if inddiff(i) == 1
            step_dur = step_dur+1;
        else
            temp_ind = indbad(i)-step_dur;
            temp_dur = step_dur+1;
            padding = 15;
            interval = [temp_ind temp_ind+temp_dur]; % interval = [ind_start ind_end]
                if interval(2) == size(d,2); interval(2) = size(d,2)-1; end;
            if interval(1)-padding < 1
                interval(1) = padding+1;
            end
            if interval(2)+padding > size(d,2)
                interval(2) = size(d,2)-padding;
            end
            d1 = ch_sum(interval(1)-padding:interval(2)+padding)';
            for Ich = 1:size(d,1)   %Run through all the channels to find correlation for the selected interval
                d2 = d(Ich, interval(1)-padding:interval(2)+padding)';
                if corr(d1,d2) >= corr_thr
                    stepmat(Ich,interval(1):interval(2)) = 1;
                end
            end
        end
     end
end