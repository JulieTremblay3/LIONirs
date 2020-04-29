function [d] = interpolate_bad(d,ind_dur_ch)
% Interpolate bad intervals for filtering.
% Inputs : matrix of data channel x sample, ind_dur_ch matrix which contains only the intervals to
% interpolate
% Outputs : matrix of data with interpolated intervals
if ~isempty(ind_dur_ch)    
    maxpoint  = ind_dur_ch(:,1)+ind_dur_ch(:,2);
    badind = find(maxpoint>size(d,2));
    if ~isempty(badind)
        disp(['Warning  marker : ' num2str(badind') ' are out of range'])
        ind_dur_ch(badind,2)=size(d,2)- ind_dur_ch(badind,1);
    end
    
    for i = 1:size(ind_dur_ch,1)
        try
            ch = ind_dur_ch(i,3);
            ind = ind_dur_ch(i,1);
            y1 = d(ch,ind);
            if isnan(y1)&ind==1
                y1 = nanmean(d(ch,:));
            end
            while isnan(y1) && ind ~= 1 %To avoid using padding NaN
                ind = ind-1;
                y1 = d(ch,ind);
            end
            
            dur = ind_dur_ch(i,2);
            y2 = d(ch,dur+ind_dur_ch(i,1));
            
            while isnan(y2) && dur+ind_dur_ch(i,1) ~= size(d,2)
                dur = dur+1;
                y2 = d(ch,dur+ind_dur_ch(i,1));
            end
            if isnan(y2)
                y2 = nanmean(d(ch,:));
            end
            
            dur = dur + ind_dur_ch(i,1)-ind;

            interval = ind:ind+dur;

            %y = ax + b
            a = (y2-y1)/dur;
            b = y1 - a*ind;
            interp = a.*interval + b;

            d(ch,interval) = interp;
        catch
            disp('Failed to interpolate.');
        end
    end
end







