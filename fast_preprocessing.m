function [dout,ok] = fast_preprocessing(d,  preprocess, fs)
%1  check for std deviation of the ch
%2  check for minimal intensity
%3  ensure that both wavelenght are valid

d0 = d;
NC = size(d,1);
samp_length = size(d,2);
ok = ones(NC,1);

%%% Normalization
for Idx = 1:NC
    meansub = mean(d(Idx,:)); %Normalize using mean over the all data
    d0(Idx,:) = d(Idx,:)./meansub;
end

%%% Rejection of noisy channels (based on norm or raw data, displays raw)
if preprocess.DCampenable
    for Idx = 1:NC/2 %For 830 nm channels     
        if  preprocess.STDmenu==1 %based on the median of random event 
            for i = 1:100
                ind = round(rand(1)*(samp_length - 10*fs - 1))+1;
                cf(i) = std(d0(Idx,ind:ind+round(10*fs)))/mean(d0(Idx,ind:ind+round(10*fs)));
            end
            cfgen = median(cf);
                if cfgen > preprocess.STDamp
                    ok(Idx) = 0;
                    ok(Idx+NC/2) = 0;
                end       
        elseif preprocess.STDmenu==2  %based on the whole recording segment      
            cf = std(d0(Idx,:))/mean(d0(Idx,:));
            if cf > preprocess.STDamp
                ok(Idx) = 0;
                ok(Idx+NC/2) = 0;
            end        
        end
    end
    for Idx = NC/2+1:NC %For 690 nm channels
        if ok(Idx) ~= 0
            for i = 1:100
                ind = round(rand(1)*(samp_length - 10*fs - 1))+1;
                cf(i) = std(d0(Idx,ind:ind+round(10*fs)))/mean(d0(Idx,ind:ind+round(10*fs)));
            end
            if  preprocess.STDmenu==1
            cfgen = median(cf);
                if cfgen > preprocess.STDamp   
                    ok(Idx) = 0;
                    ok(Idx-NC/2) = 0;                 
                end
            elseif preprocess.STDmenu==2
                cf = std(d0(Idx,:))/mean(d0(Idx,:));
                if cf > preprocess.STDamp
                    ok(Idx) = 0;
                    ok(Idx-NC/2) = 0;
                end
            end
        end
    end
end

%%% Rejection of channels of too low amplitude on based raw data.
if preprocess.STDampenable
    for Idx = 1:NC/2 %For 830 nm channels 
        if ok(Idx) ~= 0
         if mean(d(Idx,:))< preprocess.DCamp          
                ok(Idx) = 0;
                ok(Idx+NC/2) = 0;           
         end
        end
    end
    for Idx = NC/2+1:NC %For 690 nm channels
        if ok(Idx) ~= 0
             if mean(d(Idx,:))< preprocess.DCamp     
                ok(Idx) = 0;
                ok(Idx-NC/2) = 0;                         
            end
        end
    end
end

dout = d;


