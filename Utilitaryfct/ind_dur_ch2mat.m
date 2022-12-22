%Convert the marking indice dur channel format in a binary matrix nele x nsample 
function mat = ind_dur_ch2mat(ind_dur_ch, nsample,nch)
mat = logical(zeros(nsample,nch ));
if ~isempty(ind_dur_ch) 

for Idx = 1:nch
    mrks = find(ind_dur_ch(:,3)==Idx);
    ind = ind_dur_ch(mrks,1);
    indf = ind + ind_dur_ch(mrks,2) - 1;
    if ~isempty(ind)
        try
            for i = 1:numel(ind)
                mat(ind(i):indf(i),Idx) = 1;
            end
        catch
            disp('Noise reading problem')
        end
    end
end

%all idx ==0 mark for all channel
    mrks = find(ind_dur_ch(:,3)==0);
    ind = ind_dur_ch(mrks,1);
    indf = ind + ind_dur_ch(mrks,2) - 1;
    if ~isempty(ind)
        try
            for i = 1:numel(ind)
                mat(ind(i):indf(i),:) = 1;
            end
        catch
            disp('Noise reading problem')
        end
    end

end