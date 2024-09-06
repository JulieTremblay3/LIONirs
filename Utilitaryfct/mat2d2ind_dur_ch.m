function [ind_dur_ch] = mat2d2ind_dur_ch(varargin)
% transform a binary matrice (1 mean artefact) ch x time in 3 columns vector, 
% column 1 : indice where we get ones, column 2 : the duration of consecutive
% one and column 3 : the channel. 

matrice = varargin{1};
% matrice =  zeros(20,20);
% matrice(2:8,1) = 1;
% matrice(13:15,4) = 1;
% matrice = matrice';

ind_dur_ch = [];
a=1;
for Idx = 1:size(matrice,1)
    indbad = find(matrice(Idx,:));
    if ~isempty(indbad)
        inddiff = diff(indbad);
        inddiff(end+1) = 0;
        step_dur = 0;
        for i = 1:numel(indbad) %For all steps just found
            if inddiff(i) == 1
                step_dur = step_dur+1;
            else
                ind_dur_ch(a,1) = indbad(i)-step_dur;
                ind_dur_ch(a,2) = step_dur+1;
                ind_dur_ch(a,3) = Idx;
                step_dur = 0;
                a = a+1;
            end
        end
    end
end
 