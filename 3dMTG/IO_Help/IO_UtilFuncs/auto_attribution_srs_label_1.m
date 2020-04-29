function [priorite muxchosen muxdispo]= auto_attribution_srs_label_1(indp,typemuxsrs,muxdispo)

bb = typemuxsrs(indp); %1=droite 2=gauche 
right=find(bb==1);
left = find(bb==2);

if numel(right)==numel(left) %SI egal de gauche et de droite on prend les mux en ordre croissant
    list = find(muxdispo);
    muxchosen=list(1);
elseif numel(right)>numel(left) %Plus de droit que de gauche on prend les muxs de 1 a 8 a moins qu'il soit utilisé
        list = find(muxdispo(1:8));
        if ~isempty(list)
             muxchosen=list(1);
        else 
            list = find(muxdispo);
            muxchosen=list(1);
        end
elseif numel(right)<numel(left) %Plus de gauche que de droit on prend les muxs de 9 a 16 a moins qu'il soit utilisé
        list = find([zeros(8,1);muxdispo(9:16)]);
        if ~isempty(list)
             muxchosen=list(1);
        else 
            list = find(muxdispo);
            muxchosen=list(1);
        end
end 
muxdispo(muxchosen)=0;
%Place mux 1 2 
priorite = zeros(1,1);

if isempty(left)
    priorite(right(1)) = 1; %mux a
else 
    priorite(left(1)) = 2; %mux c
end
 