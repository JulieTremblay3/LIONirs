function [priorite muxchosen muxdispo]= auto_attribution_srs_label_4(indp,typemuxsrs,muxdispo,intsourceattribuer)

bb = typemuxsrs(indp); %1=droite 2=gauche 
right=find(bb==1);
left = find(bb==2);

if numel(right)==numel(left) %SI egal de gauche et de droite on prend les mux en ordre croissant
    if mod(intsourceattribuer,2)
        list = find(muxdispo(1:8));
    else
        list = find([zeros(8,1);muxdispo(9:16)]);
    end
    if isempty(list)
        list =find(muxdispo);
    end
    muxchosen=list(1);    
elseif numel(right)>numel(left) %Plus de droit que de gauche on prend les muxs de 1 a 8 a moins qu'il soit utilisé
        list = find(muxdispo(1:8));
        if ~isempty(list)
             muxchosen=list(1);
        else 
            list = find(muxdispo)
            muxchosen=list(1);
        end

elseif numel(right)<numel(left) %Plus de gauche que de droit on prend les muxs de 9 a 16 a moins qu'il soit utilisé
        list = find([zeros(8,1);muxdispo(9:16)]);
        if ~isempty(list)
             muxchosen=list(1);
        else 
            list = find(muxdispo)
            muxchosen=list(1);
        end
end 
muxdispo(muxchosen)=0;
%Place mux 1 2 3 
priorite = zeros(4,1);

if isempty(left)
   priorite(right(2)) = 2; %si pas de trou gauche on va placer un c a droite
else
    priorite(left(1)) = 2; %mux c
end

if isempty(right)
    priorite(left(2)) = 1; %si pas de trou droit on va placer un a a gauche
else
    priorite(right(1)) = 1; %mux a
end
    

unplaced = find(priorite==0);
if numel(unplaced)==2
    priorite(unplaced(1))=3; %mux e
    priorite(unplaced(2))=4; %mux g
end


%OLD VERSION
% [sortbb ix] = sort(bb);
% bad = find(sortbb  ~= (1:numel(bb))');
% bb(ix(bad)) = bad;

