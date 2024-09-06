function idxc = find_idx_color(ml,indml,nbcolor)
%Fonction attribuant la couleur en fonction du canal selectionné


% if (PMI{currentsub}.plotLst(ind) > numel(PMI{currentsub}.data(cf).MeasListAct)/2)
%      idxc = PMI{currentsub}.plotLst(ind) - numel(PMI{currentsub}.data(cf).MeasListAct)/2;
% else
%      idxc = PMI{currentsub}.plotLst(ind); 
% end
% 
% idxc= mod((numel(PMI{currentsub}.color)/3),idxc);
% if idxc == 0 
%      idxc=1;
% end

ind = find(ml(:,1) == ml(indml,1) & ml(:,2) == ml(indml,2));
idxc= mod(ind(1),nbcolor);
if idxc == 0 
     idxc=nbcolor;
end


