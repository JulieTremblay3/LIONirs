function plot_triggers(d,trig)
% trig is a Nx2 matrix of the trigger description (e.g. 255) in the first
% column, and indice point in the second column.
%

for trig_nb = 1:size(trig,1)
    trig_ind = trig(trig_nb,2);
    plot([trig_ind, trig_ind],[min(min(d)), max(max(d))], 'g');
end