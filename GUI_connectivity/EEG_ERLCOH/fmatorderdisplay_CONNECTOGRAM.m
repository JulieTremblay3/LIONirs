function [idlabelall,idzone,idlist]=fmatorderdisplay(xlsfile, listelectrode)
% find the indice to display the matrix in a comprehensive way 
% according to the xls label order and the electrode identification
% Use plotCOH
% imagesc(A(idlist,idlist));
% set(gca,'xtick', find(idzone));
% set(gca,'xticklabel', idlabelall);
% set(gca,'ytick', find(idzone));
% set(gca,'yticklabel', idlabelall);
% set(gca,'fontsize',14);
% axis square
% caxis([cmin,cmax])
% colormap(jet)
% title([titlelabel])
%elefile = [pathinfo,'LEFTRIGHTelePLG.ele']

%%Prepare plot matrix order
[num, txtzone]=xlsread(xlsfile);
% [num, txtzone]=xlsread([pathinfo,'ZoneEEGLeftRight_Connectogramme.xls']);
idlist  =[];idlabelall = [];idzone = [];
for izone = 1:size(txtzone,2)
    for ichzone = 3:size(txtzone,1)
        ich = txtzone{ichzone,izone};
        idch = strmatch(ich,listelectrode,'exact');
        
        idlist = [idlist, idch];
        if ichzone==3
            idzone =[idzone, izone];
            idlabelall = [idlabelall, {txtzone{1,izone}}];
            
        elseif ~isempty(idch)
            idzone =[idzone, 0];  
            idlabelall = [idlabelall, {''}];
        end
    end
   %  idlabelall = [idlabelall, {txtzone{1,izone}}];
end 

