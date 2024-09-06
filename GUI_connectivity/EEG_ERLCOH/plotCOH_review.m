function plotCOH_review(A,idhalf,idlabelall,idlist,idzone,cmin,cmax,listelectrode,titlelabel)
matgr1 = zeros(numel(listelectrode));
matgr1(idhalf)=A ;
%imagesc(matgr1(:,:));
matgr1 = matgr1 +flipud(rot90(matgr1))
 %imagesc(matgr1(:,:));
imagesc(matgr1(idlist,idlist));
set(gca,'xtick', find(idzone));
set(gca,'xticklabel', idlabelall);
set(gca,'ytick', find(idzone));
set(gca,'yticklabel', idlabelall);
set(gca,'fontsize',14);
axis square
caxis([cmin,cmax])
colormap(jet)
title([titlelabel])



