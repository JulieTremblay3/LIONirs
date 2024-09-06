function [matgr1,idzone,idlist,idlabelall] = createROImat(matgr1,idzone,idlist,idlabelall)
 nbROI = numel(find(idzone));
        AROI = zeros( nbROI);
            izonefull = idzone;
         for izone = 1:numel(idzone);
            if  izonefull(izone)==0
                izonefull(izone)=  previous ;
            else
                previous = idzone(izone);
            end
         end
         for i=1: nbROI
            for j=1:nbROI
                idi = find(izonefull==i);
                idj = find(izonefull==j);
                 AROI(i,j) = nanmean(nanmean(matgr1( idlist(idi),idlist(idj)),1),2);
            end
        end
        idlist = 1:nbROI;
        if numel(find(idzone)) <numel(idlabelall)
        idlabelall = idlabelall(find(idzone));
        end
        matgr1 = AROI; 
        idzone = 1:nbROI;
end