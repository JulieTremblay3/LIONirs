function Write_reviewlist(oHelmet)

sMtg = get_Mtg(oHelmet);  
vHoles = get_vHoles(oHelmet);
Ptot = sort([sMtg.v_pSrc,sMtg.v_pDet, sMtg.v_pEle ]);

[name path]=  uiputfile('.txt');
fid  = fopen([path,name],'w');
for indPtot = 1:numel(Ptot)
   if find(sMtg.v_pSrc==Ptot(indPtot));
       namehole = vHoles(Ptot(indPtot)).Label;
       if numel(namehole)<4
           namehole = [namehole,' '];
       end
       src_n = sMtg.v_HolesMtg(Ptot(indPtot));
       src_s = src_n2src_s(src_n);
       fprintf(fid,'%s\t\t%s\n', namehole,src_s');
   elseif find(sMtg.v_pDet==Ptot(indPtot))
       namehole = vHoles(Ptot(indPtot)).Label;
       if numel(namehole)<4
           namehole = [namehole,' '];
       end
       src_n = sMtg.v_HolesMtg(Ptot(indPtot));
       src_s = src_n2src_s(src_n);
       fprintf(fid,'%s\t\t%s\n', namehole, src_s');
   elseif find(sMtg.v_pEle==Ptot(indPtot))
       namehole = vHoles(Ptot(indPtot)).Label;
       if numel(namehole)<4
           namehole = [namehole,' '];
       end
       src_s = sMtg.v_HolesEle{Ptot(indPtot)};
       fprintf(fid,'%s\t\t%s\n', namehole,src_s);
   end
end
fclose(fid)
