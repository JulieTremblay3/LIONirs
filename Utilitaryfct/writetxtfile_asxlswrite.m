function writetxtfile_asxlswrite(filename,trialbysubjectbyele)
%CONVERT CELL INTO A TAB SPACE SEPARATED TXT FILE 
%input similar as xlswrite 
fid = fopen(filename,'w');
for iline = 1:size(trialbysubjectbyele,1)
    for jcol =  1:size(trialbysubjectbyele,2)
    if  ischar(trialbysubjectbyele{iline ,jcol})
        fprintf(fid,'%s\t',trialbysubjectbyele{iline ,jcol});
    else
       fprintf(fid,'%d\t', trialbysubjectbyele{iline ,jcol});
    end
    end
     fprintf(fid,'%s\r','');
end

fclose(fid);