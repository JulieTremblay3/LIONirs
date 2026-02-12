function writetxtfile_asxlswrite(filename,raw)
%CONVERT CELL INTO A CVS file separate by semicolon ; 
%input similar as xlsx spreadsheet in windows
fid = fopen(filename,'w');
for iline = 1:size(raw,1)
    for jcol =  1:size(raw,2)
    if  ischar(raw{iline ,jcol})
        fprintf(fid,'%s;',raw{iline ,jcol});
    else
       fprintf(fid,'%d;', raw{iline ,jcol});
    end
    end
     fprintf(fid,'%s\r','');
end

fclose(fid);