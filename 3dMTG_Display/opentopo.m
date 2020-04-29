function vColor = opentopo(filename)
%Permet d'ouvrir les images topographiques, quelle soient en .vcolor ou en .img

[dir,fil,ext]=fileparts(filename);
if strcmp(upper(ext),upper('.vcolor'))
    load(filename,'-mat'); 
elseif  strcmp(ext,'.img')
    Vtest = spm_vol(filename);
    try
        vColormat = spm_slice_vol(Vtest,spm_matrix([0 0 1]),Vtest.dim(1:2),0);
        vColor = vColormat(:);
    catch
        msgbox('File can''t not be open, verify if SPM is installed or if the file format is ok')
    end
end


