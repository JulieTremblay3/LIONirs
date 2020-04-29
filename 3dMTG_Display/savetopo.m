function savetopo(filename,vColor,format)
%Permet de sauvagarder les topos dans le format désirés
%.vcolor format matlab et .img par default. 
%Evantuellement il sera plus facile d'intégrer la compatibilité avec
%d'autre logiciel qui utilise d'autre format de donnée topographique. 
if format == 1
    [dir,fil,ext] = fileparts(filename);
    fname = fullfile(dir,[fil,'.img']);
    mat = [1 0 0 1
        0 1 0 1
        0 0 1 1
        0 0 0 1];
    rest = mod(numel(vColor),200);
    line = floor(numel(vColor)/200);
    if rest==0
        line = line; %On ajoute rien
    else
        padding = zeros(200-rest,1);
        vColor = [vColor;padding];
        line = line+1;
    end

    vColormat  = reshape(vColor,200,line);
    
    dim = [200,line,1];
    dt = [16,0]; %'float32'
    descrip = ['spm compatible vColor_',num2str(rest)];
    pinfo = [1,0,0]';
    Vn = struct('fname',fname,'mat',mat,'dim',dim,...
        'descrip',descrip,'pinfo',pinfo,'dt',dt);
    try
    Vn = spm_create_vol(Vn);    
    spm_write_plane(Vn,vColormat,1);
    catch
    fname = fullfile(dir,[fil,'.vcolor']);
    save(fname, 'vColor','-mat');
    end
else
    [dir,fil,ext] = fileparts(filename);
    fname = fullfile(dir,[fil,'.vcolor']);
    save(fname, 'vColor','-mat');
end
