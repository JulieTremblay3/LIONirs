function TrigAux2ManualTrig(Auxtrig, filename, fq)
% Create .m manual trig file to let the user modify or check the trig if necessary
% par convention trig.m
[pathstr, name, ext] = fileparts(filename{1})


%trouver le nombre de fichier de trig maximal à creer. 
trigall  = [];
for ifile=1:numel(Auxtrig)
 trigall = [trigall,Auxtrig{ifile}];
end
sorttrig =  sort(trigall(:,1));
trigid = (sorttrig(1)); 
val = find((sorttrig(2:end) - sorttrig(1:end-1)));
if ~isempty(val) 
    % A gerer plus d'une sorte de trig)
    for i=1:numel(val)
        trigid = [trigid, val(i)];
    end
 end
%Tout les trigs différents sont trouvés dans trigid

for itrig = 1:numel(trigid) %Tout les trigs différents dans des fichiers séparés
    fid=fopen([pathstr,,filesep,'ONSET' ,num2str(trigid(itrig)), 'trig.m'],'w')

    
    for ifile=1:numel(filename)
            triglst = Auxtrig{ifile};
            indicetrig=find(trigid(itrig)==triglst(:,1))
            [pathstr, name, ext] = fileparts(filename{1})
            fprintf(fid,'%s%s%s%s%s\n','filename{',num2str(ifile),'}=''',name,''';' )  %filename{1} = 'Acq_1402191153_Ali_Visual_006'
            fprintf(fid,'%s%s\n', 'Trigvalue = ', num2str(trigid(itrig)));
            fprintf(fid,'%s%s%s', 'timingfile{', num2str(ifile),'}=[')
            for isample=1:numel(indicetrig)
                fprintf(fid,'%s',num2str(triglst(indicetrig(isample),2)*1/fq));              
                if isample==numel(indicetrig)
                    fprintf(fid,'%s', ']')
                else
                    fprintf(fid,'%s\n', '');  
                end
            end
           
          %col 1 id trig, col 2 sample point 
            
    end
    
    fclose(fid)

end
