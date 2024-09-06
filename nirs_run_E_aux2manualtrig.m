function out = nirs_run_E_aux2manualtrig(job)
%Use iss auxiliairy trigger to create manual trig file.

load(job.NIRSmat{1,1});
triglist = []; %list de tout les trigs utilisés dans l'expé
lst = length(NIRS.Dt.fir.pp);
rDtp = NIRS.Dt.fir.pp(lst).p;
name = job.E_aux2manualname;
for ifile=1:length(rDtp)
    aux5 = NIRS.Dt.fir.aux5{ifile};
    [dir1,fil1,ext1] = fileparts(rDtp{ifile});
    sortaux5 = sortrows(aux5);
  for itrig=1:size(aux5,1)
      trig = sortaux5(itrig,1);
        if isempty(find(triglist==trig))
            triglist = ([triglist,trig]);
            fid = fopen(fullfile(dir1,[name,'_', num2str(trig),'trig.m']),'w');     
            fprintf(fid,'%s\n', ['Trigvalue = ',num2str(trig),';']);
            id = numel(triglist);
            for i=1:length(rDtp)
                timingfile(id).trig{i}=[];
            end
            fclose(fid); %close previous
        else
            id = find(triglist==trig);
            %fid = fopen(fullfile(dir1,[name,'_', num2str(trig),'trig.m']),'a');     

        end
            
    
                 
         timingfile(id).trig{ifile}=[timingfile(id).trig{ifile},sortaux5(itrig,2)];
         timingfile(id).label = num2str(trig);
  end %fin aux
 
end
freq = NIRS.Cf.dev.fs;
if 1% job.E_aux2manualoption==1 %Separate file
    for idtrig = 1:numel(timingfile)
        fid = fopen(fullfile(dir1,[name,'_', timingfile(idtrig).label ,'trig.m']),'a');     
        for ifile = 1:numel(timingfile(idtrig).trig)
            [dir1,fil1,ext1] = fileparts(rDtp{ifile});
            fprintf(fid,'%s\n',['filename{',num2str(ifile),'} = [''',fil1,'''];' ]);
            fprintf(fid,'%s', ['timingfile{',num2str(ifile), '} = [']);
            list = timingfile(idtrig).trig{ifile};
            for itrig = 1:(numel(list)-1)
                fprintf(fid,'%s,',num2str(list(itrig)*1/freq))
            end
            if isempty(list)
               fprintf(fid,'%s\n',['];']);
            else
                fprintf(fid,'%s\n',[num2str(list(end)*1/freq),'];']);
            end
        end
        fclose(fid)
    end
elseif 0 %OBSOLETE job.E_aux2manualoption==2 %Concatenate        
%         for idtrig = 1:numel(timingfile)
%         fid = fopen(fullfile(dir1,[name,'_', timingfile(idtrig).label ,'trig.m']),'a');
%         fprintf(fid,'%s\n%s\n%s\n',['% Be carefull trig time have been modify to be use when the file are Concatenate one after another'],...
%             '%trig file 1 = trig file 1,', '%trig file 2 = trig file 2  + (total time of file 1)...' );
%         offset = 0
%         for ifile = 1:numel(timingfile(idtrig).trig)
%             [dir1,fil1,ext1] = fileparts(rDtp{ifile});
%             fprintf(fid,'%s\n',['filename{',num2str(ifile),'} = [{''',fil1,'''}];' ]);
%             fprintf(fid,'%s', ['timingfile{',num2str(ifile), '} = [']);
%             list = timingfile(idtrig).trig{ifile}
%             if ifile>1
%                 offset = NIRS.Dt.fir.sizebloc{ifile-1}
%             end
%             for itrig = 1:(numel(list)-1)
%                 fprintf(fid,'%s,',num2str((list(itrig)+offset)*1/freq))
%             end
%             if isempty(list)
%                fprintf(fid,'%s\n',['];']);
%             else
%                 fprintf(fid,'%s\n',[num2str((list(end)+offset)*1/freq),'];']);
%             end
%         end
%         fclose(fid)
%     end
    
end
out.NIRSmat = job.NIRSmat;
