function out = nirs_run_createseedlist(job)

for filenb=1:size(job.NIRSmat,1) %JT
    NIRS = [];
    load(job.NIRSmat{filenb,1});
    
    %Modification pour chaque zone du plotLst en fonction des canaux de plot
    ML_new= [NIRS.Cf.H.C.id(2:3,:)',...
        ones(size(NIRS.Cf.H.C.id,2),1),...
        [ones(size(NIRS.Cf.H.C.id,2)/2,1);ones(size(NIRS.Cf.H.C.id,2)./2,1).*2]];

        [path, fname, extension]=fileparts(job.NIRSmat{filenb,1});
       % [file, path] = uiputfile(fullfile(path,'channellist.txt'));
        filenameChannelList = fullfile(path,'channellist.txt');
      % ChannelLabels = ConvertmlIDsrs2label(NIRS)        
    %   NIRS.Cf.dev.n = 'NIRx'
        switch NIRS.Cf.dev.n
            case 'ISS Imagent'
                fid = fopen(filenameChannelList,'w');               
                for i = 1:size(ML_new,1)/2   
                    strDet = SDDet2strboxy_ISS(ML_new(i,2));
                    strSrs = SDPairs2strboxy_ISS(ML_new(i,1));
                    fprintf(fid,'%s\t%s\n',strDet,strSrs);
                end
                fclose(fid)
            case 'NIRx'
                 fid = fopen(filenameChannelList,'w');               
                    for i = 1:size(ML_new,1)/2   
                    strDet = SDDet2strboxy(ML_new(i,2));
                    strSrs = SDPairs2strboxy(ML_new(i,1));
                    fprintf(fid,'%s\t%s\n',strDet,strSrs);
                    end
                    fclose(fid);
            case 'NIRSx'
                 fid = fopen(filenameChannelList,'w');               
                    for i = 1:size(ML_new,1)/2   
                    strDet = SDDet2strboxy(ML_new(i,2));
                    strSrs = SDPairs2strboxy(ML_new(i,1));
                    fprintf(fid,'%s\t%s\n',strDet,strSrs);
                    end
                    fclose(fid);
            case 'NIRS FILE HOMER'
                    fid = fopen(filenameChannelList,'w');               
                    for i = 1:size(ML_new,1)/2   
                    strDet = SDDet2strboxy(ML_new(i,2));
                    strSrs = SDPairs2strboxy(ML_new(i,1));
                    fprintf(fid,'%s\t%s\n',strDet,strSrs);
                    end
                    fclose(fid);
             otherwise
                    fid = fopen(filenameChannelList,'w');               
                    for i = 1:size(ML_new,1)/2   
                    strDet = SDDet2strboxy(ML_new(i,2));
                    strSrs = SDPairs2strboxy(ML_new(i,1));
                    fprintf(fid,'%s\t%s\n',strDet,strSrs);
                    end
                    fclose(fid);
        end
    disp(['Channel list save: ', filenameChannelList])
end
out.NIRSmat = job.NIRSmat;
