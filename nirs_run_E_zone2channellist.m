function nirs_run_E_zone2channellist(job)
for f=1:numel(job.zonemat)
    load(job.zonemat{f},'-mat');
    [filepath,name,ext] = fileparts(job.zonemat{f});
    fid = fopen(fullfile(filepath,[name,'.txt']),'w');
    disp(['Create file:',fullfile(filepath,[name,'.txt'])])
    for izone = 1:numel(zone.label)
        fprintf(fid,'Label: %s\n',zone.label{izone});
        fprintf(fid,'RGBcolor: %s\n',num2str(zone.color(izone,:)));
        plotLst = zone.plotLst{izone};
        for iplotlst = 1:numel(plotLst);
            ich = plotLst(iplotlst);
            ML_new = zone.ml;
            switch job.m_zonematformat
                case 2 %'ISS Imagent'
                    strDet = SDDet2strboxy_ISS(ML_new( ich,2));
                    strSrs = SDPairs2strboxy_ISS(ML_new( ich,1));
                    fprintf(fid,'%s\t%s\n',strDet,strSrs);
                case 1 %'NIRx'                     
                    strDet = SDDet2strboxy(ML_new( ich,2));
                    strSrs = SDPairs2strboxy(ML_new( ich,1));
                    fprintf(fid,'%s\t%s\n',strDet,strSrs);                   
            end
        end
    end
    fclose(fid);
end