function out = nirs_run_E_channellist2zone(job)
    load(job.NIRSmat{1});
    ml  = [NIRS.Cf.H.C.id(2:3,:)',...
        ones(size(NIRS.Cf.H.C.id,2),1),...
        [ones(size(NIRS.Cf.H.C.id,2)/2,1);ones(size(NIRS.Cf.H.C.id,2)./2,1).*2]];
    SD.SrcPos = NIRS.Cf.H.S.r.o.mm.p';
    SD.DetPos = NIRS.Cf.H.D.r.o.mm.p';
    SD.Lambda = NIRS.Cf.dev.wl;
    SD.nSrcs  = numel(SD.SrcPos)/3;
    SD.nDets  = numel(SD.DetPos)/3;
    SD.SrcAmp = ones(SD.nSrcs ,2);
    SD.DetAmp = ones(SD.nDets,2);
    for id =1:max(size(ml))
        idsrc=ml(id,1);
        iddet=ml(id,2);
        SrcPos=SD.SrcPos(idsrc,:);
        DetPos=SD.DetPos(iddet,:);
        pos(id,1) = (DetPos(1)-SrcPos(1))/2+SrcPos(1);
        pos(id,2) = (DetPos(2)-SrcPos(2))/2+SrcPos(2);
        pos(id,3) = (DetPos(3)-SrcPos(3))/2+SrcPos(3);
        pos(id,4) = sqrt( (DetPos(1)-SrcPos(1))^2 + (DetPos(2)-SrcPos(2))^2+ (DetPos(3)-SrcPos(3))^2);           %'distance entre les canaux
    end
     zone.ml = ml;
     zone.SD = SD;
     zone.pos = pos;
    
for f=1:numel(job.zonetxt)
fid = fopen(job.zonetxt{f});      
[pathzone,zonename,ext]=fileparts(job.zonetxt{f});
chlist = textscan(fid, '%s%s%s%s');
fclose(fid);
izone =0;
%zonemodel = load('NIRS1020.zone','-mat')
    DetL= chlist{1};
    SrsL= chlist{2};
    color2 = chlist{3};
    color3 = chlist{4};
    for ilist=1:numel(DetL)
        if strcmp(DetL(ilist), 'Label:') 
            izone = izone +1;
             zone.label{1,izone} =   SrsL{ilist};
              zone.plotLst{izone} = [];
              zone.plot{izone} = [];
        elseif strcmp(DetL(ilist),'RGBcolor:')
             zone.color(izone,:) = [str2num( SrsL{ilist}), str2num(color2{ilist}), str2num(color3{ilist}) ];
        else
                switch  job.m_zonematformat
                case 2 %'ISS Imagent'
                    SDdetL = StrBoxy2SDDet_ISS(DetL{ilist});
                    SDsrsL = StrBoxy2SDPairs(SrsL{ilist});
                case 1 %'NIRx' 
                    SDdetL = StrBoxy2SDDet(DetL{ilist}); 
                     tmp = SrsL{ilist};
                     SDsrsL =str2num(tmp(2:end));
                     
                end
               zone.plotLst{izone} = [zone.plotLst{izone},find( ml(:,1)==SDsrsL & ml(:,2)==SDdetL & ml(:,4)==1 )];
               zone.plot{izone} =    [zone.plot{izone}; [SDsrsL, SDdetL] ];
    end

    end
   
    [filepath,name,ext] =  fileparts(job.NIRSmat{1});
    save(fullfile(filepath,[zonename,'.zone']),'zone');
    disp(['SAVE : ', fullfile(filepath,[zonename,'.zone'])])

end
out.NIRSmat = job.NIRSmat;
