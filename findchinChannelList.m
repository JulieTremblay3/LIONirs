function  [ listHBOch, listHBRch, listnameHbO, listnameHbR, zonelist]= findchinChannelList(NIRS, ChannelListfile, listgood)
       
        ML_new= [NIRS.Cf.H.C.id(2:3,listgood)',...
        ones(numel(listgood),1),...
        NIRS.Cf.H.C.wl(listgood)'];
        fid = fopen(ChannelListfile);
        try
            chlist = textscan(fid, '%s%s');
        catch
            disp(['Error Channel list : ', ChannelListfile,' could not be open'])
        end
        fclose(fid);
        DetL= chlist{1};
        SrsL= chlist{2};
        name =  DetL{1};
        
        if numel(name)>1
            if sum(cell2mat(strfind(DetL, 'D0'))) %imagen do not exporte as D0 use directly D1... to make the difference between systems
                Devicename = 'NIRx';
            else
                Devicename  = 'ISS Imagent';
            end
        else
            Devicename  = 'ISS Imagent';
        end
        
        list = zeros(numel(chlist{1}),1);
        for i=1:numel(chlist{1})
            if strcmp(Devicename,'NIRx')
                SDdetL = StrBoxy2SDDet(DetL{i});
                 tmp = SrsL{i};
                 SDsrsL = str2num(tmp(2:end));
            elseif strcmp(Devicename,'nirs')
                 SDdetL = StrBoxy2SDDet(DetL{i});
                 tmp = SrsL{i};
                 SDsrsL = str2num(tmp(2:end));
            else %ISS
                SDdetL = StrBoxy2SDDet_ISS(DetL{i});
                SDsrsL = StrBoxy2SDPairs(SrsL{i});
            end
            zonelist{i,1} = [DetL{i} ' ' SrsL{i}];
            L1 = find(ML_new(:,1)==SDsrsL & ML_new(:,2)==SDdetL & ML_new(:,4)==1);
            L2 = find(ML_new(:,1)==SDsrsL & ML_new(:,2)==SDdetL & ML_new(:,4)==2);
%             L1 = find(ML(:,1)==SDsrsL & ML(:,2)==SDdetL & ML(:,4)==1);
%             L2 = find(ML(:,1)==SDsrsL & ML(:,2)==SDdetL & ML(:,4)==2);
            if isempty(L1)
                listHBOch(i,1)= nan;
                listHBRch(i,1)= nan;
                sprintf(['check ', DetL{i},' ', SrsL{i}]);
                listname{i,1} = [DetL{i} ' ' SrsL{i}];
                listnameHbO{i,1}= [DetL{i} ' ' SrsL{i},'HbO','NO FOUND'];
                listnameHbR{i,1}= [DetL{i} ' ' SrsL{i},'HbR','NO FOUND'];
            else
                listHBOch(i,1)= L1(1);
                listHBRch(i,1)= L2(1);
                listname{i,1} = [DetL{i} ' ' SrsL{i}];
                listnameHbO{i,1} = [DetL{i} ' ' SrsL{i},'HbO'];
                listnameHbR{i,1} = [DetL{i} ' ' SrsL{i},'HbR'];
            end
           
                     
        end