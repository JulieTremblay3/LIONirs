path = 'D:\Fluence\Review\VIDEOADULTE\'

namesubject = { '33_FVE019'
     '34_FVE25' 
     '35_FVE21'
     '36_FVE20'
     '37_FVE29'
     '38_FVE310'
     '40_FVE312'
     '42_FPC01_OB'
     '43_FPC02'};

% namesubject = { '31_FVE27'
%      '32_FVE10' 
%      '35_FVE21'
%      '36_FVE20'
%      '37_FVE29'
%      '38_pFVE310'
%      '39_FVE311'
%      '40_FVE312'};


%  '19_FVE32',
%     '20_FVE36'
%     '21_FVE22'
%     '22_FVE15'
%     '23_FVE16'
%     '24_FVE13'
%     '25_FVE31'
%      '26_FVE18'
%      '27_FVE09'
%      '28_FVE11'
%      '29_FVE12'
%      '30_FVE08'
for itime = 1:40
    all = [];
    labeltime =  fixdecimal2string(itime,4,2);
    %labeltime =num2str(itime);
    for iname=1:numel(namesubject);
            load([path,namesubject{iname,1},'.aviHbOTime   ', labeltime,'(s).vcolor'],'-mat');
            %ind = find(vColor==0);
            %vColor(ind)=nan;
            all = [all,vColor];
    end
    vColor = nanmean(all')';
    %id = find(isnan(vColor));
    %vColor(id)=0;
    [nameout]= [path,'HbOAvg',labeltime,'(s).vColor']
    
     save([nameout],'vColor','-mat')
%      for iname=1:numel(namesubject);
%             load([path,namesubject{iname,1},'.aviHbRTime   ', num2str(itime),'(s).vcolor'],'-mat');
%             all = [all,vColor];
%     end
%     vColor = mean(all')';
%     [nameout]= [path,'HbRAvg',num2str(itime),'(s).vColor']
%     save([nameout],'vColor','-mat')
end


