function d1=selectimagetype(echantillon_time,type,dconc,varargin)
if isempty(varargin)
    optionGUI =1 % GUI homer 
else
    if numel(varargin)>=1
        optionGUI = varargin{1} % GUI video %3
    end
end

if optionGUI==1
    guiHOMER = getappdata(0,'guiHOMER');
elseif optionGUI==2
    guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
elseif optionGUI==3
    guiHOMER = getappdata(0,'gui_SPMvideo');
end

HOMERhandles = guihandles(guiHOMER);
PMI = get(guiHOMER,'UserData');
global currentsub
cf = PMI{currentsub}.currentFile; 
if type==1 %dConc
      d1 = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,:,dconc); 
elseif type==2 %dConc mean
    pstart = find(PMI{currentsub}.data(cf).HRF.tHRF<PMI{currentsub}.imgAvgStart);
    pstop = find(PMI{currentsub}.data(cf).HRF.tHRF<(PMI{currentsub}.imgAvgStart+PMI{currentsub}.imgAvgLen));
    echantillon_time=pstart(end):pstop(end);
    d1 = [mean(PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,:,dconc),1)]';
elseif type==3 %T
      cfA = get(HOMERhandles.Lisa_popupA,'value');
      d1 = PMI{currentsub}.data(cfA).HRF.t_bas_act(echantillon_time,:,dconc);
elseif type==4 %-T
      cfA = get(HOMERhandles.Lisa_popupA,'value');
      d1 = -PMI{currentsub}.data(cfA).HRF.t_bas_act(echantillon_time,:,dconc);      
     
%       pstart = find(PMI{currentsub}.data(cf).HRF.tHRF<5);
%      pstop = find(PMI{currentsub}.data(cf).HRF.tHRF<15);
%      echantillon_time=pstart(end):pstop(end);
%      d1 =-[mean(PMI{currentsub}.data(cf).HRF.t_bas_act(echantillon_time,:,dconc),1)]';
     
      
elseif type==5 %P
       cfA = get(HOMERhandles.Lisa_popupA,'value');
       d1 = PMI{currentsub}.data(cfA).HRF.p_bas_act(echantillon_time,:,dconc);
elseif type==6 %Beta GLM 
     num_stim = 1;
     for ch = 1:size(PMI{currentsub}.data(cf).dConc,2)    
           if isfield(PMI{currentsub}.data(cf).GML{ch,dconc},'mh')
            d1(:,ch,1) = PMI{currentsub}.data(cf).GML{ch,dconc}.mh(num_stim);
           else
               d1(:,ch,1)=0;
           end
     end
     echantillon_time = 1;
elseif type==7 % T Beta
          num_stim = 1;
         for ch = 1:size(PMI{currentsub}.data(cf).HRF.AvgC,2)    
               if isfield(PMI{currentsub}.data(cf).GML{ch,dconc},'mh')
                   d1(:,ch,1) = PMI{currentsub}.data(cf).GML{ch,dconc}.mh(num_stim)/ sqrt(PMI{currentsub}.data(cf).GML{ch,dconc}.vh(num_stim));
               else
                   d1(:,ch,1)=0;
               end
         end
         echantillon_time = 1;
%      num_stim = 1; 
%      for ch = 1:size(PMI{currentsub}.data(cf).HRF.AvgC,2)    
%          if isfield(PMI{currentsub}.data(cf).GML{ch,1},'coefcorrfreq')
%             d1(:,ch,1) = PMI{currentsub}.data(cf).GML{ch,1}.coefcorrfreq;
%          end
%      end
%      echantillon_time = 1;
elseif type==8 %User
     [name,path]= uigetfile('.mat');
     d1 = load('-mat',[path name]); 
     d1 = d1.A;
%      load('G:\VIEWLACIE_GRISJT\MC\18janMC\datadan\dan107b_BGmean1035_830.DOT','-mat')
%      d1 = Data(1:end/2);
%      % d1 = Data(end/2+1:end);
elseif type==9 %User
     [name,path]= uigetfile({'*.vcolor';'*.img'}); 
     [path1, name1, ext] = fileparts([path name])
     d1= opentopo([path name]);
%      if strcmp(upper(ext),'.VCOLOR')         
%         d1 = load('-mat',[path name]); 
%      elseif strcmp(ext,'.img')
%         Vtest = spm_vol([path name]);
%         try
%         d1.vColor = spm_slice_vol(Vtest,spm_matrix([0 0 1]),Vtest.dim(1:2),0);
%         catch
%             msgbox('File can''t not be open, verify if SPM is installed')
%         end
%          
%      end
%      d1 = d1.vColor;
elseif type==10
    d1HbO = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,:,1)-PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time-1,:,1);
    d1 = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,:,1);
    d1HbR = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,:,2)-PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time-1,:,2);
%     indBad = find(d1HbO > d1HbR );  
    indBad = find(d1HbO>0 | d1HbR<0)    
    if ~isempty(indBad)
     d1(indBad) = 0;
    end       
elseif type==11 %soustraction    
    d1 = PMI{currentsub}.dataop.dconc(echantillon_time,:,dconc);
elseif type==12 %diff
    d1 = (PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time+1,:,dconc) - PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,:,dconc))./...
        (PMI{currentsub}.data(cf).t(2) -PMI{currentsub}.data(cf).t(1)) ; 
elseif type==13 %po_Conc
    d1 = 1-PMI{currentsub}.data(cf).HRF.P_Conc(echantillon_time,:,dconc);
elseif type==14 %nbavg
    d1 = PMI{currentsub}.data(cf).HRF.navg;
elseif type==15 %Conc/stdBaseline
    d1 = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,:,dconc)./PMI{currentsub}.data(cf).HRF.STDBaseline(1,:,dconc)
    %d1 = PMI{currentsub}.data(cf).HRF.STDBaseline(1,:,dconc); 
elseif type== 16
        time = 0:20;
        threshold =5;
        answer = inputdlg('Enter time, threshold and sign > or < ','dlg_title',3,{'0:20','5','>'});
        answer=answer{1};
        time  = str2num(answer(1,:));
        threshold = str2num(answer(2,:));
        sign = answer(3,:);
        signplus = findstr(sign,'>');
        start = find(PMI{currentsub}.data(cf).HRF.tHRF<time(1));
        start=start(end);
        stop = find(PMI{currentsub}.data(cf).HRF.tHRF<time(end));
        stop = stop(end);
        echantillon_time=start:stop;
        dConc = PMI{currentsub}.data(cf).HRF.AvgC(echantillon_time,:,dconc)./(ones(numel(echantillon_time),1)*PMI{currentsub}.data(cf).HRF.STDBaseline(1,:,dconc)); ;
        idlist = [];
        for i = 1:size(dConc,2);
            if isempty(signplus)
                idt = find(dConc(:,i)<=threshold);
            else
                idt = find(dConc(:,i)>=threshold);
            end
             if isempty(idt)                
                 idlist = [idlist,i];
                 d1(i)=0;
             else
                 d1(i)= idt(1); 
             end
        end
%       d1 = log(d1);
        %figure; plot(PMI{currentsub}.data(cf).HRF.tHRF(start+d1),'x');
        d1 = PMI{currentsub}.data(cf).HRF.tHRF(start+d1);
        d1 = max(d1)-d1;
        d1(idlist)=0;
        figure
        hold on
        plot(PMI{currentsub}.data(cf).HRF.tHRF(echantillon_time), dConc)
        plot([PMI{currentsub}.data(cf).HRF.tHRF(echantillon_time(1)),...
        PMI{currentsub}.data(cf).HRF.tHRF(echantillon_time(end))],[threshold,threshold])
%         figure;plot(d1,'x')
%         A = d1;
%         [name,path] = uiputfile('.mat','Save time')
%         if name == 0
%              
%         else
%             save([path,name],'A','-mat')
%         end

        echantillon_time = 4;
elseif type== 15
end