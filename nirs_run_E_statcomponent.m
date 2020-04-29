function out = nirs_run_E_statcomponent(job)

if isfield(job.c_statcomponent,'b_TtestOneSample')
    AllC = [];
    for i=1:numel(job.c_statcomponent.b_TtestOneSample.f_component)
        load(job.c_statcomponent.b_TtestOneSample.f_component{i},'-mat');
        [dir1,file1,ext1]=fileparts(job.c_statcomponent.b_TtestOneSample.f_component{i})
        AllC = [AllC,A];
    end
    A = A';

     for ich=1:size(AllC,1)
       % [h,p,ci,stats] = ttest(AllC(ich,:)); %do not suport nan 
        tval(ich) =  nanmean(AllC(ich,:))./        (nanstd(AllC(ich,:)./sqrt(sum(~isnan(AllC(ich,:))))));
        df = sum(~isnan(AllC(ich,:)))-1
        if job.c_statcomponent.b_TtestOneSample.m_TtestOneSample == 1
            pval(ich) = 2 * tcdf(-abs(tval(ich)), df);
            option = 'twotail'
        elseif job.c_statcomponent.b_TtestOneSample.m_TtestOneSample == 2
             pval(ich) = tcdf(+tval(ich), df);   
            option = 'lefttail'
        elseif job.c_statcomponent.b_TtestOneSample.m_TtestOneSample == 3
            pval(ich) = tcdf(-tval(ich), df);
            option = 'righttail'
        end
        mval(ich) = nanmean(AllC(ich,:));
     end
     A = mval;
     save(fullfile(dir1,['ONESAMPLE_Mean n=',num2str(df+1),'.mat']),'A','zonelist','option')     
     A = tval;
     save(fullfile(dir1,['ONESAMPLE_Tmap n=',num2str(df+1),'.mat']),'A','zonelist','option')

     A = mval.*double(pval<0.05);
     save(fullfile(dir1,['ONESAMPLE__mean05unc.mat']), 'A','zonelist','option')
     A = mval.*double(pval<0.01);
     save(fullfile(dir1,['ONESAMPLE__mean01unc.mat']),'A','zonelist','option')
     A = mval.*double(pval<0.001);
     save(fullfile(dir1,['ONESAMPLE__mean001unc.mat']),'A','zonelist','option')
    
     
     [FDR,Q] = mafdr(pval)         
     A = mval.*double(Q<0.05);
     save(fullfile(dir1,['ONESAMPLE__mean05fdr.mat']),'A','zonelist','option')
    [FDR,Q] = mafdr(pval)
     A = mval.*double(Q<0.01);
     save(fullfile(dir1,['ONESAMPLE__mean01fdr.mat']),'A','zonelist','option')
    [FDR,Q] = mafdr(pval)
     A = mval.*double(Q<0.001);
     save(fullfile(dir1,['ONESAMPLE__mean001fdr.mat']),'A','zonelist','option')

     
elseif isfield(job.c_statcomponent,'b_TtestUnpaired')
    
    
    
        if job.c_statcomponent.b_TtestUnpaired.m_TtestOneSample == 1
            option = 'twotail'
        elseif job.c_statcomponent.b_TtestUnpaired.m_TtestOneSample == 2  
            option = 'lefttail'
        elseif job.c_statcomponent.b_TtestUnpaired.m_TtestOneSample == 3
            option = 'righttail'
        end
    
      AllG1 = [];
    for i=1:numel(job.c_statcomponent.b_TtestUnpaired.f_componentG1)
        load(job.c_statcomponent.b_TtestUnpaired.f_componentG1{i},'-mat');
        [dir1,file1,ext1]=fileparts(job.c_statcomponent.b_TtestUnpaired.f_componentG1{i})
        AllG1 = [AllG1,A];
    end 
    
     AllG2 = [];
    for i=1:numel(job.c_statcomponent.b_TtestUnpaired.f_componentG2)
        load(job.c_statcomponent.b_TtestUnpaired.f_componentG2{i},'-mat');
        [dir1,file1,ext1]=fileparts(job.c_statcomponent.b_TtestUnpaired.f_componentG2{i})
        AllG2 = [AllG2,A];
    end
    
     for ich=1:size(AllG1,1)
        stats = testt(AllG1(ich,:),AllG2(ich,:)); %do not suport nan 
        tval(ich) =   stats.tvalue;
        df = stats.tdf
        pval(ich)=stats.tpvalue;      
        mval(ich) =  nanmean(AllG1(ich,:)) - nanmean(AllG2(ich,:));
     end
     A = nanmean(AllG1);
     save(fullfile(dir1,['TWOSAMPLE_Mean G1 n=',num2str(df+1),'.mat']),'A','zonelist','option')
     A = nanmean(AllG2);
     save(fullfile(dir1,['TWOSAMPLE_Mean G2 n=',num2str(df+1),'.mat']),'A','zonelist','option')
     
     
     A = mval;
     save(fullfile(dir1,['TWOSAMPLE_Mean G1-G2 n=',num2str(df+1),'.mat']),'A','zonelist','option')

     A = tval;
     save(fullfile(dir1,['TWOSAMPLE_Tmap n=',num2str(df+1),'.mat']),'A','zonelist','option')
     
     A = tval.*double(pval<0.05);
     save(fullfile(dir1,['TWOSAMPLE_Tmap_05unc.mat']), 'A','zonelist','option')
     A = mval.*double(pval<0.05);
     save(fullfile(dir1,['TWOSAMPLE__mean05unc.mat']),'A','zonelist','option')
     
     A = mval.*double(pval<0.01);
     save(fullfile(dir1,['TWOSAMPLE__mean01unc.mat']),'A','zonelist','option')
     A = mval.*double(pval<0.001);
     save(fullfile(dir1,['TWOSAMPLE__mean001unc.mat']),'A','zonelist','option')
    
     
     [FDR,Q] = mafdr(pval)         
     A = mval.*double(Q<0.05);
     save(fullfile(dir1,['TWOSAMPLE__mean05fdr.mat']),'A','zonelist','option')
    [FDR,Q] = mafdr(pval)
     A = mval.*double(Q<0.01);
     save(fullfile(dir1,['TWOSAMPLE__mean01fdr.mat']),'A','zonelist','option')
    [FDR,Q] = mafdr(pval)
     A = mval.*double(Q<0.001);
     save(fullfile(dir1,['TWOSAMPLE_mean001fdr.mat']),'A','zonelist','option')
     
     

end


% G1 = load('C:\data\Analysed\fcNIRSepilepsie\P25\GLMexport\Sem\D1matrix allGLMGLMP25Global.zoneHbO033.mat')
% G2 = load('C:\data\Analysed\fcNIRSepilepsie\P25\GLMexport\Sym\D1matrix allGLMGLMP25Global.zoneHbO033.mat')
% 
% for ich = 1:size(G1.A,1)
% [h,p,ci,stats] =  ttest(G1.A(ich,:)',G2.A(ich,:)')
% 
% tval(ich,1)  = stats.tstat;
% pval(ich,1) = p;
% [p,h,stats] = ranksum(G1.A(ich,:)',G2.A(ich,:)')
% pranksum(ich,1) = p
% end
% 
% A = tval
% save('TTESTPAIRED','A')
% 
% A = nanmean(G1.A,2)- nanmean(G2.A,2)
% save('MeanG1moinsG2','A')
% A = A.*(pval<0.1)
% save('MeanG1moinsG2pval','A')
% 
% A = nanmean(G2.A,2)- nanmean(G1.A,2)
% save('MeanG2moinsG1','A')

try
    out = {srsfile};
catch
    out = 'no dependancy';
end