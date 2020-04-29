function out = nirs_run_E_statmatrix(job)

if isfield(job.c_statmatrix,'b_TtestOneSamplematrix')
    AllC = [];
    xlslistfile = job.c_statmatrix.b_TtestOneSamplematrix.f_matrix{1};
    [raw, txt, info]=xlsread([xlslistfile]);
    id =1
    for isubject=2:size(info,1)
        if info{ isubject,4}==1
            MAT = load(fullfile(info{isubject,1},[ info{isubject,2},'.mat']));
            MATall( id ,:,:)=MAT.meancorr;
            id = id + 1;            
        end
    end
    meanall = squeeze(nanmean(MATall,1));
   refval  = 0.1; 
  
   tval = squeeze((nanmean(MATall(:,:,:),1)-refval)./(nanstd(MATall(:,:,:),1)./sqrt(sum(~isnan(MATall(:,:,:)),1))  ));
   dfall = squeeze(sum(~isnan(MATall(:,:,:)),1))-1;
   for i=1:size(tval,1)
       for j=1:size(tval,2) 
           try
        pval(i,j) = tcdf(-tval(i,j),  dfall(i,j));
           catch
                pval(i,j) = 1;
           end   

       end
   end
   infonew = [{'Dir'},{'File'},{'Zone'},{'GR'}];
   file = 'tval'
   ZoneList = MAT.ZoneList;
   matcorr =  tval;
   meancorr = tval;
  totaltrialgood = mean(dfall(:))
  save(fullfile(info{isubject,1},'tval.mat'),'ZoneList','matcorr','meancorr','totaltrialgood');
   new = [{info{isubject,1}},{file}, {info{isubject,3}},{1} ];
  infonew = [infonew;new];
  
   file = 'tvalp05'
   matcorr = tval.*double(pval<0.05);
   meancorr = tval.*double(pval<0.05);
    save(fullfile(info{isubject,1},[file,'.mat']),'ZoneList','matcorr','meancorr','totaltrialgood');     
   new = [{info{isubject,1}},{file}, {info{isubject,3}},{1} ];
  infonew = [infonew;new];
   [FDR,Q] = mafdr(pval(:));
    Q = reshape(Q,size(pval));
    infonew = [infonew;new];
    
      file = 'meanp05'
   matcorr = meanall.*double(pval<0.05);
   meancorr = meanall.*double(pval<0.05);
    save(fullfile(info{isubject,1},[file,'.mat']),'ZoneList','matcorr','meancorr','totaltrialgood');     
   new = [{info{isubject,1}},{file}, {info{isubject,3}},{1} ];
    infonew = [infonew;new];
    
   file = 'tvalFDR05'
   matcorr = tval.*double(Q<0.05);
   meancorr = tval.*double(Q<0.05);
    save(fullfile(info{isubject,1},[file,'.mat']),'ZoneList','matcorr','meancorr','totaltrialgood');     
   new = [{info{isubject,1}},{file}, {info{isubject,3}},{1} ];
  infonew = [infonew;new];

  
   [FDR,Q] = mafdr(pval(:));
    Q = reshape(Q,size(pval));
   file = 'tvalFDR01'
   matcorr = tval.*double(Q<0.01);
   meancorr = tval.*double(Q<0.01);
    save(fullfile(info{isubject,1},[file,'.mat']),'ZoneList','matcorr','meancorr','totaltrialgood');     
   new = [{info{isubject,1}},{file}, {info{isubject,3}},{1} ];
  infonew = [infonew;new];

  
[filepath,name,ext] = fileparts(xlslistfile)
xlswrite(fullfile(filepath,[name,'SimpleTtest.xlsx']),infonew);
out='test';
end

