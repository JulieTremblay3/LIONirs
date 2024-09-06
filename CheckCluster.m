function CheckCluster()




lstwavelet = {'C:\data\Data_NIRS\BebeResting\WAVMAT\WAV_C07_001_HBO.mat',
    'C:\data\Data_NIRS\BebeResting\WAVMAT\WAV_C11_001_HBO.mat',
    'C:\data\Data_NIRS\BebeResting\WAVMAT\WAV_C15_001_HBO.mat',
    'C:\data\Data_NIRS\BebeResting\WAVMAT\WAV_C16_001_HBO.mat'
    'C:\data\Data_NIRS\BebeResting\WAVMAT\WAV_C19_001_HBO.mat'
    'C:\data\Data_NIRS\BebeResting\WAVMAT\WAV_C21_001_HBO.mat'
    'C:\data\Data_NIRS\BebeResting\WAVMAT\WAV_C26_001_HBO.mat'
    'C:\data\Data_NIRS\BebeResting\WAVMAT\WAV_C28_001_HBO.mat'
    'C:\data\Data_NIRS\BebeResting\WAVMAT\WAV_C30_001_HBO.mat'};
for ifile = 1
    filewavelet = lstwavelet{ifile};
    twin= 0; %set 0 for all
    fwin = 1:114;
    tic
    load(filewavelet )
    
    
    
    
    [MATCORR,Args] = RUNCOH(filewavelet, twin,fwin);
    
    
    [dir1,file1,ext1]=fileparts(filewavelet);
    fileCOH = fullfile(dir1,['COH',file1(4:end),ext1]);
    fileCluster = fullfile(dir1,['CLU',file1(4:end),ext1]);
    Args.nbcluster = 10;
    Args.twin = twin;
    Args.fwin = fwin;
    C = clusterMATCORR(MATCORR,Args.nbcluster);

    save(fileCOH,'MATCORR','C','Args','-v7.3')
    save(fileCluster,'C','Args')
    saveas(gcf,[fileCluster,'.fig'],'fig')
end
%  %
%         figure
%         plot(distsum)
%         xlabel('Nb Cluster')
%         ylabel('Error')
%         %
%         fs = 1/Args.dt;
%         t = 1/fs:1/fs:(1/fs*size(MATCORR,2))
%         f = 1./Args.period
%                 
%         figure;imagesc( t,log2(Args.period),MATCORR(:,:,20))
%                             set(gca,'YLim',log2([min(period),max(period)]), ...
%                             'YDir','reverse', ...
%                             'YTick',log2(Yticks(:)), ...
%                             'YTickLabel',num2str(Fticks'), ...
%                             'layer','top')
%                             ylabel('Frequency')
%                             caxis([0 1])
%         imagesc(MATCORR(:,:,20))
%          axis xy
%         for icluster = 2:10
%             subplot(2,5,icluster)
%             IDtF =reshape(C{icluster}.IDX,size(MATCORR,1),size(MATCORR,2));
%             imagesc(t,f,IDtF);
%             axis xy
%         end
% %save(fullfile(pathout,[filloutput,'_HBR','_COH FFT','.mat']),'ZoneList','matcorr','meancorr');
% toc
%         id = 1;
%         NC=       68
%         nele = 68
%         matid = zeros(nele,nele);
%         for ielex=2:nele
%             ielex
%             ieley = 1;
%             while ieley < ielex %| ieley==1
%                 %                     option{id}.ele1 =  listelectrode{ielex};
%                 %                     option{id}.ele2 =listelectrode{ieley};
% 
%                 option{id}.matposition = [ielex,ieley];
%                 matid(ielex,ieley)=id;
%                 ieley = ieley + 1;
%                 id = id + 1;
% 
%             end
%         end
%         idhalf = find(matid); %index dans la matrice 99x99
%         idoption = matid(idhalf); %index dans la matrice option
% 
%         figure
%         for icl = 1:4
%             subplot(2,2,icl)
%             matgr1 = zeros(nele,nele);
%             matgr1(idhalf)=C(icl,:)
%             matgr1 = matgr1 +flipud(rot90(matgr1))
%             imagesc(matgr1);
%         end
%     end
%         

end

function [MATCORR,Args] = RUNCOH(filewaveletdecompositon, twin, fwin)
load(filewaveletdecompositon)
 sinv=1./(Args.scale');     
 if twin(1) ==0
     twin = 1:size(TFR,2);
 end
 id = 1;
 for i=2:size(TFR,3)
     j = 1;
     while j<i
         TFRx =  TFR(fwin,twin,i);
         TFRy =  TFR(fwin,twin,j);
         sX =  sTFR(fwin,twin,i);
         sY =  sTFR(fwin,twin,j);
         Wxy= TFR(:,:,i).*conj(TFR(:,:,j));       
         sWxyt=smoothwavelet(sinv(:,ones(1,size(Wxy,2))).* Wxy,Args.dt,Args.period,Args.Dj,Args.scale);
         sWxy = sWxyt(fwin,twin);
        % Rsq=abs(sWxy).^2./(sX.*sY);
         MATCORR(:,:,id) =abs(sWxy).^2./(sX.*sY);
         j = j + 1;
         id = id +1;
     end
 end   
    
  

end
function  C = clusterMATCORR(MATCORR,ncluster)
        %%%%%%%%%%%%% MATCORR %%%%%
        %Try to find several cluster of connectivity matrix in the time
        %frequency space
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        X = reshape(MATCORR, size(MATCORR,1)*size(MATCORR,2),size(MATCORR,3));
        clear MATCORR
        for icluster= 1:ncluster
            icluster
            [C{icluster}.IDX,C{icluster}.mean,sumd,D] = kmeans(abs(X),icluster);
            distsum(icluster) = sum(sumd);
            C{icluster}.sumd = sumd;
        end
         h = figure;
        plot(distsum)
        xlabel('Nb Cluster')
        ylabel('Error')
    end