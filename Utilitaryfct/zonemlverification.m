function [zone, varargout]=zonemlverification(zone,ML)
%V�rification des champs zone.plot (combinaison source detecteur des canaux de la zone)
%et zone.plotLst (Indice ml) de chaque zone, si jamais il y a des erreurs
%N.B. le m�me sujet peux avoir un ml diff�rent si les distances de pruning
%utilis� lors de la conversion des donn�es �tait diff�rent 
%[zone] = zonemlverification(zone,ml) %retoure les zones avec correction du
%plotLst dans les donn�es actuelles s'il y a lieu
%[zone,flag]= zonemlverification(zone,ml) retourne les zones avec correction du
%plotLst dans les donn�es actuelles s'il y a lieu et une valeur 0 si aucune
%modification n'a �t� effectu� 1 si une modification a �t� effectu�.
flag = 0;
nout = max(nargout,1) - 1;
 for i = 1:numel(zone.plotLst)
            plotLst =  zone.plotLst{i};
            plotLstnew = [];
            plotold = zone.plot{i};
            plotnew = [];
            for indplot = 1:numel(plotLst)
                srs = plotold(indplot,1); %srs
                det = plotold(indplot,2); %det                 
                newid = find(ML(:,1) == srs & ML(:,2) == det );
                if ~isempty(newid)
                    plotLstnew=[plotLstnew;newid(1)];
                    plotnew = [plotnew;plotold(indplot,:)];
                end
            end
            zone.plotLst{i} = plotLstnew;
            zone.plot{i}=plotnew;
            if numel(plotLstnew)==numel(plotLst)
            else
                  flag = 1
            end
 end
for k = 1:nout
   varargout{k} = flag;
end
        
       