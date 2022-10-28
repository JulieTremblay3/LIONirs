function [h,atab,ctab,stats]=aoctoolpermutation(x,y,group,  nperm)
%Apply aoctool statistic and built a null empirical distribution permuting
%group. 
% aoctool(x,y,group)see matlab help
% return one additionnal column in the atab using permutation distribution instead of normal distribution. 
for i=1:nperm
    % when apply on many link use the highest value of all comparison in
    % this case all link will be control and the maximal f value use for the
    % empirical distribution. Keep here the maximal value for the empirical distribution.  
    % Theroriticaly best control, If the maximum is too strick you can use the maximal 
    % Ref Lage-Castellanos, A., Martínez-Montes, E., Hernández-Cabrera, J.A., Galán, L., 2010. False discovery rate and permutation test: An evaluation in ERP data analysis. Statistics in Medicine 29, 63–74. https://doi.org/10.1002/sim.3784

    idperm=randperm(numel(group));
    [h,atab,ctab,stats] = aoctool(x(idperm),y,group(idperm),0.05,'x','y','Groupe','off');
    fdist_groupe(i) = atab{2,5};
    fdist_cov(i) = atab{3,5};
    fdist_GroupebyCoV(i) = atab{4,5};
end


%reference real groupe
[h,atab,ctab,stats] = aoctool(x,y,group,0.05,'x','y','Groupe','off');
f_groupe = atab{2,5};
f_cov = atab{3,5};
f_GroupebyCoV = atab{4,5};

[fecdf,xecdf] = ecdf( [fdist_groupe,fdist_cov,fdist_GroupebyCoV]);
figure
plot(xecdf,fecdf)
atabperm{1,1} = 'Prob>F permutation'; 
atabperm{2,1} = 1-fecdf(sum(xecdf<=f_groupe));
atabperm{3,1} = 1-fecdf(sum(xecdf<=f_cov));
atabperm{4,1} =  1-fecdf(sum(xecdf<= f_GroupebyCoV));
atabperm{5,1} =  [];
atab = [atab,atabperm];

