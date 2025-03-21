function id = Col_ID(alltitlecolumn, columnlabeltofind)
    for idcol=1:numel(alltitlecolumn)
        if strcmp(deblank(alltitlecolumn{idcol}),deblank(columnlabeltofind))
            id = idcol;
            break
        end
    end 