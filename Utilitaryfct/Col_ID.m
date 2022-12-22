function id = Col_ID(alltitlecolumn, columnlabeltofind)
    for idcol=1:numel(alltitlecolumn)
        if strcmp(alltitlecolumn{idcol},columnlabeltofind)
            id = idcol;
            break
        end
    end