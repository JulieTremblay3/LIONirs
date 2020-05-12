 function [num,txt,raw] =readtxtfile_asxlsread(filename)
%Opem tab separated txt file and return input similar to xlsread function
irow = 1;
num = [];
txt = [];
raw = [];
fid = fopen(filename);
while ~feof(fid)
icol = 1;

    line = fgetl(fid);
    while ~isempty(line)
        if icol==4
            1
        end
    [token, line] = strtok(line,char(9));
    numtoken=str2num(token);
    if isempty(numtoken) | numel(numtoken)>1
        raw{irow,icol} = token;
    else
        raw{irow,icol} = str2num(token);
    end
        icol = icol + 1;
    end 
 irow =  irow +1;
end
%ensure they do not have empty space fill with nan instead
for  irow = 1:size(raw,1)
    for icol = 1:size(raw,2)
       if isempty(raw{irow ,icol})
           raw{irow ,icol} = nan;
       end
    end
end

    [num, txt] = xlsreadSplitNumericAndText(raw)
    fclose(fid)
 end


function [numericData, textData] = xlsreadSplitNumericAndText(data)
    % xlsreadSplitNumericAndText parses raw data into numeric and text arrays.
    %   [numericData, textData] = xlsreadSplitNumericAndText(DATA) takes cell
    %   array DATA from spreadsheet and returns a double array numericData and
    %   a cell string array textData.
    %
    %   See also XLSREAD, XLSWRITE, XLSFINFO.
    %   Copyright 1984-2015 The MathWorks, Inc.

    
    % ensure data is in cell array
    if ischar(data)
        data = cellstr(data);
    elseif isnumeric(data) || islogical(data)
        data = num2cell(data);
    end
    
    % Check if raw data is empty
    if isempty(data)
        % Abort when all data cells are empty.
        textData = {};
        numericData = [];
        return
    end
    
    % Initialize textData as an empty cellstr of the right size.
    textData = cell(size(data));
    textData(:) = {''};
    
    % Find non-numeric entries in data cell array
    isTextMask = cellfun('isclass',data,'char');
    
    anyText = any(isTextMask(:));
    
    % Place text cells in text array
    if anyText
        textData(isTextMask) = data(isTextMask);
    else
        textData = {};
    end
    
    % place NaN in empty numeric cells
    if anyText
        data(isTextMask)={NaN};
    end
    % Each mask we make is the same dimensions of the data. This can be
    % very big so conserve memory
    clear isTextMask;
    
    % Excel returns COM errors when it has a #N/A field.
    isErrorMask = strcmp(textData, 'ActiveX VT_ERROR: ');

    if any(isErrorMask(:))
        textData(isErrorMask) = {'#N/A'};
    end
    % Clear the mask to save memory
    clear isErrorMask;  
    
    % Trim the leading and trailing empties from textData
    emptyTextMask = cellfun('isempty', textData);
    [rowStart, colStart] = getCorner(emptyTextMask, 'first');
    [rowEnd, colEnd] = getCorner(emptyTextMask, 'last');
    textData = textData(rowStart:rowEnd, colStart:colEnd);
    
    % Clear the mask to save memory
    clear emptyTextMask;

    % Convert cell array to numeric array through concatenating columns then
    % rows.
    cols = size(data,2);
    tempDataColumnCell = cell(1,cols);
    % Concatenate each column first
    for n = 1:cols
        tempDataColumnCell{n} = cat(1, data{:,n});
    end
  
    % Now concatenate the single column of cells into a numeric array.
    numericData = cat(2, tempDataColumnCell{:});
    % Clear the mask to save memory
    clear tempDataColumnCell;
   
    isNaNMask = isnan(numericData);
    if all(isNaNMask(:))
        numericData = [];
        return;
    end
    
    % Trim all-NaN leading and trailing rows and columns from numeric array
    [rowStart, colStart] = getCorner(isNaNMask, 'first');
    [rowEnd, colEnd] = getCorner(isNaNMask, 'last');
    numericData = numericData(rowStart:rowEnd, colStart:colEnd);
    
    % Trim the original data
    data = data(rowStart:rowEnd, colStart:colEnd);
    
    % Determine how much of the trimmed data is logical and restore logical
    % type if all values were logical.
    if all(cellfun('islogical', data))
        numericData = logical(numericData);
    end
    
    % Ensure numericArray is 0x0 empty.
    if isempty(numericData)
        numericData = [];
    end
end

%--------------------------------------------------------------------------
function [row, col] = getCorner(mask, firstlast)
    isLast = strcmp(firstlast,'last');

    % Find first (or last) row that is not all true in the mask.
    row = find(~all(mask,2), 1, firstlast);
    if isempty(row)
        row = emptyCase(isLast, size(mask,1));
    end
    
    % Find first (or last) column that is not all true in the mask.
    col = find(~all(mask,1), 1, firstlast);
    % Find returns empty if there are no rows/columns that contain a false value.
    if isempty(col)
        col = emptyCase(isLast, size(mask,2));
    end    
end

%--------------------------------------------------------------------------
function dim = emptyCase(isLast, dimSize)
    if isLast
        dim = dimSize;
    else
        dim = 1;
    end
end


