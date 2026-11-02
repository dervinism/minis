function writeEventsTSV(filePath, headers, dataMatrix)
% WRITEEVENTSTSV writes a header row followed by numeric data rows to a
% tab-separated (.tsv) file that opens directly in Excel, with no manual
% copy-paste required.
%
%   WRITEEVENTSTSV(filePath, headers, dataMatrix) writes HEADERS (a cell
%   array of strings, one per column) as the first row, then one row per
%   row of DATAMATRIX, all separated by tabs. Numeric values are written
%   with 14 significant digits.

if ~isempty(dataMatrix) && numel(headers) ~= size(dataMatrix, 2)
    error('writeEventsTSV:columnMismatch', ...
        'Number of headers (%d) does not match number of data columns (%d).', ...
        numel(headers), size(dataMatrix, 2));
end

fid = fopen(filePath, 'wt+');
if fid == -1
    error('writeEventsTSV:fopen', 'Could not open file for writing: %s', filePath);
end

nCols = numel(headers);
headerFormat = [repmat('%s\t', 1, nCols-1) '%s\n'];
fprintf(fid, headerFormat, headers{:});

if ~isempty(dataMatrix)
    dataFormat = [repmat('%.14g\t', 1, size(dataMatrix,2)-1) '%.14g\n'];
    fprintf(fid, dataFormat, dataMatrix');
end

fclose(fid);
end
