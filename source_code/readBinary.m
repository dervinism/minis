 function data = readBinary(fileName, format)
% The function reads short binary files.
%
% Input: fileName is the name of the binary file to read.
%        format is the data format (i.e., 'int16', 'double', etc).
%        Currently supported types are 'int16' and 'double'.
%
% Output: data is the data vector read from fileName.

chunkSize = 1000000;

fid = []; data = [];

d = dir(fileName);
if strcmpi(format, 'int16')
    nSampsTotal = d.bytes/2;
elseif strcmpi(format, 'double')
    nSampsTotal = d.bytes/8;
end
nChunksTotal = ceil(nSampsTotal/chunkSize);

try
  fid = fopen(fileName, 'r');
  
  chunkInd = 1;
  while 1
    fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
    dat = fread(fid, [1 chunkSize], ['*' format]);
    if ~isempty(dat)
      data = [data dat]; %#ok<*AGROW>
    else
      break
    end
    chunkInd = chunkInd+1;
  end
  
  fclose(fid);
  
catch me
  if ~isempty(fid)
    fclose(fid);
  end
  
  rethrow(me)
  
end
