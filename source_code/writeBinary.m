function outputFilename = writeBinary(data, fileName, format)
% The function writes a short binary file.
%
% Input: data is the data vector to write.
%        fileName is the name of the binary file when saved.
%        format is optional. If not supplied, binary file will be saved in
%           the same format as the data.
%
% Output: outputFilename is the name of the binary file when saved.

chunkSize = 1000000;

fidOut = [];

if nargin == 3
    if strcmp(format, 'int8')
        data = int8(data);
    elseif strcmp(format, 'uint8')
        data = uint8(data);
    elseif strcmp(format, 'int16')
        data = int16(data);
    elseif strcmp(format, 'uint16')
        data = uint16(data);
    elseif strcmp(format, 'int32')
        data = int32(data);
    elseif strcmp(format, 'uint32')
        data = uint32(data);
    elseif strcmp(format, 'single')
        data = single(data);
    elseif strcmp(format, 'int64')
        data = int64(data);
    elseif strcmp(format, 'uint64')
        data = uint64(data);
    elseif strcmp(format, 'double')
        data = double(data);
    else
        error(['Unsupported data format. '...
            'Supported formats are double, single, int8, int16, int32, int64, uint8, uint16, uint32, uint64']);
    end
end

w = whos('data');
if nargin < 3
    format = w.class;
end
if strcmp(format, 'int8') || strcmp(format, 'uint8')
    nSampsTotal = w.bytes;
elseif strcmp(format, 'int16') || strcmp(format, 'uint16')
    nSampsTotal = w.bytes/2;
elseif strcmp(format, 'int32') || strcmp(format, 'uint32') || strcmp(format, 'single')
    nSampsTotal = w.bytes/4;
elseif strcmp(format, 'int64') || strcmp(format, 'uint64') || strcmp(format, 'double')
    nSampsTotal = w.bytes/8;
else
    error(['Unsupported data format. '...
        'Supported formats are double, single, int8, int16, int32, int64, uint8, uint16, uint32, uint64']);
end

nChunksTotal = ceil(nSampsTotal/chunkSize);

try
    outputFilename  = fileName;
    fidOut = fopen(outputFilename, 'w');
    
    chunkInd = 1;
    while 1
        fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
        inds = (1:chunkSize) + (chunkInd-1)*chunkSize;
        if inds(1) > numel(data)
            break
        elseif inds(end) > numel(data)
            inds = inds(1):numel(data);
        end
        dat = data(inds);
        fwrite(fidOut, dat, format);
        chunkInd = chunkInd+1;
    end
    
    fclose(fidOut);
    
catch me
    if ~isempty(fidOut)
        fclose(fidOut);
    end
    
    rethrow(me)
    
end