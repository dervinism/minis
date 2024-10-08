function [SDmean, SDmin, SDmax] = estimateSDs(sweepData, excludedIndices, sweepSize)



nFiles = floor(size(sweepData,1)/sweepSize);
fileData = reshape(sweepData(1:sweepSize*nFiles,:)', sweepSize*size(sweepData,2), nFiles)';
fileExcludedIndices = reshape(excludedIndices(1:sweepSize*nFiles,:)', sweepSize*size(excludedIndices,2), nFiles)';
if nFiles < size(sweepData,1)/sweepSize
    fileData = [fileData; reshape(sweepData(end-sweepSize+1:end,:)', 1, sweepSize*size(sweepData,2))];
    fileExcludedIndices = [fileExcludedIndices; reshape(excludedIndices(end-sweepSize+1:end,:)', 1, sweepSize*size(excludedIndices,2))];
    nFiles = nFiles + 1;
end

SD = zeros(1,nFiles);
for iFile = 1:nFiles
    SD(iFile) = std(fileData(iFile,~fileExcludedIndices(iFile,:)),1);
end
SDmean = mean(SD);
if isempty(SDmean)
    SDmean = NaN;
end
SDmin = min(SD);
if isempty(SDmin)
    SDmin = NaN;
end
SDmax = max(SD);
if isempty(SDmax)
    SDmax = NaN;
end